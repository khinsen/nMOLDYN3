import operator
from os import system
from shutil import copyfile
from string import *
from tempfile import mktemp
from time import asctime,localtime,time
                                                 
from Scientific import N
from Scientific.IO.NetCDF import NetCDFFile
from Scientific.IO.TextFile import TextFile
from Scientific.Functions.Interpolation import InterpolatingFunction
from Scientific.Threading.TaskManager import TaskManager

from MMTK import *
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput

from nMOLDYN.Tests.nMOLDYN_Reference.misc import timePrepare, BincohList, OnesList, BcohList,\
                                                 logfileInit, logfileUpdate, cleanUp, getReference, \
                                                 uniqueGroupPairs,  getVariables
                         
from nMOLDYN.Tests.nMOLDYN_Reference.calc import AutoCorrelationFunction, getMeanSquareDisplacement,\
                                                 GaussianWindow, getAngularVelocity, acf, diff,\
                                                 sphericalHarmonics, CorrelationFunction, rotMatrix, \
                                                 WignerFunctions, histogram, qMatrix, \
                                                 fft, inverse_fft, AutoCorrelationFunction_2

try:
    import threading
except ImportError:
    threading = None

threading = None

def MeanSquareDisplacement(traj,atoms,timeInfo,weightList=None,normVect=None):

    if weightList is None:
        weightList = OnesList(traj.universe,atoms.atomList())
    tim = timePrepare(traj,timeInfo)
    MSD = N.zeros((len(tim)),N.Float)

    file, filename = logfileInit('msd')

    if threading is None:
        nthreads = 0
    else:
        nthreads = 1

    if nthreads > 0:
        tasks = TaskManager(nthreads)

    def oneatom(series, weight, MSD, lock=None):
        dsq      = series[:,0]**2+series[:,1]**2+series[:,2]**2
        sum_dsq1 = N.add.accumulate(dsq)
        sum_dsq2 = N.add.accumulate(dsq[::-1])
        sumsq    = 2.*sum_dsq1[-1]
        msd      = sumsq - N.concatenate(([0.], sum_dsq1[:-1])) \
                   - N.concatenate(([0.], sum_dsq2[:-1]))
        Sab      = 2.*N.add.reduce(acf(series, 0, lock), 1)
        msd_term = (msd-Sab)/(len(series)-N.arange(len(series)))*weight
        if lock is not None: lock.acquire()
        N.add(MSD, msd_term, MSD)
        if lock is not None: lock.release()

    i = 0
    normFact = len(atoms.atomList())
    for at in atoms.atomList():
        series   = traj.readParticleTrajectory(at,first=timeInfo[0], \
                        last=timeInfo[1],skip=timeInfo[2]).array
        if normVect is not None: series = series * normVect[N.NewAxis,:]
        if nthreads == 0:
            oneatom(series, weightList[at], MSD)
        else:
            tasks.runTask(oneatom, (series, weightList[at], MSD))
        i = logfileUpdate(file, i, normFact)

    if nthreads > 0:
        tasks.terminate()

    cleanUp(file,filename)
    return tim[:len(MSD)], MSD


def CoherentScatteringFunction(traj,qVect,ncDataFN,
                               atoms=None,bcoh=None,timeInfo=None,nsteps=None):

    file, filename = logfileInit('csf')

    if atoms is None: atoms = traj.universe
    mask = atoms.booleanMask().array
    if bcoh is None:
        bcoh = BcohList(traj.universe,atoms)
    tim  = timePrepare(traj,timeInfo)
    if nsteps is None:
        nsteps = len(tim)
    nsteps = min(nsteps, len(tim))
    Fcoh = N.zeros((nsteps,len(qVect[1])),N.Float)
    bcoh = N.repeat(bcoh.array, mask)[:, N.NewAxis]

    comment = ''
    try:
        comment = comment + 'Trajectory: ' + traj.filename + '\n'
    except AttributeError:
        pass
    
    file2       = NetCDFFile(ncDataFN, 'w', comment)
    file2.title = 'Coherent Scattering Function'
    file2.createDimension('Q', None)
    file2.createDimension('TIME', nsteps)
    SF      = file2.createVariable('Fqt-total', N.Float, ('Q','TIME'))
    Time    = file2.createVariable('time', N.Float, ('TIME',))
    Qlength = file2.createVariable('q', N.Float, ('Q',))

    Time[:] = tim[:nsteps]    
 
    for j in range(len(qVect[1])):
        qarray = N.array(map(lambda v: v.array, qVect[1][j]))
        qarray = traj.universe._boxToRealPointArray(qarray)
        qVect[1][j] = N.transpose(qarray)

    normFact = len(qVect[0])
    for j in range(len(qVect[0])):
        rho = N.zeros((len(tim), qVect[1][j].shape[1]), N.Complex)
        for i in range(len(tim)):
            conf = traj.configuration[timeInfo[0]+i*timeInfo[2]]
            conf.convertToBoxCoordinates()
            sekw     = N.repeat(conf.array, mask)
            rho[i,:] = N.add.reduce(bcoh*N.exp(1j*N.dot(sekw, qVect[1][j])))
        rho        = AutoCorrelationFunction(rho)[:nsteps]
        Fcoh[:,j]  = N.add.reduce(rho,1)/qVect[1][j].shape[1]
        SF[j,:]    = Fcoh[:,j]
        Qlength[j] = qVect[0][j]
        file2.flush()
        gj = logfileUpdate(file,j,normFact)
        
    cleanUp(file,filename)
    file2.close()
    return 


def CoherentScatteringAR(trajectory, sf_filename, memory_filename,
                         dsf_filename, order, qvect, atoms=None, weights=None,
                         timeInfo=None, nsteps=None, nfreq=None,
                         precision=None):

    from Scientific.Signals.Models import AutoRegressiveModel, \
                                          AveragedAutoRegressiveModel

    logfile, logfilename = logfileInit('ar_csf')
    if atoms is None: atoms = trajectory.universe
    mask = atoms.booleanMask().array
    if weights is None:
        weights = BcohList(trajectory.universe,atoms)
    weights = N.repeat(weights.array, mask)[:, N.NewAxis]
    tim  = timePrepare(trajectory,timeInfo)
    if nsteps is None:
        nsteps = len(tim)
    nsteps = min(nsteps, len(tim))
    if nfreq is None:
        nfreq = nsteps
    dt = tim[1]-tim[0]

    comment = ''
    try:
        comment = comment + 'Trajectory: ' + trajectory.filename + '\n'
    except AttributeError:
        pass
        
    sf_file = NetCDFFile(sf_filename, 'w', comment)
    sf_file.title = 'Coherent Scattering Function (AR order %d)' % order
    sf_file.createDimension('Q', None)
    sf_file.createDimension('TIME', nsteps)
    sf_file.createDimension('TIME_MEMORY', order+order/2)
    sf_file.createDimension('MODEL', order+2)
    sf_file.createDimension('FREQUENCY', nsteps+1)
    
    Time    = sf_file.createVariable('time', N.Float, ('TIME',))
    Time[:] = tim[:nsteps]
    
    Time    = sf_file.createVariable('time_memory', N.Float, ('TIME_MEMORY',))
    Time[:] = dt*N.arange(order+order/2)
    
    SF      = sf_file.createVariable('Fqt', N.Float, ('Q','TIME'))
    Qlength_SF = sf_file.createVariable('q', N.Float, ('Q',))
    
    model_r1 = sf_file.createVariable('ar_coefficients_real', N.Float, ('Q', 'MODEL'))
    model_i1 = sf_file.createVariable('ar_coefficients_imag', N.Float, ('Q', 'MODEL'))
       
    SFMEM = sf_file.createVariable('Fqt_memory', N.Float, ('Q','TIME_MEMORY'))
         
    freq_max = 1./(2*dt)
    dfreq = freq_max/nsteps
    freq = N.arange(nsteps+1)*dfreq
    Freq = sf_file.createVariable('frequency', N.Float, ('FREQUENCY',))
    Freq[:] = freq
    
    DSF = sf_file.createVariable('Sqw', N.Float, ('Q','FREQUENCY'))

    for j in range(len(qvect[1])):
        qarray = N.array(map(lambda v: v.array, qvect[1][j]))
        qarray = trajectory.universe._boxToRealPointArray(qarray)
        qvect[1][j] = N.transpose(qarray)

    normFact = len(qvect[0])
    for j in range(len(qvect[0])):
        rho = N.zeros((len(tim), qvect[1][j].shape[1]), N.Complex)
        for i in range(len(tim)):
            conf = trajectory.configuration[timeInfo[0]+i*timeInfo[2]]
            conf.convertToBoxCoordinates()
            sekw = N.repeat(conf.array, mask)
            rho[i,:] = N.add.reduce(weights*N.exp(1j*N.dot(sekw,
                                                             qvect[1][j])))
        model = AveragedAutoRegressiveModel(order, dt)
        for i in range(rho.shape[1]):
            data = rho[:, i]
            data = data - N.add.reduce(data)/len(data)
            m = AutoRegressiveModel(order, data, dt)
            model.add(m)
        parameters = N.concatenate((model.coeff[::-1],
                                  N.array([model.sigma, model.variance])))
        average = N.add.reduce(N.add.reduce(rho))/N.multiply.reduce(rho.shape)
        
        fcoh = model.correlation(nsteps).values.real + (average*N.conjugate(average)).real
        SF[j,:] = fcoh
        Qlength_SF[j] = qvect[0][j]
        model_r1[j, :] = parameters.real
        model_i1[j, :] = parameters.imag
        
        spectrum = model.spectrum(2.*N.pi*freq)
        DSF[j,:] = spectrum.values
        
        if precision is not None:
            from ARMP import MPAutoRegressiveModel
            model = MPAutoRegressiveModel(model, precision)
        mem = model.memoryFunction(order+order/2).values.real
        SFMEM[j,:] = mem
        
        gj = logfileUpdate(logfile,j,normFact)
        
    cleanUp(logfile, logfilename)
    
    sf_file.close()

    return 


def CurrentDensityAutocorrelationFunction(traj,qVect,ncDataFN,
                                          atoms=None,bcoh=None,timeInfo=None,
                                          nsteps=None):

    file, filename = logfileInit('cdaf')

    if atoms is None: atoms = traj.universe
    mask = atoms.booleanMask().array
    if bcoh is None:
        bcoh = BcohList(traj.universe,atoms)
    tim  = timePrepare(traj,timeInfo)
    if nsteps is None:
        nsteps = len(tim)
    nsteps = min(nsteps, len(tim))
    bcoh = N.repeat(bcoh.array, mask)[:, N.NewAxis]

    comment = ''
    try:
        comment = comment + 'Trajectory: ' + traj.filename + '\n'
    except AttributeError:
        pass
    
    file2       = NetCDFFile(ncDataFN, 'w', comment)
    file2.title = 'Current Density Autocorrelation Function'
    file2.createDimension('Q', None)
    file2.createDimension('TIME', nsteps)
    cdaf    = file2.createVariable('cdaf', N.Float, ('Q','TIME'))
    Time    = file2.createVariable('time', N.Float, ('TIME',))
    Qlength = file2.createVariable('q', N.Float, ('Q',))

    Time[:] = tim[:nsteps]
 
    for j in range(len(qVect[1])):
        qVect[1][j] = N.transpose(array(map(lambda v: v.array, qVect[1][j])))

    normFact = len(qVect[0])
    for j in range(len(qVect[0])):
        sum = N.zeros((len(tim),qVect[1][j].shape[1]),'D')
        for i in range(timeInfo[0],timeInfo[1],timeInfo[2]):
            sekw     = N.repeat(traj.configuration[i].array, mask)
            vel      = N.repeat(traj.velocities[i].array, mask)
            vq       = N.dot(vel, qVect[1][j])/qVect[0][j]
            sum[i,:] = N.add.reduce(bcoh*vq*N.exp(1j*N.dot(sekw, qVect[1][j])))
        sum        = AutoCorrelationFunction(sum)[:nsteps]
        cdaf[j,:]  = N.add.reduce(sum,1)/qVect[1][j].shape[1]
        Qlength[j] = qVect[0][j]
        file2.flush()
        gj = logfileUpdate(file,j,normFact)
        
    cleanUp(file,filename)
    file2.close()
    return 


def IncoherentScatteringFunction(traj,qVect,ncDataFN,
                                 atoms=None,bincoh=None,timeInfo=None,
                                 nsteps=None):

    file, filename = logfileInit('isf')

    if atoms is None: atoms = traj.universe
    if bincoh is None: bincoh = BincohList(traj.universe,atoms)
    tim    = timePrepare(traj,timeInfo)
    if nsteps is None:
        nsteps = len(tim)
    nsteps = min(nsteps, len(tim))
    Fincoh = N.zeros((nsteps,len(qVect[1])),N.Float)
 
    comment = ''
    try:
        comment = comment + 'Trajectory: ' + traj.filename + '\n'
    except AttributeError:
        pass
    
    file2       = NetCDFFile(ncDataFN, 'w', comment)
    file2.title = 'Incoherent Scattering Function'
    file2.createDimension('Q', None)
    file2.createDimension('TIME', nsteps)
    SF      = file2.createVariable('Fqt-total', N.Float, ('Q','TIME'))
    Time    = file2.createVariable('time', N.Float, ('TIME',))
    Qlength = file2.createVariable('q', N.Float, ('Q',))
    Time[:]=tim[:nsteps]
 
    for j in range(len(qVect[1])):
        qVect[1][j] = N.transpose(N.array(map(lambda v: v.array, qVect[1][j])))

    normFact = len(qVect[0])
    for j in range(len(qVect[0])):
        for at in atoms.atomList():
            sekw        = traj.readParticleTrajectory(at,first=timeInfo[0],\
                               last=timeInfo[1],skip=timeInfo[2]).array
            series      = N.exp(1j*N.dot(sekw, qVect[1][j]))
            res         = AutoCorrelationFunction(series)[:nsteps]
            Fincoh[:,j] = Fincoh[:,j] + N.add.reduce(bincoh[at]*res,1)/\
                                        qVect[1][j].shape[1]
        SF[j,:]    = Fincoh[:,j]
        Qlength[j] = qVect[0][j]
        file2.flush()
        gj = logfileUpdate(file,j,normFact)

    cleanUp(file,filename)
    file2.close()
    return 


def IncoherentScatteringFunctionGaussian(traj,qVect,ncDataFN,
                                         atoms=None,bincoh=None,timeInfo=None,
                                         nsteps=None):

    file, filename = logfileInit('isfg')

    if atoms is None: atoms = traj.universe
    if bincoh is None: bincoh = BincohList(traj.universe,atoms)
    tim    = timePrepare(traj,timeInfo)
    qun    = N.array(map(lambda v: v*v, qVect[0]))
    if nsteps is None:
        nsteps = len(tim)
    nsteps = min(nsteps, len(tim))
    Fincoh = N.zeros((nsteps,len(qVect[1])),N.Float)

    comment = ''
    try:
        comment = comment + 'Trajectory: ' + traj.filename + '\n'
    except AttributeError:
        pass
    
    file2       = NetCDFFile(ncDataFN, 'w', comment)
    file2.title = 'Incoherent Scattering Function (Gaussian Approximation)'
    file2.createDimension('Q', None)
    file2.createDimension('TIME', nsteps)
    SF      = file2.createVariable('Fqt-total', N.Float, ('Q','TIME'))
    Time    = file2.createVariable('time', N.Float, ('TIME',))
    Qlength = file2.createVariable('q', N.Float, ('Q',))
    Time[:]=tim[:nsteps]

    i = 0
    normFact = len(atoms.atomList())
    for at in atoms.atomList():
        sekw = traj.readParticleTrajectory(at,first=timeInfo[0],\
                    last=timeInfo[1],skip=timeInfo[2]).array
        res  = getMeanSquareDisplacement(sekw)[:nsteps]
        for qi in range(len(qun)):
            Fincoh[:, qi] = Fincoh[:, qi] + N.exp(-res*qun[qi]/6.)*bincoh[at]
        i = logfileUpdate(file,i,normFact)
    for i in range(len(qun)):
        SF[i, :] = Fincoh[:, i]
        Qlength[i] = N.sqrt(qun[i])

    cleanUp(file,filename)
    file2.close()
    return 


def IncoherentScatteringAR(trajectory, sf_filename, memory_filename,
                           dsf_filename, order, qvect,
                           atoms=None, weights=None, timeInfo=None,
                           nsteps=None, nfreq=None, precision=None):

    from Scientific.Signals.Models import AutoRegressiveModel, \
                                          AveragedAutoRegressiveModel

    logfile, logfilename = logfileInit('ar_isf')
    if atoms is None: atoms = trajectory.universe
    if weights is None: weights = BincohList(trajectory.universe, atoms)
    tim    = timePrepare(trajectory,timeInfo)
    if nsteps is None:
        nsteps = len(tim)
    nsteps = min(nsteps, len(tim))
    dt = tim[1]-tim[0]
 
    comment = ''
    try:
        comment = comment + 'Trajectory: ' + trajectory.filename + '\n'
    except AttributeError:
        pass
    
    sf_file = NetCDFFile(sf_filename, 'w', comment)
    sf_file.title = 'Incoherent Scattering Function (AR order %d)' % order
    sf_file.createDimension('Q', len(qvect[0]))
    sf_file.createDimension('TIME', nsteps)
    sf_file.createDimension('MODEL', order+2)

    sf_file.createDimension('FREQUENCY', nsteps + 1)
    sf_file.createDimension('TIME_MEMORY', order+order/2)
    
    Time    = sf_file.createVariable('time', N.Float, ('TIME',))
    Time[:] = tim[:nsteps]

    Time    = sf_file.createVariable('time_memory', N.Float, ('TIME_MEMORY',))
    Time[:] = dt*N.arange(order+order/2)
    
    SF      = sf_file.createVariable('Fqt', N.Float, ('Q','TIME'))
    Qlength_SF = sf_file.createVariable('q', N.Float, ('Q',))

    freq_max = 1./(2*dt)
    dfreq = freq_max/nsteps
    freq = N.arange(nsteps+1)*dfreq
    Freq = sf_file.createVariable('frequency', N.Float, ('FREQUENCY',))
    Freq[:] = freq
    
    DSF = sf_file.createVariable('Sqw', N.Float, ('Q','FREQUENCY'))

    SFMEM = sf_file.createVariable('Fqt_memory', N.Float, ('Q','TIME_MEMORY'))
    
    model_r = sf_file.createVariable('ar_coefficients_real', N.Float, ('Q', 'MODEL'))
    model_i = sf_file.createVariable('ar_coefficients_imag', N.Float, ('Q', 'MODEL'))
 
    for j in range(len(qvect[1])):
        qarray = N.array(map(lambda v: v.array, qvect[1][j]))
        qvect[1][j] = N.transpose(qarray)

    natoms = atoms.numberOfAtoms()
    normFact = len(qvect[0])*natoms
    for j in range(len(qvect[0])):

        Qlength_SF[j] = qvect[0][j]
        SF[j,:] = 0.
        
        DSF[j,:] = 0.

        SFMEM[j,:] = 0.

        tot_model = AveragedAutoRegressiveModel(order, dt)
        cf = 0.
        ds = 0.
        atom_count = 0
        for atom in atoms.atomList():
            if weights[atom] != 0.:
                at = trajectory.readParticleTrajectory(atom,first=timeInfo[0],\
                                                       last=timeInfo[1],
                                                       skip=timeInfo[2]).array
                rho = N.exp(1j*N.dot(at, qvect[1][j]))
                rho_av = N.add.reduce(rho)/rho.shape[0]
                model = AveragedAutoRegressiveModel(order, dt)
                for i in range(rho.shape[1]):
                    data = rho[:, i] - rho_av[i]
                    m = AutoRegressiveModel(order, data, dt)
                    model.add(m, weights[atom])
                    tot_model.add(m, weights[atom])
                rho_av_sq = (rho_av*N.conjugate(rho_av)).real
                average = N.add.reduce(rho_av_sq)/rho.shape[1]
                
                cf = cf + weights[atom]* (model.correlation(nsteps).values.real + average)
                SF[j,:] = cf

                ds = ds + weights[atom]*(model.spectrum(2.*N.pi*freq).values)
                DSF[j,:] = ds

            atom_count = atom_count + 1
            gj = logfileUpdate(logfile,j*natoms+atom_count,normFact)

        parameters = N.concatenate((tot_model.coeff[::-1],
                                  N.array([tot_model.sigma,
                                         tot_model.variance])))

        if precision is not None:
            from ARMP import MPAutoRegressiveModel
            tot_model = MPAutoRegressiveModel(tot_model, precision)
        mem = tot_model.memoryFunction(order+order/2).values.real
        SFMEM[j,:] = mem

        model_r[j, :] = parameters.real
        model_i[j, :] = parameters.imag

    cleanUp(logfile, logfilename)
    sf_file.close()

    return 


def ElasticIncoherentStructureFactor(traj,qVect,
           atoms=None,bincoh=None,timeInfo=None):

    file, filename = logfileInit('eisf')

    tim  = timePrepare(traj,timeInfo)
    EISF = N.zeros((len(qVect[1]),),N.Float)
    if atoms is None: atoms = traj.universe
    if bincoh == None: bincoh = BincohList(traj.universe,atoms)
    if timeInfo is None: timeInfo = (0,len(traj)-1,1)
    for j in range(len(qVect[1])):
        qVect[1][j] = N.transpose(N.array(map(lambda v: v.array, qVect[1][j])))

    i = 0
    normFact = len(atoms.atomList())
    for at in atoms.atomList():
        tr = traj.readParticleTrajectory(at,first=timeInfo[0],\
                  last=timeInfo[1],skip=timeInfo[2]).array
        for j in range(len(qVect[0])):
            a       = N.add.reduce(N.exp(1j*N.dot(tr,qVect[1][j])))/len(tim)
            EISF[j] = EISF[j] +  N.add.reduce(bincoh[at]*(a*N.conjugate(a)).real)/\
                                                        qVect[1][j].shape[1]
        i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    return qVect[0], EISF


def VelocityAutocorrelationFunction_vel(traj,atoms=None,
            normVect=None,weightList=None,timeInfo=None):

    file, filename = logfileInit('vacf_vel')

    if atoms is None: atoms = traj.universe
    if weightList is None:
        weightList = OnesList(traj.universe,atoms.atomList())
    tim  = timePrepare(traj,timeInfo)
    vacf = N.zeros((len(tim)),N.Float)

    i = 0
    normFact = len(atoms.atomList())
    for at in atoms.atomList():
        tr= traj.readParticleTrajectory(at,first=timeInfo[0],\
                 last=timeInfo[1],skip=timeInfo[2],variable='velocities').array
        if normVect is None:
            res = N.add.reduce(AutoCorrelationFunction(tr), -1)/3.
        else:
            res = AutoCorrelationFunction(N.dot(tr, normVect))
        #res = res/res[0]
        vacf = vacf+res*weightList[at]
        i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    return tim, vacf


def VelocityAutocorrelationFunction_xyz(traj,atoms=None, 
            normVect=None,weightList=None,timeInfo=None,diffScheme='fast'):

    file, filename = logfileInit('vacf_xyz')

    if atoms is None: atoms = traj.universe
    if weightList is None:
        weightList = OnesList(traj.universe,atoms.atomList())        
    tim  = timePrepare(traj,timeInfo)
    dt   = tim[1]-tim[0]
    vacf = N.zeros(len(tim),N.Float)

    i = 0
    normFact = len(atoms.atomList())
    velo = N.zeros((len(tim),3),'d')
    for at in atoms.atomList():
        tr = traj.readParticleTrajectory(at,first=timeInfo[0],\
                  last=timeInfo[1],skip=timeInfo[2]).array
        for ix in range(3): velo[:,ix] = diff(tr[:,ix],diffScheme,dt)
        if normVect is None:
            res = N.add.reduce(AutoCorrelationFunction(velo), -1)/3.
        else:
            res = AutoCorrelationFunction(N.dot(velo, normVect))
        #res  = res/res[0]
        vacf = vacf+res*weightList[at]
        i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    return tim, vacf


def DensityOfStates_vel(traj,alpha=1.,atoms=None,
           normVect=None,weightList=None,timeInfo=None):    

    file, filename = logfileInit('dos_vel')
    if atoms is None: atoms = traj.universe    
    if weightList is None:
        weightList = OnesList(traj.universe,atoms.atomList())    
    tim      = timePrepare(traj,timeInfo)
    dos      = N.zeros(len(tim),N.Float)
    TimeStep = tim[1]-tim[0]
 
    i = 0
    normFact = len(atoms.atomList())
    for at in atoms.atomList():
        velo = traj.readParticleTrajectory(at,first=timeInfo[0],\
               last=timeInfo[1],skip=timeInfo[2],variable='velocities').array
        if normVect is None:
            res = N.add.reduce(AutoCorrelationFunction(velo), -1)/3.
        else: res = AutoCorrelationFunction(N.dot(velo, normVect))
        
        dos = dos + fft(GaussianWindow(res,alpha)).real[:len(tim)] * \
              weightList[at]
        i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    frequencies = N.arange(len(tim))/(2.*len(tim)*TimeStep)
    return frequencies, 0.5*TimeStep*dos


def DensityOfStates_xyz(traj,alpha=1.,atoms=None,
           normVect=None,weightList=None,timeInfo=None,diffScheme='fast'):

    file, filename = logfileInit('dos_xyz')

    if atoms is None: atoms = traj.universe      
    if weightList is None:
        weightList = OnesList(traj.universe,atoms.atomList()) 
    tim      = timePrepare(traj,timeInfo)
    dos      = N.zeros(len(tim),N.Float) # -2
    TimeStep = tim[1]-tim[0]

    i = 0
    normFact = len(atoms.atomList())
    velo = N.zeros((len(tim),3),'d')
    for at in atoms.atomList():
        tr= traj.readParticleTrajectory(at,first=timeInfo[0],\
                 last=timeInfo[1],skip=timeInfo[2]).array
        for ix in range(3): velo[:,ix] = diff(tr[:,ix],diffScheme,TimeStep)
        if normVect is None:
            res = N.add.reduce(AutoCorrelationFunction(velo), -1)/3.
        else:
            res = AutoCorrelationFunction(N.dot(velo, normVect))
        dos = dos+fft(GaussianWindow(res,alpha)).real[:len(tim)]*weightList[at]
        i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    frequencies = N.arange(len(tim))/(2.*len(tim)*TimeStep)
                                #-2
    return frequencies, 0.5*TimeStep*dos


def AngularVelocityAutocorrelationFunction(traj,groups,refGroup,normVect=None,
                                           timeInfo=None,diffScheme='fast'):

    file, filename = logfileInit('avacf')

    tim   = timePrepare(traj,timeInfo)
    avacf = N.zeros(len(tim),N.Float)

    ngr = 0
    for key in groups.keys(): ngr = ngr + len(groups[key])

    print 'Group names: ',groups.keys()
    print 'Total number of groups is: ',ngr
    
    i = 0
    normFact = ngr
    for key in groups.keys():
        for  gr in groups[key]:
            try:
                if refGroup[key]:
                    group, refs = getReference(gr, refGroup[key])
                else:
                    group = gr
                    refs = None
            except:
                print 'Warning: Exception raised for',gr,i,key
                refs = None
                group = gr
            angvel = getAngularVelocity(traj,group,timeInfo,reference=refs,
                                        diffScheme=diffScheme)
            if normVect is None:
                res = N.add.reduce(AutoCorrelationFunction(angvel), -1)/3.
            else:
                res = AutoCorrelationFunction(N.dot(angvel,normVect))
            avacf = avacf+res/res[0]
            i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    return tim, avacf/ngr


def MethylAngularVelocityAutocorrelationFunction(traj,atoms,groups,refGroup,normVect=None,
                                           timeInfo=None,diffScheme='fast'):

    file, filename = logfileInit('methyl_avacf')

    tim   = timePrepare(traj,timeInfo)
    avacf = N.zeros(len(tim),N.Float)

    print atoms
       
    for key in groups.keys():
        remove=[]
        for  gr in groups[key]:
             if gr[0] not in atoms: 
               remove.append(gr)
        for r in remove:
            groups[key].remove(r)
      
       
    ngr = 0
    for key in groups.keys(): ngr = ngr + len(groups[key])

    print 'Group names: ',groups.keys()
    print 'Total number of groups is: ',ngr
    
    i = 0
    normFact = ngr
    for key in groups.keys():
        for  gr in groups[key]:
            try:
                if refGroup[key]:
                    group, refs = getReference(gr, refGroup[key])
                else:
                    group = gr
                    refs = None
            except:
                print 'Warning: Exception raised for',gr,i,key
                refs = None
                group = gr
            c_coordinate = traj.readParticleTrajectory(gr[0],
               first=timeInfo[0],last=timeInfo[1],
               skip=timeInfo[2])
            h_coordinate = traj.readParticleTrajectory(gr[1],
               first=timeInfo[0],last=timeInfo[1],
               skip=timeInfo[2])
            
            mu = N.zeros([len(c_coordinate),3], N.Float)
            for index in range(len(c_coordinate)):
               mu[index] = (c_coordinate[index]-h_coordinate[index]).normal()
            
            res = AutoCorrelationFunction_2(mu)
            avacf = avacf+res/res[0]
            i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    return tim, 1.5*avacf/ngr-0.5


def MethylAlignedAngularVelocityAutocorrelationFunction(traj,atoms,groups,refGroup,normVect=None,
                                           timeInfo=None,diffScheme='fast'):

    file, filename = logfileInit('methyl_aligned_avacf')

    tim   = timePrepare(traj,timeInfo)
    avacf = N.zeros(len(tim),N.Float)

    print atoms
       
    for key in groups.keys():
        remove=[]
        for  gr in groups[key]:
             if gr[0] not in atoms: 
               remove.append(gr)
        for r in remove:
            groups[key].remove(r)
      
       
    ngr = 0
    for key in groups.keys(): ngr = ngr + len(groups[key])

    print 'Group names: ',groups.keys()
    print 'Total number of groups is: ',ngr
    
    i = 0
    normFact = ngr
    for key in groups.keys():
        for  gr in groups[key]:
            try:
                if refGroup[key]:
                    group, refs = getReference(gr, refGroup[key])
                else:
                    group = gr
                    refs = None
            except:
                print 'Warning: Exception raised for',gr,i,key
                refs = None
                group = gr
            c_coordinate = traj.readParticleTrajectory(gr[0],
               first=timeInfo[0],last=timeInfo[1],
               skip=timeInfo[2])
            h_coordinate = traj.readParticleTrajectory(gr[1],
               first=timeInfo[0],last=timeInfo[1],
               skip=timeInfo[2])
            c_parent = gr[0].bondedTo()[0]  
            c_parent_coordinate = traj.readParticleTrajectory(c_parent,
               first=timeInfo[0],last=timeInfo[1],
               skip=timeInfo[2])
            
            mu = N.zeros([len(c_coordinate),3], N.Float)
            mu_axis = N.zeros([len(c_coordinate),3], N.Float)
            for index in range(len(c_coordinate)):
               mu[index] = (h_coordinate[index]-c_coordinate[index]).normal()
               mu_axis[index] = (c_coordinate[index]-c_parent_coordinate[index]).normal()
               if(index>0):
                  a_b = N.dot(mu_axis[index],mu_axis[0])
                  a_c = N.dot(mu_axis[index],mu[index])
                  b_c = N.dot(mu_axis[0],mu[index])
                  denom = a_b**2-1
                  project_1 = (a_b*b_c-a_c)/denom
                  project_2 = (a_b*a_c-b_c)/denom
                  project_perp = mu[index] - project_1*mu_axis[index] - project_2*mu_axis[0]
                  
                  m_perp = N.dot(mu_axis[index], mu_axis[0])*mu_axis[0]
                  m_par = mu_axis[index] - m_perp
                  axis_new = m_perp - m_par
                  
                  mu[index] = project_1*mu_axis[0] + project_2*axis_new + project_perp
                  
                  
            
            res = AutoCorrelationFunction_2(mu)
            avacf = avacf+res/res[0]
            i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    return tim, 1.5*avacf/ngr-0.5


def MethylAxisAngularVelocityAutocorrelationFunction(traj,atoms,groups,refGroup,normVect=None,
                                           timeInfo=None,diffScheme='fast'):

    file, filename = logfileInit('methyl_axis_avacf')

    tim   = timePrepare(traj,timeInfo)
    avacf = N.zeros(len(tim),N.Float)

    print atoms
       
    for key in groups.keys():
        remove=[]
        for  gr in groups[key]:
             if gr[0] not in atoms: 
               remove.append(gr)
        for r in remove:
            groups[key].remove(r)
      
       
    ngr = 0
    for key in groups.keys(): ngr = ngr + len(groups[key])

    print 'Group names: ',groups.keys()
    print 'Total number of groups is: ',ngr
    
    i = 0
    normFact = ngr
    for key in groups.keys():
        for  gr in groups[key]:
            try:
                if refGroup[key]:
                    group, refs = getReference(gr, refGroup[key])
                else:
                    group = gr
                    refs = None
            except:
                print 'Warning: Exception raised for',gr,i,key
                refs = None
                group = gr
            c_coordinate = traj.readParticleTrajectory(gr[0],
               first=timeInfo[0],last=timeInfo[1],
               skip=timeInfo[2])
            c_parent = gr[0].bondedTo()[0]  
            c_parent_coordinate = traj.readParticleTrajectory(c_parent,
               first=timeInfo[0],last=timeInfo[1],
               skip=timeInfo[2])
            
            mu = N.zeros([len(c_coordinate),3], N.Float)
            for index in range(len(c_coordinate)):
               mu[index] = (c_coordinate[index]-c_parent_coordinate[index]).normal()
            
            res = AutoCorrelationFunction_2(mu)
            avacf = avacf+res/res[0]
            i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    return tim, 1.5*avacf/ngr-0.5

def SpectraOfAngularVACF(traj,groups,refGroup,normVect=None,diffScheme='fast',
                         timeInfo=None,alpha=5.):

    file, filename = logfileInit('savacf')
    ngr = 0
    for key in groups.keys(): ngr = ngr + len(groups[key])

    print 'Group names: ',groups.keys()
    print 'Total number of groups is: ',ngr

    tim   = timePrepare(traj,timeInfo)
    savacf = N.zeros(len(tim),N.Float)
    TimeStep = tim[1] - tim[0]

    total = 0.0
    i = 0
    normFact = ngr
    for key in groups.keys():
        for gr in groups[key]:
            weight = gr.mass()
            try:
                if refGroup[key]:
                    group, refs = getReference(gr,refGroup[key])
                else:
                    group = gr
                    refs = None
            except:
                print 'Warning: Exception raised for',gr,i,key
                refs = None
                group = gr
            angvel = getAngularVelocity(traj,gr,timeInfo,reference=refs,
                                        diffScheme=diffScheme)
            if normVect is None:
                res = N.add.reduce(AutoCorrelationFunction(angvel), -1)/3.
            else: res = AutoCorrelationFunction(N.dot(angvel,normVect))
            savacf = savacf+\
                     weight*fft(GaussianWindow(res,alpha)).real[:len(tim)]
            total = total+weight
            i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    frequencies = N.arange(len(tim))/(2.*len(tim)*TimeStep)
    return frequencies, 0.5*TimeStep*savacf/(ngr*total)

 
def AngularTrajectory(traj,newTrajFN,groups,refGroup,timeInfo=None):
    """ perform fitting and store QUATERNIONS, CMS and FIT. qTrajectory
    object is inherited from Trajectory class with overwritten
    readRigidBodyTrajectory and readParticleTrajectory method. """
 
    file, filename = logfileInit('at')

    traj_new = Trajectory(traj.universe,newTrajFN,"w")
    traj_new.close()

    if timeInfo is None: timeInfo = (0,len(traj)-1,1)
    tim = timePrepare(traj,timeInfo)
    dt  = tim[1]-tim[0]

    ngr = 0
    for key in groups.keys(): ngr = ngr + len(groups[key])

    nc = NetCDFFile(newTrajFN,"a")
    nc.variables['step'][:len(tim)] = N.arange(len(tim))
    btype = None
    try: btype = traj.box_size.typecode()
    except: pass

    if btype:
        box = nc.createVariable('box_size',btype,('step_number',\
                                                  'box_size_length'))
        box[:len(tim)] = getVariables(traj,'box_size')\
                         [timeInfo[0]:timeInfo[1]:timeInfo[2]].astype(btype)

    nc.createDimension('quaternion_length',4)
    nc.createDimension('group_number',ngr)
    nc.createDimension('start',1)
    xyz_atom = nc.createVariable('configuration',N.Float32,\
                                ('start','atom_number','xyz'))
    xyz_cms = nc.createVariable('cms',N.Float32,\
                                ('step_number','group_number','xyz'))
    quatern = nc.createVariable('quaternion',N.Float32,\
                         ('step_number','group_number','quaternion_length'))
    time_new = nc.createVariable('time',N.Float32,('step_number',))
    fit = nc.createVariable('fit',N.Float32,('step_number',
                                                   'group_number'))
    time_new[:len(tim)] = reshape(traj.time[timeInfo[0]:timeInfo[1]:\
                          timeInfo[2]],(len(tim),)).astype(N.Float32)
    xyz_atom[0] = traj.configuration[0].array.astype(N.Float32)

    i = 0
    normFact = ngr
    grdesc = []
    for key in groups.keys():        
        for gr in groups[key]:
            # rbt = traj.readRigidBodyTrajectory(gr,first=timeInfo[0],\
            #           last=timeInfo[1],skip=timeInfo[2])
            # print rbt.cms[:10]
            try:
                if refGroup[key]:
                    group, refs = getReference(gr,refGroup[key])
                else: group, refs = gr, None
            except:
                print 'Warning: Exception raised for',gr,i,key
                refs, group = None, gr
            grdesc.append(map(lambda x: x.index,group.atomList()))
            rbt = traj.readRigidBodyTrajectory(group,first=timeInfo[0],\
                       last=timeInfo[1],skip=timeInfo[2],reference=refs)
            xyz_cms[:,i,:] = rbt.cms.astype(N.Float32)
            quatern[:,i,:] = rbt.quaternions.astype(N.Float32)
            fit[:,i] = rbt.fit.astype(N.Float32)
            i = logfileUpdate(file,i,normFact)

    grdesc = str(grdesc)
    nc.createDimension('group_description_length',len(grdesc))
    group_description = nc.createVariable('group_description','c',
                                          ('group_description_length',))
    group_description[:] = grdesc
    cleanUp(file,filename)
    nc.close()
    return


def RigidBodyTrajectory(traj,newTrajFN,groups,refGroup,timeInfo=None,
                        remove_translation=0):
    """ replace configurations of groups accumulated during MD sim.
    with their rigid reference structure """

    block_size = 100
    file, filename = logfileInit('rbt')
    ngr = 0
    for key in groups.keys(): ngr = ngr + len(groups[key])

    print 'Group names: ',groups.keys()
    print 'Total number of groups is: ',ngr

    if timeInfo is None: timeInfo = (0,len(traj),1)
    if timeInfo[1] is None: timeInfo = (timeInfo[0],len(traj),timeInfo[2])
    tim = timePrepare(traj,timeInfo)
    TimeStep  = tim[1]-tim[0]
    traj_new = Trajectory(traj.universe,newTrajFN,'w',block_size=block_size)

    snapshot = SnapshotGenerator(traj.universe,
                                 actions = [TrajectoryOutput(traj_new, None,
                                 timeInfo[0],timeInfo[1],timeInfo[2])])
    snapshot()
    traj_new.close()

    nc_new = NetCDFFile(newTrajFN, "a")
    nsteps = len(tim)
    step_info = traj.step[timeInfo[0]:timeInfo[1]:timeInfo[2]].astype('l')
    for i in range(nsteps):
        nc_new.variables['step'][i / block_size, i % block_size] = step_info[i]

    nsteps_major = (nsteps+block_size-1) / block_size
    nsteps_minor = nsteps % block_size
    totsteps = block_size*nsteps_major
    if not nc_new.variables.has_key('time'):
	nc_new.createVariable('time', N.Float32,
                              ('step_number', 'minor_step_number'))
    time_data = N.zeros((totsteps,), N.Float32)
    time_data[:nsteps] = traj.time[timeInfo[0]:timeInfo[1]:
                                   timeInfo[2]].astype(N.Float32)
    time_data.shape = nsteps_major, block_size
    nc_new.variables['time'][:, :] = time_data

    try:
        btype = traj.box_size.typecode()
    except AttributeError:
        btype = None
    if btype:
        box_size = N.zeros((totsteps, 3), btype)
        box_size[:nsteps] = getVariables(traj,'box_size') \
                            [timeInfo[0]:timeInfo[1]:timeInfo[2]].astype(btype)
        nc_new.variables['box_size'].shape
        for j in range(nsteps_major):
            nc_new.variables['box_size'][j, :, :] = \
                     N.transpose(box_size[j*block_size:(j+1)*block_size, :])

    i = 0 
    normFact = ngr + len(tim)
    nc_new.variables['configuration'].shape
    for key in groups.keys():
        for gr in groups[key]:
            try:
                if refGroup[key]:
                    group, refs = getReference(gr,refGroup[key])
                else:
                    group = gr
                    refs = None
            except:
                refs = None
                group = gr
            rbt = traj.readRigidBodyTrajectory(object=group,first=timeInfo[0],
                       last=timeInfo[1],skip=timeInfo[2],reference=refs)
            D = rotMatrix(rbt.quaternions)
            offset = traj.universe.contiguousObjectOffset([group],
                                             traj.configuration[timeInfo[0]])
            #Dt = transpose(D,(0,2,1))
            cms = group.centerOfMass(refs).array
            for atom in group.atomList():
                if refs is None:
                    xyz = atom.position().array
                    if offset is not None:
                        xyz = xyz + offset[atom].array
                else:
                    xyz = refs[atom].array
                xyz = xyz - cms
                xyz_r = N.zeros((totsteps, 3), N.Float32)
                if remove_translation:
                    xyz_r[:nsteps] = (N.dot(D, xyz)+rbt.cms[0,:]).astype(N.Float32)
                else:
                    xyz_r[:nsteps] = (N.dot(D, xyz)+rbt.cms).astype(N.Float32)
                for j in range(nsteps_major):
                    nc_new.variables['configuration'][j, atom.index, :, :] = \
                         N.transpose(xyz_r[j*block_size:(j+1)*block_size, :])
            i = logfileUpdate(file,i,normFact)

    nc_new.close()
    traj_new = Trajectory(None,newTrajFN,'r')
    nc_new = NetCDFFile(newTrajFN,'a')
    for i in range(len(traj_new)):
        conf = traj_new.configuration[i]
        traj_new.universe.setConfiguration(conf)
        traj_new.universe.foldCoordinatesIntoBox()
        conf_new = traj_new.universe.configuration()
        nc_new.variables['configuration'][i/block_size, :, :, i%block_size] \
                                       = conf_new.array.astype('f')
        i = logfileUpdate(file,i,normFact)
    nc_new.close()
    traj_new.close()
    cleanUp(file,filename)
    
    return


def ReorientationalCorrelationFunction(traj,jmn,groups,refGroup,timeInfo=None):

    file, filename = logfileInit('rcf')

    ngr = 0
    for key in groups.keys(): ngr = ngr + len(groups[key])

    print 'Group names: ',groups.keys()
    print 'Total number of groups is: ',ngr

    tim = timePrepare(traj,timeInfo)
    res = N.zeros((len(tim)),N.Float)

    i = 0
    normFact = ngr
    for key in groups.keys():
        for  gr in groups[key]:
            try:
                if refGroup[key]:
                    group, refs = getReference(gr,refGroup[key])
                else:
                    group = gr
                    refs = None
            except:
                refs = None
                group = gr
            sh  = sphericalHarmonics(traj,jmn,group,refs,timeInfo)
            res = res + CorrelationFunction(sh[0],sh[1])
            i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)
    return tim, res*4*N.pi/ngr


def DigitalFilter(traj,newname,groups=None,timeInfo=None,filterSet=None):
 
    file, filename = logfileInit('df')
    if  newname == traj.filename: newname = traj.filename + ".df"

    gj = []
    if groups is None: groups = traj.universe._objects.objects
    else:
        for key in groups.keys(): gj = gj + groups[key]
        groups = gj

    keep = Collection()
    for object in groups: keep.addObject(object)
    traj_new = Trajectory(keep,newname,"w") 
    alist = tuple(keep.atomList())
    traj_new.close()

    if timeInfo is None: timeInfo = (0,len(traj)-1,1)
    tim      = timePrepare(traj,timeInfo)
    dos      = N.zeros((len(tim)-2),N.Float)
    TimeStep = tim[1]-tim[0]
    frequencies = N.arange(len(tim))/(2*len(tim)*TimeStep)
    F = N.zeros(len(frequencies),'d')
    if filterSet is None: # by default we take only 5%
        filterSet = (0,frequencies[len(frequencies)/20])
    fmin = len(frequencies)-len(compress(greater_equal(frequencies,\
                                filterSet[0]),frequencies))
    if filterSet[1] is None:
        fmax = len(frequencies)
    else:
        fmax = len(frequencies)-len(compress(greater_equal(frequencies,\
                                    filterSet[1]),frequencies))
    F[fmin:fmax] = 1.0
    filter = N.zeros(2*len(F)-2,'d')
    filter[:len(F)] = F
    filter[len(F):] = F[-2:0:-1]

    nc = NetCDFFile(traj.filename,"r")
    ctype = nc.variables['configuration'].typecode()
    ttype = nc.variables['time'].typecode()
    btype = None
    try: 
        btype = nc.variables['box_size'].typecode()
        bsize = nc.dimensions['box_size_length']
    except: pass
    nc.close()

    nc = NetCDFFile(newname,"a")
    nc.variables['step'][:len(tim)] = N.arange(len(tim))
    coord = nc.createVariable('configuration',ctype,\
                              ('step_number','atom_number','xyz'))
    time_new = nc.createVariable('time',ttype,('step_number',))
    if btype:
        box = nc.createVariable('box_size',btype,('step_number',\
                                                  'box_size_length'))
        box[:len(tim)] = getVariables(traj,'box_size')\
               [timeInfo[0]:timeInfo[1]:timeInfo[2]].astype(btype)

    atom_index = nc.createVariable('atom_index','i',('atom_number',))
    time_new[:len(tim)] = traj.time[timeInfo[0]:timeInfo[1]:timeInfo[2]]\
                              .astype(ttype)
    
    xold = N.zeros(len(frequencies)*2-2,N.Float)
    xnew = N.zeros((len(frequencies),3),N.Float)
  #  before = N.zeros(len(frequencies),N.Float)
  #  after = N.zeros(len(frequencies),N.Float)

    i = 0
    normFact = len(alist) + len(tim)
    for at in range(len(alist)): 
        xyz = traj.readParticleTrajectory(alist[at],first=timeInfo[0],\
                   last=timeInfo[1],skip=timeInfo[2]).array
        for ia in range(3):
            x = xyz[:,ia]
            xold[:len(x)] = x
            xold[len(x):] = x[-2:0:-1]
            xfft = fft(xold*N.sqrt(alist[at]._mass))
          #  before = before + abs(xfft)[:len(x)]**2
          #  after = after + abs(filter*xfft)[:len(x)]**2
            xnew[:,ia] = inverse_fft(filter*xfft)[:len(x)].real/\
                         N.sqrt(alist[at]._mass)
        coord[:,at,:] = xnew.astype(ctype)
        atom_index[at] = alist[at].index
        i = logfileUpdate(file,i,normFact)

 #   dos_before = TimeStep*TimeStep*frequencies*frequencies*before
 #   dos_after  = TimeStep*TimeStep*frequencies*frequencies*after
    #
    # todo: (1) filter function smoothing, (2) dos normalization,
    # (3) dos weightening
    #
    nc.close()
    #
    # wrapping atoms, whose trajectory was made continuoues,
    # back into box
    # MINOR_STEP: use Trajectory.SnaphotGenerator...
    #
    traj_new = Trajectory(None,newname,'r')
    nc = NetCDFFile(newname,'a')
    for i in range(len(traj_new)):
        conf = traj_new.configuration[i]
        traj_new.universe.setConfiguration(conf)
        traj_new.universe.foldCoordinatesIntoBox()
        conf_new = traj_new.universe.configuration()
        nc.variables['configuration'][i] = conf_new.array.astype('f')
        i = logfileUpdate(file,i,normFact)
    nc.close()
    traj_new.close()
    cleanUp(file,filename)
    
    return


class AngularTrajectoryCache:

    def __init__(self, trajectory, size, first, last, skip):
        self.trajectory = trajectory
        self.size = size
        self.first = first
        self.last = last
        self.skip = skip
        self.data = {}
        self.keys = []
        self.count = 0

    def __call__(self, group, reference):
        key = (group, reference)
        data = self.data.get(key, None)
        if data is None:
            if len(self.keys) == self.size:
                del self.data[self.keys[0]]
                del self.keys[0]
            data = self.trajectory.readRigidBodyTrajectory(group,
                                                           self.first,
                                                           self.last,
                                                           self.skip,
                                                           reference)
            self.count = self.count + 1
            self.keys.append(key)
            self.data[key] = data
        return data


def MolecularPairCorrelationFunction(traj,groups_A,refGroup_A,groups_B=None,
                                     refGroup_B=None,symbols=(0,0,0,0,0),
                                     timeInfo=None,rLimit=None,width=.02,
                                     todo=None):

    l1, l2, m1, m2, n1 = symbols
    file, filename = logfileInit('mpcf')

    if timeInfo is None: timeInfo = (0,len(traj)-1,1)
    tim   = timePrepare(traj,timeInfo)
    TimeStep = tim[1] - tim[0]

    if refGroup_A is None:
        typs = getTypes(traj.universe)
        refGroup_A = {}
        for ityp in typs.keys(): refGroup_A[ityp] = None
        
    if not refGroup_B: refGroup_B = refGroup_A
    if not groups_B: groups_B = groups_A

    def validateRefs(reference,group,traj,timeInfo):

        for ref in group.keys():
            if reference[ref] is None: # the first configuration as a reference
                if traj.variables().count('quaternion') == 1: iref = 0
                else: iref = timeInfo[0]
                traj.universe.setConfiguration(\
                    traj.universe.contiguousObjectConfiguration(\
                    group[ref],traj.configuration[iref]))
                reference[ref] = Collection(group[ref]).atomList()
        return reference
    
    refGroup_A = validateRefs(refGroup_A,groups_A,traj,timeInfo)
    refGroup_B = validateRefs(refGroup_B,groups_B,traj,timeInfo)

    if not groups_B: groups_B = groups_A

    gA =  gB = []
    for objects in groups_A.values(): gA = gA + objects
    for objects in groups_B.values(): gB = gB + objects
    nA = len(gA)
    nB = len(gB)
    #ngr = len(gA) + len(gB)

    print 'Group names: %s (a), %s (b)' % (str(groups_A.keys()),
                                           str(groups_B.keys()))
    print '# groups: %d (a), %d (b)' % (len(gA),len(gB))

    ez = N.array([0.,0.,1.])

    i = 0
    gj = 1.8*traj.universe.largestDistance()
    if rLimit is None:
        if gj: rLimit = gj
        else:  rLimit = 1.5
    bins = N.arange(0.,rLimit,width)
    done = []
    mpcf = N.zeros(len(bins)-1,N.Float)
    box = getVariables(traj,'box_size')\
          [timeInfo[0]:timeInfo[1]:timeInfo[2]].astype(N.Float)
    hbox = 0.5*box
    done = []
    if todo is None:
        todo = uniqueGroupPairs(groups_A,groups_B,refGroup_A,refGroup_B)
    normFact = N.add.reduce(map(len, todo.values()))

    cache_size = 11
    cache = AngularTrajectoryCache(traj, cache_size, timeInfo[0], timeInfo[1],
                                   timeInfo[2])

    group_list = todo.keys()
    group_list.sort()
    diagonal = 0
    while 1:
        done = 1
        for row in range(0, len(group_list)):
            ga, ra = group_list[row]
            partner_list = todo[(ga,ra)]
            if len(partner_list) < diagonal:
                continue
            done = 0
            rbt_A = cache(ga, ra)
            for column in range(diagonal, min(len(partner_list),
                                              diagonal+cache_size-1)):
                gb, rb = partner_list[column]
                rbt_B = cache(gb, rb)
                Rab = rbt_B.cms - rbt_A.cms
                Rab = fmod(Rab+hbox, box)-hbox
                Rab = fmod(Rab-hbox, box)+hbox
                R = N.sqrt(N.add.reduce(Rab*Rab,-1))
                nab = Rab/R[:,N.NewAxis]
                nq = nab + ez
                norm = N.sqrt(N.add.reduce(nq*nq, -1))
                mask = greater(norm, 1.e-10)
                imask = 1-mask
                nq = mask[:,N.NewAxis]*(nq/(norm+imask)[:,N.NewAxis])
                nq = nq + imask[:,N.NewAxis]*ez[N.NewAxis,:]
                #
                Qa = N.zeros((len(tim), 4), N.Float)
                qa = rbt_A.quaternions
                Qa[:, 0] = N.add.reduce(qa[:, 1:]*nq, -1)
                Qa[:, 1] = -qa[:, 0]*nq[:, 0] \
                           -qa[:, 2]*nq[:, 2] \
                           +qa[:, 3]*nq[:, 1]
                Qa[:, 2] = - qa[:, 0]*nq[:, 1] \
                           + qa[:, 1]*nq[:, 2] \
                           - qa[:, 3]*nq[:, 0]
                Qa[:, 3] = - qa[:, 0]*nq[:, 2] \
                           - qa[:, 1]*nq[:, 1] \
                           + qa[:, 2]*nq[:, 0]
                Qb = N.zeros((len(tim), 4), N.Float)
                qb = rbt_B.quaternions
                Qa[:, 0] = N.add.reduce(qa[:, 1:]*nq, -1)
                Qb[:, 1] = -qb[:, 0]*nq[:, 0] \
                           -qb[:, 2]*nq[:, 2] \
                           +qb[:, 3]*nq[:, 1]
                Qb[:, 2] = - qb[:, 0]*nq[:, 1] \
                           + qb[:, 1]*nq[:, 2] \
                           - qb[:, 3]*nq[:, 0]
                Qb[:, 3] = - qb[:, 0]*nq[:, 2] \
                           - qb[:, 1]*nq[:, 1] \
                           + qb[:, 2]*nq[:, 0]
                #
                A1 = WignerFunctions(l1,m1,n1,Qa)
                B1 = WignerFunctions(l2,m2,-n1,Qb)
                A2 = WignerFunctions(l1,m1,n1,Qb)
                B2 = WignerFunctions(l2,m2,-n1,Qa)
                G = (A1*B1+A2*B2).real

                indices = argsort(R)
                R = take(R, indices)
                G = take(G, indices)
                indices = N.concatenate((searchsorted(R, bins), [len(R)]))
                h = map(lambda i1, i2, w=G: N.add.reduce(operator.getslice(w,
                                                                         i1,
                                                                         i2)),
                        indices[:-1], indices[1:])
                mpcf = mpcf + h[:-1]
                i = logfileUpdate(file,i,normFact)
        if done: break
        diagonal = diagonal + cache_size-1

    #print cache.count, "trajectory reads"
    volume = N.add.reduce(box[:,0]*box[:,1]*box[:,2])/len(tim)
    pairs  = normFact*len(tim)
    r = 0.5*(bins[1:]+bins[:-1])
    shell_volume = (4.*N.pi/3.)*(bins[1:]**3-bins[:-1]**3)
    factor = (volume * (2*l1+1) * (2*l2+1)) / \
             (shell_volume * nA*nB * len(tim))
    mpcf = mpcf*factor

    cleanUp(file,filename)
    return r, mpcf


def ScatteringFunctionFFT(ncResultsFN,ncDataFN,nfreq,alpha=5.):
    """ perform fft on previously calculated scattering functions """

    file, filename = logfileInit('sf_fft')

    ncin = NetCDFFile(ncResultsFN,'r')
    tim = ncin.variables['time']
    TimeStep = tim[1] - tim[0]
    frequencies = N.arange(len(tim))/(2.*len(tim)*TimeStep)
    sample = max(1, len(frequencies)/nfreq)
    nfreq = (len(frequencies)+sample-1)/sample
    q = ncin.variables['q'][:]
    fqt = ncin.variables['Fqt-total'][:,:]
    ncin.close()
    
    ncout = NetCDFFile(ncDataFN,'a')
    ncout.title = 'Dynamic Structure Factor'
    ncout.createDimension('FREQUENCY',nfreq)
    SF    = ncout.createVariable('Sqw-total',N.Float, ('Q','FREQUENCY'))
    FREQ  = ncout.createVariable('frequency',N.Float, ('FREQUENCY',))

    FREQ[:] = frequencies[::sample]
    
    i = 0
    normFact = len(fqt)
    for qlen in range(normFact):
        SF[qlen] = 0.5*TimeStep * \
                   fft(GaussianWindow(fqt[qlen],alpha)).real[:nfreq*sample:sample]
        i = logfileUpdate(file,i,normFact)

    cleanUp(file,filename)

    ncout.close()
    return


## def RadialDistributionFunction(traj,groups_A,groups_B=None,timeInfo=None,
##                                 rLimit=None,width=.02):

##     file, filename = logfileInit('rdf')
##     if not groups_B: groups_B = groups_A
##     gA =  gB = []
##     for key in groups_A.keys(): gA = gA + groups_A[key]
##     for key in groups_B.keys(): gB = gB + groups_B[key]
##     ngr = len(gA) + len(gB)

##     print 'Group names: ',groups_A.keys(),groups_B.keys()
##     print 'Total number of groups ',ngr

##     if timeInfo is None: timeInfo = (0,len(traj)-1,1)
##     tim   = timePrepare(traj,timeInfo)
##     TimeStep = tim[1] - tim[0]

##     i = 0
##     normFact = len(gA)*(len(gB)-1)/2
##     gj = traj.universe.largestDistance()
##     if rLimit is None:
##         if gj: rLimit = gj
##         else:  rLimit = 1.5
##     bins = arange(0.,rLimit,width)
##     done = []
##     rdf = N.zeros(len(bins)-1,N.Float)
##     box = getVariables(traj,'box_size')\
##           [timeInfo[0]:timeInfo[1]:timeInfo[2]].astype(N.Float)
##     done = []
##     todo = uniqueGroupPairs(gA,gB)

##     for igA in todo.keys():
##         rbt_A = traj.readRigidBodyTrajectory(object=igA,first=timeInfo[0],
##                           last=timeInfo[1],skip=timeInfo[2],reference=None)
##         for igB in todo[igA]:
##             rbt_B = traj.readRigidBodyTrajectory(object=igB,
##                               first=timeInfo[0],last=timeInfo[1],
##                               skip=timeInfo[2],reference=None)
##             Rab = (rbt_B.cms - rbt_A.cms).astype(N.Float)
##             # PBC for Rab
##             for ia in range(3):
##                 iR = Rab[:,ia]
##                 iB = box[:,ia]/2.
##                 Rab[:,ia] = where(greater_equal(iR,iB),iR-2.*iB,iR)
##                 Rab[:,ia] = where(less_equal(iR,-iB),iR+2.*iB,iR)
##             R = sqrt(add.reduce(Rab*Rab,-1))
##             dax = histogram(R,bins)[:-1]
##             rdf = rdf + dax
##             i = logfileUpdate(file,i,normFact)

##     pairs = normFact*len(tim)
##     volume = add.reduce(box[:,0]*box[:,1]*box[:,2])/len(tim)
##     factor = 4./3. * pi * pairs * (bins[1:]**3 - bins[:-1]**3)/volume
##     rdf = rdf/factor
##     return bins[:-1],rdf


def AutoRegressiveAnalysis(memory, correlation, spectrum, msd,
                           traj,atoms,timeInfo,order,nsteps,nfreq,
                           weightList=None, precision=None):

    from Scientific.Signals.Models import AutoRegressiveModel, \
                                          AveragedAutoRegressiveModel
    if weightList is None:
        weightList = OnesList(traj.universe,atoms.atomList())
    t = traj.time[timeInfo[0]:2*timeInfo[2]:timeInfo[2]]
    dt = t[1]-t[0]
    file, filename = logfileInit('ar-vel')

    model = AveragedAutoRegressiveModel(order, dt)
    i = 0
    normFact = len(atoms.atomList())
    for at in atoms.atomList():
        series = traj.readParticleTrajectory(at,first=timeInfo[0],
                        last=timeInfo[1],skip=timeInfo[2],
                        variable='velocities').array
        for j in range(3):
            model.add(AutoRegressiveModel(order, series[:, j], dt),
                      weightList[at])
        i = logfileUpdate(file, i, normFact)

    cleanUp(file,filename)
    return evaluate_model(model, nsteps, nfreq, memory, correlation,
                          spectrum, msd, precision)


def AutoRegressiveAnalysisXYZ(memory, correlation, spectrum, msd,
                              traj,atoms,timeInfo,order,nsteps,nfreq,
                              weightList=None, precision=None,
                              diffScheme='fast'):

    from Scientific.Signals.Models import AutoRegressiveModel, \
                                          AveragedAutoRegressiveModel
    if weightList is None:
        weightList = OnesList(traj.universe,atoms.atomList())
    t = traj.time[timeInfo[0]:timeInfo[1]:timeInfo[2]]
    
    nsteps = len(t)
    nfreq = len(t)
    
    dt = t[1]-t[0]
    file, filename = logfileInit('ar-xyz')

    model = AveragedAutoRegressiveModel(order, dt)
    i = 0
    normFact = len(atoms.atomList())
    for at in atoms.atomList():
        series = traj.readParticleTrajectory(at,first=timeInfo[0],
                        last=timeInfo[1],skip=timeInfo[2]).array
        for j in range(3):
            v = diff(series[:,j],diffScheme,dt)
            model.add(AutoRegressiveModel(order, v, dt), weightList[at])
        i = logfileUpdate(file, i, normFact)

    cleanUp(file,filename)
    return evaluate_model(model, nsteps, nfreq, memory, correlation,
                          spectrum, msd, precision)

def evaluate_model(model, nsteps, nfreq, memory, correlation, spectrum, msd,
                   precision):
        
    if correlation:
        correlation = model.correlation(nsteps)
        c = correlation.values
        c_init = N.absolute(c[0])
        for i in range(len(c)):
            if N.absolute(c[i])/c_init < 1.e-10:
                break
        if i < len(c):
            correlation = InterpolatingFunction((correlation.axes[0][:i],),
                                                c[:i])
    if spectrum:
        omega_max = 1.1*N.pi/model.delta_t
        omega = omega_max*N.arange(nfreq)/float(nfreq)
        spectrum = model.spectrum(omega)
        spectrum = InterpolatingFunction((omega/(2.*N.pi),), spectrum.values)
    if msd:
        msd = msd_eval(model, nsteps)
    if memory:
        friction = model.frictionConstant()
        if precision is not None:
            from ARMP import MPAutoRegressiveModel
            model = MPAutoRegressiveModel(model, precision)
        try:
            memory = model.memoryFunction(model.order+model.order/2)
        except OverflowError:
            null = N.zeros((0,), N.Float)
            memory = InterpolatingFunction((null,), null)
        memory.friction = friction
    return model, memory, correlation, spectrum, msd

def msd_eval(model, nsteps):
    poles = model.poles()
    cpoles = N.conjugate(poles)
    coeff0 = N.conjugate(model.coeff[0])
    beta = N.zeros((model.order,), N.Complex)
    for i in range(model.order):
        pole = poles[i]
        beta[i] = -(model.sigsq*pole**(model.order-1)/coeff0) / \
                  (N.multiply.reduce((pole-poles)[:i]) *
                   N.multiply.reduce((pole-poles)[i+1:]) *
                   N.multiply.reduce(pole-1./cpoles) *
                   model.variance)
    beta = beta/N.sum(beta)
    msd = N.zeros((nsteps,), N.Float)
    n = N.arange(nsteps)
    for i in range(model.order):
        pole = poles[i]
        msd = msd + (beta[i]*((pole**n-1.)*pole/(1-pole)**2
                              + n/(1-pole))).real
    msd = 6. * model.delta_t**2 * model.variance * msd
    return InterpolatingFunction((model.delta_t*n,), msd)
    
def writeARParameters(model, filename):
    file = TextFile(filename, 'w')
    file.write('Variance: %20.15f\n' % model.variance)
    file.write('Sigma: %20.15f\n' % model.sigma)
    for i in range(model.order):
        file.write('A%d: %20.15f\n' % (i+1, model.coeff[-i-1]))
    file.close()
