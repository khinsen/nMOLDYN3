import sys
sys.path.insert(0,"C:\\Datas\\Eclipse_Projects\\Python\\nMOLDYN\\trunk")

from getopt import *
import os
from os import _exit
from string import strip
from sys import argv, stderr
import tempfile
from time import sleep,asctime,localtime,time

from Scientific import N

from MMTK import Units
from MMTK.Proteins import defineAminoAcidResidue

from nMOLDYN.Core.IO import ASCII_To_NetCDF

from nMOLDYN.Tests.nMOLDYN_Reference.core import AngularVelocityAutocorrelationFunction, \
                                                 AutoRegressiveAnalysis, \
                                                 AutoRegressiveAnalysisXYZ, \
                                                 CoherentScatteringAR, \
                                                 CoherentScatteringFunction, \
                                                 DensityOfStates_xyz, \
                                                 ElasticIncoherentStructureFactor, \
                                                 IncoherentScatteringAR, \
                                                 IncoherentScatteringFunction, \
                                                 IncoherentScatteringFunctionGaussian, \
                                                 MeanSquareDisplacement, \
                                                 ReorientationalCorrelationFunction, \
                                                 ScatteringFunctionFFT, \
                                                 SpectraOfAngularVACF, \
                                                 VelocityAutocorrelationFunction_xyz
                                                 
from nMOLDYN.Tests.nMOLDYN_Reference.misc import inputFileRead,saveText,saveNetCDF,ghostBusters, tup2dict,getTypes,getChemicalObjects,getProgress

control = ['ar-vel','ar-xyz','dos-vel','dos-xyz','msd','vacf-vel','vacf-xyz',
           'csf','arcsf','eisf','isf','arisf','cdaf-vel','cdaf-xyz','df',
           'isfg','at','avacf','rcf','input=','progress','savacf','fft','mpcf',
           'rbt','rbrt','methyl_avacf', 'methyl_axis_avacf',
           'methyl_aligned_avacf']
optlist = tuple(getopt(argv[1:],'',control)[0])
optlist = tup2dict(optlist)

start_time = time()
print """
 --------------------------
  pMoldyn Python Version
 --------------------------
"""
print " Started at: ",asctime(localtime(start_time))

if optlist:
    data = None
    if optlist.has_key('--input'):
        input = optlist['--input']
        file = open(input,'r')
        print ' Processing input file: \n'
        print ' ---'
        for line in file.readlines(): print line[:-1]
        print '\n ---'
        file.close()
        data = inputFileRead(input)
        print 'Done ...\n ---\n'
    elif optlist.has_key('--progress'): pass
    else:
        print ' Input file is needed '
        _exit(0)

    if data is not None:
        uL       = data.units_length
        uF       = data.units_frequency
        output   = data.output_files
        title    = data.title
        if not data.trajectory: scatFunc = data.results_file
        else:
            traj     = data.trajectory
            qVect    = data.q_vector_set
            types    = None
            if data.groups:    groups = data.groups
            else: groups = [{'Universe': traj.universe._objects.objects}]
            if data.reference: refGroup = data.reference
            else: refGroup = None
            if data.atoms: atoms = data.atoms
            else: atoms = traj.universe
            weights  = data.weights
            timeInfo = data.time_info
            if data.projection_vector: normVect = N.array(data.projection_vector)
            else: normVect = None
            if data.rotation_coefficients: rotCoef = data.rotation_coefficients
            else: rotCoef = (0.0,0.0,0.0)
            if data.symbols: symbols = data.symbols
            else: symbols = (0,0,0,0,0)
            if data.filter_window: setFilter = data.filter_window
            if data.differentiation: diffScheme = data.differentiation
            else: diffScheme = 'fast'
            if data.ft_window is None: alpha = 10.
            else: alpha = 100./data.ft_window
            print " Input file processed in : %.1f sec." % (time()-start_time)
    elif optlist.has_key('--progress'): pass
    else:
        sys.stderr.write('Problem with input file\n')
        raise SystemExit
    
    if optlist.has_key('--dos-vel'):
        try:
            output = output['dos']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        DOS = DensityOfStates_vel(traj,alpha=alpha,atoms=atoms,
                                  normVect=normVect,
                                  weightList=weights,timeInfo=timeInfo)
        sample = max(1, len(DOS[0])/data.frequency_points)
        dos_scaled = []
        dos_scaled.append(N.array(DOS[0])[::sample]/uF)
        dos_scaled.append(N.array(DOS[1])[::sample]/(uL*uL))
        saveText(data=N.transpose(N.array(dos_scaled)),filename=output,title=title)
    elif optlist.has_key('--dos-xyz'):
        try:
            output = output['dos']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        DOS = DensityOfStates_xyz(traj,alpha=alpha,atoms=atoms,
                                  normVect=normVect,
                                  weightList=weights,timeInfo=timeInfo)
        sample = max(1, len(DOS[0])/data.frequency_points)
        dos_scaled = []
        dos_scaled.append(N.array(DOS[0])[::sample]/uF)
        dos_scaled.append(N.array(DOS[1])[::sample]/(uL*uL))
        saveText(data=N.transpose(N.array(dos_scaled)),filename=output,title=title)
        
        ASCII_To_NetCDF(output, \
                        os.path.join(tempfile.gettempdir(),'CDOS_Reference.nc'), \
                        varname = ["frequency", "dos-total"], \
                        skiprows = 3, bycolumn = True)
        
    elif optlist.has_key('--vacf-vel'):
        try:
            output = output['vacf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        VACF = VelocityAutocorrelationFunction_vel(traj,
                       atoms=atoms,normVect=normVect,
                       weightList=weights,timeInfo=timeInfo)
        saveText(data=N.transpose(N.array(VACF)[:data.time_steps]),
                 filename=output,title=title)
    elif optlist.has_key('--vacf-xyz'):
        try:
            output = output['vacf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        VACF = VelocityAutocorrelationFunction_xyz(traj,
                       atoms=atoms,normVect=normVect,
                       weightList=weights,timeInfo=timeInfo)
        saveText(data=N.transpose(N.array(VACF))[:data.time_steps],
                 filename=output,title=title)

        ASCII_To_NetCDF(output, \
                        os.path.join(tempfile.gettempdir(),'CVACF_Reference.nc'), \
                        varname = ["time", "vacf-total"], \
                        skiprows = 3, bycolumn = True)
        
    elif optlist.has_key('--msd'):
        try:
            output = output['msd']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        MSD = MeanSquareDisplacement(traj,atoms,timeInfo,
                                     weightList=weights,normVect=normVect)
        msd_scaled = []
        msd_scaled.append(N.array(MSD[0])[:data.time_steps])
        msd_scaled.append(N.array(MSD[1])[:data.time_steps]/(uL*uL))
        saveText(data=N.transpose(N.array(msd_scaled)),filename=output,title=title)

        ASCII_To_NetCDF(output, \
                        os.path.join(tempfile.gettempdir(),'MSD_Reference.nc'), \
                        varname = ["time", "msd-total"], \
                        skiprows = 3, bycolumn = True)
        
    elif optlist.has_key('--ar-xyz') or optlist.has_key('--ar-vel'):
        order = data.ar_order
        precision = data.ar_precision
        ara = output.get('ara', None)

        if optlist.has_key('--ar-xyz'):
            model, mem, corr, spect, msdd = \
                 AutoRegressiveAnalysisXYZ(True, True, True, True,
                                           traj,atoms,timeInfo,order,
                                           data.time_steps,
                                           data.frequency_points,
                                           weightList=weights,
                                           precision=precision,
                                           diffScheme=diffScheme)
        else:
            model, mem, corr, spect, msdd = \
                 AutoRegressiveAnalysis(True, True, True, True,
                                        traj,atoms,timeInfo,order,
                                        data.time_steps,
                                        data.frequency_points,
                                        weightList=weights,
                                        precision=precision)
                                                                                  
        corr.values = corr.values.real/(uL*uL)
        saveText(data = corr, filename = ara, title = title)
        ASCII_To_NetCDF(ara, \
                        os.path.join(tempfile.gettempdir(),'ARA_Reference.nc'), \
                        varname = ["time_vacf", "vacf"], dimname = "vacf_dim", \
                        skiprows = 3, bycolumn = True)
                
        spect.values = spect.values.real/(uL*uL)
        spect.axes = (spect.axes[0]/uF,)
        saveText(data = spect, filename = ara, title = title)
        ASCII_To_NetCDF(ara, \
                        os.path.join(tempfile.gettempdir(),'ARA_Reference.nc'), \
                        varname = ["frequency", "dos"], dimname = "dos_dim", \
                        skiprows = 3, bycolumn = True, append = True)
        
        msdd.values = msdd.values.real/(uL*uL)
        saveText(data = msdd, filename = ara, title = title)
        ASCII_To_NetCDF(ara, \
                        os.path.join(tempfile.gettempdir(),'ARA_Reference.nc'), \
                        varname = ["time_msd", "msd"], dimname = "msd_dim",\
                        skiprows = 3, bycolumn = True, append = True)
        
        mem.values = mem.values.real/(uL*uL)
        saveText(data = mem, filename = ara, title = title + ('\n# Friction constant: %20.15f'
                                      % mem.friction))
        ASCII_To_NetCDF(ara, \
                        os.path.join(tempfile.gettempdir(),'ARA_Reference.nc'), \
                        varname = ["time_memory", "memory_function"], dimname = "mem_dim",\
                        skiprows = 4, bycolumn = True, append = True)
        
    elif optlist.has_key('--csf'):
        try:
            file = output['csf']
        except KeyError:
            sys.stderr.write('No csf output file specified.\n')
            raise SystemExit
        CSF = CoherentScatteringFunction(traj,qVect=qVect,ncDataFN=file,
                                         atoms=atoms,bcoh=weights,
                                         timeInfo=timeInfo,
                                         nsteps=data.time_steps)
        
        ScatteringFunctionFFT(file,output['csf'],data.frequency_points,alpha=alpha)
        
    elif optlist.has_key('--arcsf'):
        CoherentScatteringAR(traj, output.get('csf', None),
                             output.get('memory', None),
                             output.get('fft', None), data.ar_order,
                             qVect, atoms, weights, timeInfo,
                             data.time_steps, data.frequency_points,
                             data.ar_precision)
        
    elif optlist.has_key('--arisf'):
        IncoherentScatteringAR(traj, output.get('isf', None),
                               output.get('memory', None),
                               output.get('fft', None), data.ar_order,
                               qVect, atoms, weights, timeInfo,
                               data.time_steps, data.frequency_points,
                               data.ar_precision)        
    elif optlist.has_key('--eisf'):
        try:
            output = output['eisf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        EISF = ElasticIncoherentStructureFactor(traj,qVect=qVect,
                                                atoms=atoms,bincoh=weights,
                                                timeInfo=timeInfo)
        saveText(data=N.transpose(N.array(EISF)),filename=output,title=title)
        
        ASCII_To_NetCDF(output, \
                        os.path.join(tempfile.gettempdir(),'EISF_Reference.nc'), \
                        varname = ["q", "eisf-total"], \
                        skiprows = 3, bycolumn = True)
        
    elif optlist.has_key('--isf'):
        try:
            file = output['isf']
        except KeyError:
            sys.stderr.write('No isf output file specified.\n')
            raise SystemExit
        ISF = IncoherentScatteringFunction(traj,qVect=qVect,
                                           ncDataFN=file,atoms=atoms,
                                           bincoh=weights,timeInfo=timeInfo,
                                           nsteps=data.time_steps)

        ScatteringFunctionFFT(file,output['isf'],data.frequency_points,alpha=alpha)
        
    elif optlist.has_key('--isfg'):
        try:
            file = output['isf']
        except KeyError:
            sys.stderr.write('No isf output file specified.\n')
            raise SystemExit
        ISFG = IncoherentScatteringFunctionGaussian(traj,qVect=qVect,
                                                    ncDataFN=file,
                                                    atoms=atoms,bincoh=weights,
                                                    timeInfo=timeInfo,
                                                    nsteps=data.time_steps)

        ScatteringFunctionFFT(file,output['isf'],data.frequency_points, alpha=alpha)
        
    elif optlist.has_key('--cdaf-vel'):
        try:
            file = output['cdaf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        CDAF = CurrentDensityAutocorrelationFunction(traj,qVect=qVect,
                                                     ncDataFN=file,
                                                     atoms=atoms,bcoh=weights,
                                                     timeInfo=timeInfo,
                                                     nsteps=data.time_steps)
    elif optlist.has_key('--avacf'):
        try:
            output = output['avacf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        AVACF = AngularVelocityAutocorrelationFunction(traj,groups[0],
                               refGroup[0],normVect=normVect,
                               timeInfo=timeInfo,diffScheme=diffScheme)
        saveText(data=N.transpose(N.array(AVACF)[:data.time_steps]),
                 filename=output,title=title)

        ASCII_To_NetCDF(output, \
                        os.path.join(tempfile.gettempdir(),'AVACF_Reference.nc'), \
                        varname = ["time", "avacf-total"], \
                        skiprows = 3, bycolumn = True)
        
    elif optlist.has_key('--methyl_avacf'):
        try:
            output = output['methyl_avacf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        AVACF = MethylAngularVelocityAutocorrelationFunction(traj,atoms,groups[0],
                               refGroup[0],normVect=normVect,
                               timeInfo=timeInfo,diffScheme=diffScheme)
        saveText(data=N.transpose(N.array(AVACF)[:data.time_steps]),
                 filename=output,title=title)
    elif optlist.has_key('--methyl_aligned_avacf'):
        try:
            output = output['methyl_aligned_avacf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        AVACF = MethylAlignedAngularVelocityAutocorrelationFunction(traj,atoms,groups[0],
                               refGroup[0],normVect=normVect,
                               timeInfo=timeInfo,diffScheme=diffScheme)
        saveText(data=N.transpose(N.array(AVACF)[:data.time_steps]),
                 filename=output,title=title)
    elif optlist.has_key('--methyl_axis_avacf'):
        try:
            output = output['methyl_axis_avacf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        AVACF = MethylAxisAngularVelocityAutocorrelationFunction(traj,atoms,groups[0],
                               refGroup[0],normVect=normVect,
                               timeInfo=timeInfo,diffScheme=diffScheme)
        saveText(data=N.transpose(N.array(AVACF)[:data.time_steps]),
                 filename=output,title=title)
        
    elif optlist.has_key('--savacf'):
        try:
            output = output['savacf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        SAVACF = SpectraOfAngularVACF(traj,normVect=normVect,
                       diffScheme=diffScheme,groups=groups[0],
                       timeInfo=timeInfo,refGroup=refGroup[0],alpha=alpha)
        sample = max(1, len(SAVACF[0])/data.frequency_points)
        scaled = []
        scaled.append(N.array(SAVACF[0])[::sample]/uF)
        scaled.append(N.array(SAVACF[1])[::sample])
        saveText(data=N.transpose(N.array(scaled)),filename=output,title=title)

        ASCII_To_NetCDF(output, \
                        os.path.join(tempfile.gettempdir(),'ADOS_Reference.nc'), \
                        varname = ["frequency", "ados-total"], \
                        skiprows = 3, bycolumn = True)
        
    elif optlist.has_key('--rcf'):
        try:
            output = output['rcf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        rcf = ReorientationalCorrelationFunction(traj,jmn=rotCoef,
                       groups=groups[0],refGroup=refGroup[0],timeInfo=timeInfo)
        saveText(data=N.transpose(N.array(rcf))[:data.time_steps],
                 filename=output,title=title)

        ASCII_To_NetCDF(output, \
                        os.path.join(tempfile.gettempdir(),'RCF_Reference.nc'), \
                        varname = ["time", "rcf-total"], \
                        skiprows = 3, bycolumn = True)
        
    elif optlist.has_key('--df'):
        try:
            output = output['trajectory']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        gj = DigitalFilter(traj,output,groups=groups[0],
                           timeInfo=timeInfo,filterSet=setFilter)
    elif optlist.has_key('--rbt'):
        try:
            output = output['trajectory']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        gj = RigidBodyTrajectory(traj,output,groups[0],refGroup[0],
                                 timeInfo=timeInfo,remove_translation=0)
    elif optlist.has_key('--rbrt'):
        try:
            output = output['trajectory']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        gj = RigidBodyTrajectory(traj,output,groups[0],refGroup[0],
                                 timeInfo=timeInfo,remove_translation=1)
    elif optlist.has_key('--at'):
        try:
            output = output['trajectory']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        gj = AngularTrajectory(traj,output,groups[0],refGroup[0],
                               timeInfo=timeInfo)
    elif optlist.has_key('--mpcf'):
        try:
            output = output['mpcf']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        if len(groups) == 2:
            groups_A = groups[0];  refGroup_A = refGroup[0]
            groups_B = groups[1];  refGroup_B = refGroup[1]
        else:
            groups_A = groups[0];  refGroup_A = refGroup[0]
            groups_B = None;  refGroup_B = None
        bins, mpcf = MolecularPairCorrelationFunction(traj,groups_A,refGroup_A,
                     groups_B=groups_B,refGroup_B=refGroup_B,
                     symbols=symbols,timeInfo=timeInfo)
        saveText(data=N.transpose(N.array([bins/uL,mpcf])),
                 filename=output,title=title)
    elif optlist.has_key('--fft'):
        try:
            output = output['fft']
        except KeyError:
            sys.stderr.write('No output file specified.\n')
            raise SystemExit
        gj = ScatteringFunctionFFT(scatFunc,output,data.frequency_points,
                                   alpha=alpha)
    elif optlist.has_key('--progress'): getProgress()
else:
    print """
 Usage: pMoldyn calc_mode commands
 where 

       calc_mode can be one from the list below:

       --dos-vel             : Density Of States from velocities
       --dos-xyz             : Density Of States from coordinates
       --vacf-vel            : Velocity Autocorrelation Function 
                               from velocities
       --vacf-xyz            : Velocity Autocorrelation Function 
                               from coordinates
       --ar-vel              : Autoregressive analysis from velocity
       --ar-xyz              : Autoregressive analysis from coordinates
       --msd                 : Mean Square Displacement
       --avacf               : Angular Velocity Autocorrelation Function
       --savacf              : Spectra of Angular VACF
       --rcf                 : Rotational Correlation Function
       --df                  : Digital Filter
       --rbt                 : Rigid Body Trajectory
       --rbt                 : Rigid Body Rotation Trajectory
       --at                  : Angular Trajectory
       --csf                 : Coherent Scattering Function
       --arcsf               : Autoregressive analysis of coherent scattering
       --isf                 : Incoherent Scattering Function
       --isfg                : Incoherent Scattering Function 
                               in Gaussian Approximation
       --arisf               : Autoregressive analysis of incoherent scattering
       --fft                 : FFT of Scattering Functions
       --eisf                : Elastic Incoherent Structure Factor
 
       commands which influence calculations to be performed
       can be defined either in an external file:

       --input filename      : all settings are read from an input file

       To check a status of running calculations use only:

       --progress
"""

end_time = time()
print " Finished at:",asctime(localtime(end_time))
print " Wall clock time used by this task was: %.1f sec." % \
      (end_time-start_time)








