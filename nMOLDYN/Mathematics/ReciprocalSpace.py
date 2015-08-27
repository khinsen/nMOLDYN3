"""This modules implements the class that generates a set of Q-vectors within a set of q shells.

Classes:
    * QVectors: the class that actually performs the q vectors generation.
"""

# The python distribution modules
import copy
from random import sample, shuffle, uniform

import numpy

# The Scientific modules
from Scientific import N
from Scientific.Geometry import Tensor, Vector

# The nMOLDYN modules
from nMOLDYN.Core.Error import Error
from nMOLDYN.Core.Logger import LogMessage
from nMOLDYN.Mathematics.Analysis import pgcd
from nMOLDYN.Mathematics.Geometry import randomVector

def Generate_QShells(qShells):
    
    genQShells = []
            
    if isinstance(qShells, (list, tuple)):
        
        genQShells = qShells
        
    elif isinstance(qShells, str):
        
        for qInter in [v.strip() for v in qShells.split(";")]:
            
            qMin, qMax, dQ = [float(v) for v in qInter.split(":")]
            
            temp = N.arange(qMin, qMax + 0.001*dQ, dQ).tolist()
            
            if temp[-1] > qMax:
                del temp[-1]
                
            genQShells.extend(temp)
                    
    else:
        raise Error("%s: invalid format for qshells." % genQShells)
    
    genQShells = [v for v in sorted(set(genQShells)) if v > 0.0]
    
    if len(genQShells) == 0:
        raise Error("%s triggered emptry qshells generation." % genQShells)
    
    return genQShells
        
class QVectors(object):
    
    def __init__(self, universe, qVectorsDict):
        
        if not isinstance(qVectorsDict, dict):
            raise Error("%s: not a valid python dictionnary." % qVectorsDict)
        
        self.qVectorsDict = copy.copy(qVectorsDict)

        # A copy of the input universe.
        self.universe = universe
                
        if self.universe.reciprocalBasisVectors() is None:
            
            # The direct basis.
            self.dirBasis = None

            # The reciprocal basis.
            self.recBasis = None
            
        else:
            
            # The direct basis.
            self.dirBasis = Tensor(N.array([N.array(v) for v in self.universe.basisVectors()], typecode = N.Float))/(2.0*N.pi)

            # The reciprocal basis.
            self.recBasis = (2.0*N.pi)*Tensor(N.array([N.array(v) for v in self.universe.reciprocalBasisVectors()], typecode = N.Float))

            self.transRecBasis = self.recBasis.transpose()
             
        # These tree qttributes will the "output" for the class.
        self.qRadii = []
        self.qVectors = []
        self.hkls = []                

        self.qGeometry = self.qVectorsDict.get("qgeometry", None) 
    
        if self.qGeometry in ["spatial", "planar", "axial"]:

            qShells = self.qVectorsDict.get("qshellvalues", None)        
            self.qShells = Generate_QShells(qShells)

            try:
                self.qShellWidth = float(self.qVectorsDict.get("qshellwidth"))
            except:
                raise Error("Invalid value for qshellwidth parameter.")
        
            try:
                self.qVectorsPerShell = int(self.qVectorsDict.get("qvectorspershell"))
            except:
                raise Error("Invalid value for qvectorspershell parameter.")

            self.qGeometry = self.qVectorsDict.get("qgeometry", None) 
        
            if self.qGeometry == "spatial":
                self._spatial_qvectors()
            
            elif self.qGeometry == "planar":
                self._planar_qvectors()
            
            elif self.qGeometry == "axial":
                self._axial_qvectors()
            
            else:
                raise Error("%s is not a valid Qvector geometry." % self.qGeometry)
        
        elif self.qGeometry == "userdefined":
            self._userdefined_qvectors()
        
        else:
            raise Error("%s is not a valid Qvector geometry." % self.qGeometry)

        self.qRadii = N.array(self.qRadii, typecode = N.Float)

        self.qvectors_statistics()

        self._display_qvectors()

    def _explicit_qvectors(self, qMin, qMax, indRange):
        
        qVects = []
        hkls = []
        
        dimen = len(self.qDirections)

        for ind in indRange:

            qVect = Vector()            
            for i in range(dimen):
                qVect += ind[i]*self.qDirections[i]
        
            if qMin < qVect.length() <= qMax:
                
                if not qVect in qVects:
                    
                    qVects.append(qVect)
                    
                    hkl = Vector([0,0,0])
                    for i in range(dimen):
                        hkl += ind[i]*self.hkl_dir[i]
                        
                    hkls.append(hkl)                                        

                    if len(qVects) >= self.qVectorsPerShell:
                        break

        LogMessage('info', '%d explicit Q vectors generated for shell [%s,%s].' % (len(qVects), qMin, qMax), ['file','console'])
        
        return qVects, hkls

    def _random_qvectors(self, qMin, qMax):

        nRedundant = 0

        qVects = []

        hkls = []

        while(len(qVects) < self.qVectorsPerShell):
            
            # A random vector is generated.
            qVect = randomVector(None)*uniform(qMin, qMax)

            # Case of a periodic universe. The q vector is rounded to the closest q vectors pointing on a reciprocal node.
            if self.recBasis is not None:

                # This is a shortcut to get hkl from h=q.a, k=q.b, l=q.c.
                hkl = Vector([round(v) for v in N.array(self.dirBasis*qVect, typecode = N.Float)])
                
                # The approximated reciprocal space vector.
                qVect = self.transRecBasis*hkl

            else:

                hkl = qVect
                
            # Checks whether the q vector falls inside the q shell.
            if qMin < qVect.length() <= qMax:
                
                # If the q vector is already stored in |generatedQVectors| list skip it.
                if qVect in qVects:
                        
                    # Increment the redundancy counter.
                    nRedundant += 1
                        
                    # After 5000 trials to generate a new q vector, gives up and returns 
                    # the generated q vectors list in its current state.
                    if nRedundant == 5000:
                        
                        # Displays on the console and file loggers that no more q vectors could be generated for that shell.
                        LogMessage('warning', 'Give up to generate a new q vector for shell with radius [%s,%s] after 5000 trials.' % (qMin, qMax), ['file','console'])
                        break

                # Otherwise append it to the |qVects| and |hkls| list.
                else:
                    
                    qVects.append(qVect)
                    hkls.append(hkl)
                        
                    # A new q vector could be found, so reset the redundancy counter to 0.
                    nRedundant = 0
                    
            # If not, retries a new q vector generation
            else:
                
                # Increment the redundancy counter.
                nRedundant += 1
                
                # After 5000 trials to generate a new q vector, gives up and returns the q vector list in its current state.
                if nRedundant == 5000:

                    # Displays on the console and file loggers that no more q vectors could be generated for that shell.
                    LogMessage('warning', 'Give up to generate a new q vector for shell with radius [%s,%s] after 5000 trials.' % (qMin, qMax), ['file','console'])
                    break

        LogMessage('info', '%d random Q vectors generated for shell [%s,%s].' % (len(qVects), qMin, qMax), ['file','console'])
        
        return qVects, hkls
    
    def _spatial_qvectors(self):

        # The reciprocal cell volume.
        recCellVolume = Vector(self.recBasis[0])*(Vector(self.recBasis[1]).cross(Vector(self.recBasis[2])))
        
        for q in self.qShells:

            qMin = max(q - 0.5*self.qShellWidth, 0.0)

            qMax = q + 0.5*self.qShellWidth

            if self.recBasis is None:
                
                qVects, hkls = self._random_qvectors(qMin, qMax)
            
            else:
                                
                self.hkl_dir = [Vector([1,0,0]), Vector([0,1,0]), Vector([0,0,1])]

                # The volume of the q shell.
                recShellVolume = 4.0*N.pi*q**2 * (qMax - qMin)
            
                # Estimation of the number of elementary cells contained in the q shell as the ratio
                # between the volume of the q shell and the volume of the elementary cell.
                nCells = recShellVolume/recCellVolume

                # If the number of estimated cells is > 500, then proceed to a random generation.
                if nCells > 500:
                    
                    # Displayss that it will be a random geneation.
                    qVects, hkls = self._random_qvectors(qMin, qMax)
                
                # Otherwise, the number is considered small enough for an explicit q vectors generation.
                else:

                    self.qDirections = [Vector(v) for v in self.recBasis]

                    indMax = [int(qMax/v.length()) + 2 for v in self.qDirections]

                    grid = N.transpose(numpy.mgrid[-indMax[0]:indMax[0],-indMax[1]:indMax[1],-indMax[2]:indMax[2]])

                    indRange = grid.reshape(grid.size/3,3).tolist()
                    
                    shuffle(indRange)
                    
                    # Displays that it will be an exhaustive generation.
                    qVects, hkls = self._explicit_qvectors(qMin, qMax, indRange)

            if qVects:
                
                self.qRadii.append(q)
                self.qVectors.append(qVects)
                self.hkls.append(hkls)
                
    def _planar_qvectors(self):
                        
        hkl_dir = self._generate_hkls()

        self.hkl_dir = []
    
        factor = pgcd(hkl_dir[0])
        self.hkl_dir.append(Vector(hkl_dir[0])/factor)

        factor = pgcd(hkl_dir[1])
        self.hkl_dir.append(Vector(hkl_dir[1])/factor)

        self.qDirections = []
        
        self.qDirections.append(self.transRecBasis*self.hkl_dir[0])
        self.qDirections.append(self.transRecBasis*self.hkl_dir[1])
        
        if  self.qDirections[0].cross(self.qDirections[1]).length() <= 0:
            raise Error("The q plan must be defined by two non-colinear q-directions.")

        for q in self.qShells:
            
            qMin = max(q - 0.5*self.qShellWidth, 0.0)

            qMax = q + 0.5*self.qShellWidth

            indMax = [int(qMax/v.length()) + 2 for v in self.qDirections]

            grid = N.transpose(numpy.mgrid[-indMax[0]:indMax[0],-indMax[1]:indMax[1]])

            indRange = (-grid).reshape(grid.size/2,2).tolist() + grid.reshape(grid.size/2,2).tolist()

            shuffle(indRange)

            qVects, hkls = self._explicit_qvectors(qMin, qMax, indRange)

            if qVects:
                
                self.qRadii.append(q)
                self.qVectors.append(qVects)
                self.hkls.append(hkls)
        
    def _axial_qvectors(self):
                        
        hkl_dir = self._generate_hkls()

        factor = pgcd(hkl_dir[0])

        self.hkl_dir = [Vector(hkl_dir[0])/factor]

        self.qDirections = [self.transRecBasis*self.hkl_dir[0]]
        
        qLength = self.qDirections[0].length()

        if qLength <= 0.0:
            raise Error("The %s node gave a null direction in the reciprocal space." % self.qDirection)

        for q in self.qShells:
            
            qMin = max(q - 0.5*self.qShellWidth, 0.0)

            if qLength > qMin:
                continue

            qMax = q + 0.5*self.qShellWidth

            indMin = int(qMin/qLength)
            indMax = int(qMax/qLength) + 1

            indRange = range(-indMax, -indMin + 1) + range(indMin, indMax + 1)

            shuffle(indRange)
            indRange = [[v] for v in indRange]

            qVects, hkls = self._explicit_qvectors(qMin, qMax, indRange)

            if qVects:
                
                self.qRadii.append(q)                  
                self.qVectors.append(qVects)
                self.hkls.append(hkls)

    def _userdefined_qvectors(self):
        
        hkls = self._generate_hkls()
                
        if self.universe.reciprocalBasisVectors() is None:

            qVects = [Vector(v) for v in hkls]

        else:

            qVects = [self.transRecBasis*Vector(hkl) for hkl in hkls]

        qVectsDict = {}

        for i in range(len(qVects)):

            qVect = qVects[i]
            hkl = Vector(hkls[i])

            qLength = round(qVect.length(),3)

            if qVectsDict.has_key(qLength):
                qVectsDict[qLength]["qvect"].append(qVect)
                qVectsDict[qLength]["hkl"].append(hkl)
                
            else:
                qVectsDict[qLength] = {"qvect" : [qVect], "hkl" : [hkl]}

        for q in sorted(qVectsDict.keys()):
            self.qRadii.append(q)
            self.qVectors.append(qVectsDict[q]["qvect"])
            self.hkls.append(qVectsDict[q]["hkl"])
                        
    def _generate_hkls(self):
        
        if self.recBasis is None:
            raise Error("The universe must be periodic for an anisotropic q vectors generation.")

        hkls = self.qVectorsDict.get("hkls", None)
        
        if hkls is None:
            raise Error("No hkls provided for the qvectors generation.")
        
        if self.qGeometry == 'planar':
            
            try:

                if isinstance(hkls, (list, tuple)):
                    if len(hkls) != 2:
                        raise

                elif isinstance(hkls, str):
                    hkls = [[int(vv.strip()) for vv in v.strip().split(",")] for v in hkls.split(";")]

                    if len(hkls) != 2:
                        raise
        
                else:
                    raise 
            
            except:
                raise Error("Invalid format for hkls")
        
            else:
                return hkls                    
                
        elif self.qGeometry == 'axial':
                
            try:
                if isinstance(hkls, (list, tuple)):
                    if len(hkls) != 1:
                        raise

                elif isinstance(hkls, str):
                    hkls = [[int(v.strip()) for v in hkls.split(",")]]
                    if len(hkls) != 1:
                        raise
                
                else:
                    raise
            
            except:
                raise Error("Invalid format for hkls")
        
            else:
                return hkls
                
        elif self.qGeometry == 'userdefined':
            
            if isinstance(hkls, (list,tuple)):
            
                return hkls
        
            elif isinstance(hkls, str):
            
                hs, ks, ls = [[int(vv.strip()) for vv in v.strip().split(":")] for v in hkls.split(",")]            

                grid = N.transpose(numpy.mgrid[hs[0]:hs[1]+1:hs[2],ks[0]:ks[1]+1:ks[2],ls[0]:ls[1]+1:ls[2]])
            
                hkls = grid.reshape(grid.size/3,3).tolist()
                                                        
                return hkls

            else:
                raise Error("Invalid format for hkls")
                
    def qvectors_statistics(self):
        
        # The |self.statistics| attribute is set to a 2D array that will store the number of generated q vectors 
        # per q shell and per q space octan. This may be useful to check for the homogeneity of the q vectors generation 
        # as it is usually performed randomy. 
        self.statistics = N.zeros((len(self.qVectors), 8), typecode = N.Int32)        
        # Loop over the q shell radii indexes.
        for comp in range(len(self.qRadii)):
            # The generated list of q vectors for that q shell.
            qVectors = self.qVectors[comp]
            
            # Loop over each q vector of the generated  list.
            for qVect in qVectors:
                
                ind = 4*(qVect[0] > 0) + 2*(qVect[1] > 0) + qVect[2] > 0
                
                self.statistics[comp, ind] += 1
                
    def _display_qvectors(self):
        LogMessage("info", "Q shells contents:", ["file"])
        
        for comp in range(len(self.qRadii)):
            LogMessage("info", "Q = %8.3f  h   k   l          qx       qy       qz" % self.qRadii[comp], ["file"])
            for comp1 in range(len(self.qVectors[comp])):
                LogMessage("info", "            %3d %3d %3d    %8.3f %8.3f %8.3f" % tuple(list(self.hkls[comp][comp1]) + list(self.qVectors[comp][comp1])), ["file"])
            LogMessage("info", "", ["file"])
                
if __name__ == "__main__":
    
    from MMTK.Universe import CubicPeriodicUniverse
    
    u = CubicPeriodicUniverse(10)    
    qv = QVectors(u, {"qgeometry" : "planar", "qvectorspershell" : 5, "qshellwidth" : 1.0, "qshellvalues" : "0.0:10.0:1.0", "hkls" : "1,0,0;0,1,0"})
    print "Planar 1", qv.hkls
    
    qv = QVectors(u, {"qgeometry" : "planar", "qvectorspershell" : 5, "qshellwidth" : 1.0, "qshellvalues" : "0.0:10.0:1.0", "hkls" : [[1,0,0],[0,1,0]]})
    print "Planar 2", qv.hkls

    qv = QVectors(u, {"qgeometry" : "axial", "qvectorspershell" : 5, "qshellwidth" : 1.0, "qshellvalues" : "0.0:10.0:1.0", "hkls" : "1,0,0"})
    print "Axial 1", qv.hkls

    qv = QVectors(u, {"qgeometry" : "axial", "qvectorspershell" : 5, "qshellwidth" : 1.0, "qshellvalues" : "0.0:10.0:1.0", "hkls" : [[1,0,0]]})
    print "Axial 2", qv.hkls
    
    qv = QVectors(u, {"qgeometry" : "userdefined", "hkls" : "1:3:1,1:4:1,2:5:2"})
    print "Userdefined 1", qv.hkls

