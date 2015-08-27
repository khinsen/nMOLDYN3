"""This modules contains the functions used by Pyro slave to perform analysis remotely.

Functions:
    * do_analysisPerElement: performs an analysis element-by-element.
"""

from MMTK.Trajectory import Trajectory
from nMOLDYN.Chemistry.Chemistry import hierarchizeUniverse

trajectory = None

# Define (or import) all the task handlers.
def do_analysisPerElement(analysis, element, trajname):
    """Performs the analysis element-by-element, the element being either
    an atom (atom-by-atom analysis), a frame index (frame-by-frame analysis),
    a group of atom (group-by-group analysis) or a set of q vectors.
    
    @param analysis: the selected analysis.
    @type analysis: a subclass of nMOLDYN.Analysis.Analysis.Analysis class
    
    @param element: the element on which the analysis is based.
    @type element: MMTK.Atom|integer|MMTK.Collections.Collection|nMOLDYN.Mathematics.QVectors
    
    @param trajname: a string specifying the name of the trajectory.
    @type trajname: string
    
    @return: the results of the analysis performed on one element.
    @rtype: depends on the analysis    
    """

    global trajectory
    if trajectory is None:
        trajectory = Trajectory(None, trajname)
        hierarchizeUniverse(trajectory.universe)
    return analysis.calc(element, trajectory)

if __name__ == '__main__':
    import sys
    from Scientific.DistributedComputing.MasterSlave import startSlaveProcess
    startSlaveProcess(sys.argv[1])

