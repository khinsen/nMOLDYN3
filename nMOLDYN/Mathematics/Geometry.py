"""This modules contains functions related to general mathematics (analysis, geometry, algebra ...).

Functions:
    * sphericalCoordinates: converts from x, y and z cartesian coordinates to r, theta, phi spherical coordinates.
    * basisVectors        : computes the basis vectors of the simulation cell from a set of values defining its geometry (3 distances and 3 angles).
    * randomPointInCircle : returns a vector within a circle of radius |r| and orthogonal to a given direction.
    * randomDirection2D   : returns a normalized vector generated from a unit circle orthogonal to a given direction.
    * randomPlane2D       : generates a normalized random q-vector on a plane defined by vect1, vect2.
"""

# The python distribution modules
import copy
from random import uniform

# The ScientificPython modules
from Scientific import N
from Scientific.Geometry import Tensor, Vector

# The MMTK distribution modules
from MMTK.Random import randomDirection

# The nMOLDYN modules
from nMOLDYN.Core.Error import Error
from nMOLDYN.Mathematics.Analysis import differentiate
    
def sphericalCoordinates(x, y, z):
    """Returns the r, theta and phi spherical coordinates corresponding to x, y z cartesian coordinates.

    @param x: the cartesian x.
    @type x: float   

    @param y: the cartesian y.
    @type y: float   

    @param z: the cartesian z.
    @type z: float   

    @return: the r, theta and phi spherical coordinates..
    @rtype: a list of three floats
    """

    # The spherical radius is computed
    r = N.sqrt(x**2 + y**2 + z**2)

    # The spherical theta is computed
    theta = N.arccos(z/r)

    # The spherical phi is computed
    phi = N.arctan2(y,x)

    return r, theta, phi

def changeBasis(pt, op, ip, jp, kp):
    """Returns the coordinates of a point after a change of basis.
    
    @param pt: the coordinates of the point in the old basis.
    @type pt: Scientific Vector   

    @param op: the coordinates of the new origin in the old basis.
    @type op: Scientific Vector   

    @param ip: the coordinates of the new x axis in the old basis.
    @type ip: Scientific Vector   

    @param jp: the coordinates of the new y axis in the old basis.
    @type jp: Scientific Vector   

    @param kp: the coordinates of the new z axis in the old basis.
    @type kp: Scientific Vector   

    @return: the coordinates of the point in the new basis.
    @rtype: Scientific Vector   
    """
    
    mat = Tensor(N.array([ip,jp,kp], typecode = N.Float).transpose())
    matinv = mat.inverse()

    # This is the formula of the change of basis.
    mprim = matinv * (pt - op)

    return mprim

def basisVectors(parameters):
    """Returns the basis vectors for the simulation cell from the six crystallographic parameters.
    
    @param parameters: the a, b, c, alpha, bete and gamma of the simulation cell.
    @type: parameters: list of 6 floats
    
    @return: a list of three Scientific.Geometry.Vector objects representing respectively a, b and c 
        basis vectors.
    @rtype: list
    """

    # The simulation cell parameters.
    a, b, c, alpha, beta, gamma = parameters
    
    # By construction the a vector is aligned with the x axis.
    e1 = Vector(a, 0.0, 0.0)

    # By construction the b vector is in the xy plane.
    e2 = b*Vector(N.cos(gamma), N.sin(gamma), 0.0)
    
    e3_x = N.cos(beta)
    e3_y = (N.cos(alpha) - N.cos(beta)*N.cos(gamma)) / N.sin(gamma)
    e3_z = N.sqrt(1.0 - e3_x**2 - e3_y**2)
    e3 = c*Vector(e3_x, e3_y, e3_z)

    return (e1, e2, e3)

def randomPointInCircle(r,axis):
    """Returns a vector drawn from an uniform distribution within a circle of radius |r| and 
    orthogonal to vector |axis|.
    
    @param r: the radius of the circle.
    @type r: float
    
    @param axis: the axis orthogonal to the plane where the vectors have to be generated.
    @type axis: Scientific.Geometry.Vector

    @return: a vector pointing to a random point of the circle.
    @rtype: Scientific.Geometry.Vector
    """
    
    rsq = r*r
    while 1:
        x = N.array([uniform(-r, r), uniform(-r,r), uniform(-r,r)], typecode = N.Float)
        y = Vector(x) - (Vector(x)*axis/axis.length()**2)*axis
        if y*y < rsq :
            break
    return y

def randomDirection2D(axis):
    """Returns a normalized vector drawn from an uniform distribution on the surface of a unit circle on a
    plane orthogonal to |axis|.

    @param axis: the axis orthogonal to the plane where the vectors have to be generated.
    @type axis: Scientific.Geometry.Vector object
    
    @return: A normalized vector defined in a unit disk orthogonal to |axis|
    @rtype: Scientific.Geometry.Vector object
    """
    
    r = randomPointInCircle(1.0,axis)
    return r.normal()

def randomVector(directions = None):
    """Returns a normalized random vector on an axis, a plane or in space.
    
    @param directions: if not None, a list of 2 Scientific.Vector that will define the plane on which
        the vector should be generated.
    @type directions: list of 2 Scientific.Vector or None

    @return: a normalized random vector on a plane defined by |directions| or in space (|directions| = None).
    @rtype: Scientific.Geometry.Vector object
    """
    
    # 3D case.
    if directions is None:
        return randomDirection()
    
    # 1D case.
    if len(directions) == 1:
        return directions[0]

    # 2D case.
    elif len(directions) == 2:
        # The first vector of the plane.
        v1, v2 = directions

        # v1 /\ v2.
        axis = (v1.cross(v2))
    
        if axis.length() == 0.0:
            raise Error('Your basis vectors are colinear. Can not build a plane out of them.')

        return randomDirection2D(axis)
    
def qMatrix(data):
    
    res = N.zeros((len(data),4,4), typecode = N.Float)
    # qs
    res[:,0,0]=  data[:,0] 
    res[:,1,1]=  data[:,0]
    res[:,2,2]=  data[:,0]
    res[:,3,3]=  data[:,0]
    # qx
    res[:,0,1]=  data[:,1]
    res[:,1,0]= -data[:,1]
    res[:,2,3]=  data[:,1]
    res[:,3,2]= -data[:,1]
    # qy
    res[:,0,2]=  data[:,2]
    res[:,1,3]= -data[:,2]
    res[:,2,0]= -data[:,2]
    res[:,3,1]=  data[:,2]
    # qz
    res[:,0,3]=  data[:,3]
    res[:,1,2]=  data[:,3]
    res[:,2,1]= -data[:,3]
    res[:,3,0]= -data[:,3]
        
    return res

def getAngularVelocity(trajectory, group, frameIndexes, dt, stepwise = False, referenceFrame = 0, differentiation = 1):
    """Computes the Angular Velocity Function."""
    
    nFrames = len(frameIndexes)
    
    if nFrames <= 1:
        raise Error('There must be at least two steps.')
            
    quaternions = N.zeros((nFrames, 4), typecode = N.Float)

    # Case of a moving reference.
    if stepwise:
            
        # The reference configuration is always the one of the previous frame excepted for the first frame
        # where it is set by definition to the first frame (could we think about a cyclic alternative way ?).
        for comp in range(nFrames):
                
            frameIndex = frameIndexes[comp]
                
            if comp == 0:
                previousFrame = frameIndexes[0]
                    
            else:
                previousFrame = frameIndexes[comp-1]
                    
            refConfig = trajectory.configuration[previousFrame]

            # The RBT is created just for the current step.
            rbt = trajectory.readRigidBodyTrajectory(group,\
                                                     first = frameIndex,\
                                                     last = frameIndex + 1,\
                                                     skip = 1,\
                                                     reference = refConfig)
                
            # The corresponding quaternions and cms are stored in their corresponding matrix.
            quaternions[comp,:] = copy.copy(rbt.quaternions)

    # The simplest case, the reference frame is fixed.
    # A unique RBT is performed from first to last skipping skip steps and using refConfig as the reference.
    else:
            
        # If a fixed reference has been set. We can already set the reference configuration here.
        refConfig = trajectory.configuration[referenceFrame]
        
        first = frameIndexes[0]
        last = frameIndexes[-1] + 1
        skip = frameIndexes[1] - frameIndexes[0]
        
        # The RBT is created.
        rbt = trajectory.readRigidBodyTrajectory(group, first = first, last = last, skip = skip, reference = refConfig)
            
        quaternions[:,:] = copy.copy(rbt.quaternions)

    # The quaternions derivatives.
    quaternions_dot = N.zeros((nFrames,4), typecode = N.Float)
        
    for i in range(4):
        quaternions_dot[:,i] = differentiate(quaternions[:,i], differentiation, dt)

    q = qMatrix(quaternions)

    angvel = 2.0*N.add.reduce(q*quaternions_dot[:,N.NewAxis,:], -1)
    
    # 1: refers to the vectorial part of the quaternion matrix.
    return angvel[:,1:]
