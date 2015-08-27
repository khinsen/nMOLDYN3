"""This modules stores the templates for all nMOLDYN analysis.
"""

# The nMOLDYN modules.
from nMOLDYN.Analysis.Dynamics import *
from nMOLDYN.Analysis.Master import ParallelPerAtom, ParallelPerFrame, ParallelPerGroup, ParallelPerQShell, SerialPerAtom, SerialPerFrame, SerialPerGroup, SerialPerQShell
from nMOLDYN.Analysis.NMR import *
from nMOLDYN.Analysis.Scattering import *
from nMOLDYN.Analysis.Structure import *
from nMOLDYN.Analysis.RigidBody import *
from nMOLDYN.Analysis.Trajectory import *
        
# The templates for each analysis that will be parallelized.
# By construction, those that are not currently parallelizable will have only a serial version.

# ############
# MSD
# ############
class MeanSquareDisplacement_serial(MeanSquareDisplacement, SerialPerAtom):
    pass

class MeanSquareDisplacement_parallel(MeanSquareDisplacement, ParallelPerAtom):
    pass

# ############
# RMSD
# ############
class RootMeanSquareDeviation_serial(RootMeanSquareDeviation, SerialPerAtom):
    pass

class RootMeanSquareDeviation_parallel(RootMeanSquareDeviation, ParallelPerAtom):
    pass

# ############
# VACF
# ############
class CartesianVelocityAutoCorrelationFunction_serial(CartesianVelocityAutoCorrelationFunction, SerialPerAtom):
    pass

class CartesianVelocityAutoCorrelationFunction_parallel(CartesianVelocityAutoCorrelationFunction, ParallelPerAtom):
    pass

# ############
# PACF
# ############
class CartesianPositionAutoCorrelationFunction_serial(CartesianPositionAutoCorrelationFunction, SerialPerAtom):
    pass

class CartesianPositionAutoCorrelationFunction_parallel(CartesianPositionAutoCorrelationFunction, ParallelPerAtom):
    pass

# ############
# DOS
# ############
class CartesianDensityOfStates_serial(CartesianDensityOfStates, SerialPerAtom):
    pass

class CartesianDensityOfStates_parallel(CartesianDensityOfStates, ParallelPerAtom):
    pass

# ############
# ARA
# ############
class AutoRegressiveAnalysis_serial(AutoRegressiveAnalysis, SerialPerAtom):
    pass

class AutoRegressiveAnalysis_parallel(AutoRegressiveAnalysis, ParallelPerAtom):
    pass

# ############
# PBFT
# ############
class PassBandFilteredTrajectory_serial(PassBandFilteredTrajectory, SerialPerAtom):
    pass

class PassBandFilteredTrajectory_parallel(PassBandFilteredTrajectory, ParallelPerAtom):
    pass

# ############
# DISF
# ############
class DynamicIncoherentStructureFactor_serial(DynamicIncoherentStructureFactor, SerialPerAtom):
    pass

class DynamicIncoherentStructureFactor_parallel(DynamicIncoherentStructureFactor, ParallelPerAtom):
    pass

# ############
# DISFGA
# ############
class DynamicIncoherentStructureFactorGaussianApproximation_serial(DynamicIncoherentStructureFactorGaussianApproximation, SerialPerAtom):
    pass

class DynamicIncoherentStructureFactorGaussianApproximation_parallel(DynamicIncoherentStructureFactorGaussianApproximation, ParallelPerAtom):
    pass

# ############
# EISF
# ############
class ElasticIncoherentStructureFactor_serial(ElasticIncoherentStructureFactor, SerialPerAtom):
    pass

class ElasticIncoherentStructureFactor_parallel(ElasticIncoherentStructureFactor, ParallelPerAtom):
    pass

# ############
# ROG
# ############
class RadiusOfGyration_serial(RadiusOfGyration, SerialPerAtom):
    pass

class RadiusOfGyration_parallel(RadiusOfGyration, ParallelPerAtom):
    pass

# ############
# ARDISF
# ############
class AutoRegressiveDynamicIncoherentStructureFactor_serial(AutoRegressiveDynamicIncoherentStructureFactor, SerialPerAtom):
    pass

class AutoRegressiveDynamicIncoherentStructureFactor_parallel(AutoRegressiveDynamicIncoherentStructureFactor, ParallelPerAtom):
    pass

# ############
# SSCSF
# ############
class SmoothedStaticCoherentStructureFactor_serial(SmoothedStaticCoherentStructureFactor, SerialPerFrame):
    pass

class SmoothedStaticCoherentStructureFactor_parallel(SmoothedStaticCoherentStructureFactor, ParallelPerFrame):
    pass

# ############
# HBSA
# ############
class HydrogenBondSurvivalAnalysis_serial(HydrogenBondSurvivalAnalysis, SerialPerFrame):
    pass

class HydrogenBondSurvivalAnalysis_parallel(HydrogenBondSurvivalAnalysis, ParallelPerFrame):
    pass

# ############
# PDF
# ############
class PairDistributionFunction_serial(PairDistributionFunction, SerialPerFrame):
    pass

class PairDistributionFunction_parallel(PairDistributionFunction, ParallelPerFrame):
    pass

# ############
# CN
# ############
class CoordinationNumber_serial(CoordinationNumber, SerialPerFrame):
    pass

class CoordinationNumber_parallel(CoordinationNumber, ParallelPerFrame):
    pass

# ############
# DP
# ############
class DensityProfile_serial(DensityProfile, SerialPerFrame):
    pass

class DensityProfile_parallel(DensityProfile, ParallelPerFrame):
    pass

# ############
# SFA
# ############
class ScrewFitAnalysis_serial(ScrewFitAnalysis, SerialPerFrame):
    pass

class ScrewFitAnalysis_parallel(ScrewFitAnalysis, ParallelPerFrame):
    pass

# ############
# SD
# ############
class SpatialDensity_serial(SpatialDensity, SerialPerFrame):
    pass

class SpatialDensity_parallel(SpatialDensity, ParallelPerFrame):
    pass

# ############
# GMFT
# ############
class GlobalMotionFilteredTrajectory_serial(GlobalMotionFilteredTrajectory, SerialPerFrame):
    pass

class GlobalMotionFilteredTrajectory_parallel(GlobalMotionFilteredTrajectory, ParallelPerFrame):
    pass

# ############
# RT
# ############
class ReducedTrajectory_serial(ReducedTrajectory, SerialPerFrame):
    pass

#class ReducedTrajectory_parallel(ReducedTrajectory, ParallelPerFrame):
#    pass

# ############
# COMT
# ############
class CenterOfMassTrajectory_serial(CenterOfMassTrajectory, SerialPerFrame):
    pass

class CenterOfMassTrajectory_parallel(CenterOfMassTrajectory, ParallelPerFrame):
    pass

# ############
# OPCM
# ############
class OrderParameterContactModel_serial(OrderParameterContactModel, SerialPerFrame):
    pass

class OrderParameterContactModel_parallel(OrderParameterContactModel, ParallelPerFrame):
    pass

# ############
# AC
# ############
class AngularCorrelation_serial(AngularCorrelation, SerialPerGroup):
    pass

class AngularCorrelation_parallel(AngularCorrelation, ParallelPerGroup):
    pass

# ############
# RBT
# ############
class RigidBodyTrajectory_serial(RigidBodyTrajectory, SerialPerGroup):
    pass

class RigidBodyTrajectory_parallel(RigidBodyTrajectory, ParallelPerGroup):
    pass

# ############
# RCF
# ############
class ReorientationalCorrelationFunction_serial(ReorientationalCorrelationFunction, SerialPerGroup):
    pass

class ReorientationalCorrelationFunction_parallel(ReorientationalCorrelationFunction, ParallelPerGroup):
    pass

# ############
# AVACF
# ############
class AngularVelocityAutoCorrelationFunction_serial(AngularVelocityAutoCorrelationFunction, SerialPerGroup):
    pass

class AngularVelocityAutoCorrelationFunction_parallel(AngularVelocityAutoCorrelationFunction, ParallelPerGroup):
    pass

# ############
# ADOS
# ############
class AngularDensityOfStates_serial(AngularDensityOfStates, SerialPerGroup):
    pass

class AngularDensityOfStates_parallel(AngularDensityOfStates, ParallelPerGroup):
    pass

# ############
# OP
# ############
class OrderParameter_serial(OrderParameter, SerialPerGroup):
    pass

class OrderParameter_parallel(OrderParameter, ParallelPerGroup):
    pass

# ############
# QHA
# ############
class QuasiHarmonicAnalysis_serial(QuasiHarmonicAnalysis):
    pass

# ############
# DCSF
# ############
class DynamicCoherentStructureFactor_serial(DynamicCoherentStructureFactor, SerialPerQShell):
    pass

class DynamicCoherentStructureFactor_parallel(DynamicCoherentStructureFactor, ParallelPerQShell):
    pass

# ############
# SCSF
# ############
class StaticCoherentStructureFactor_serial(StaticCoherentStructureFactor, SerialPerQShell):
    pass

class StaticCoherentStructureFactor_parallel(StaticCoherentStructureFactor, ParallelPerQShell):
    pass

# ############
# ARDCSFAR
# ############
class AutoRegressiveDynamicCoherentStructureFactor_serial(AutoRegressiveDynamicCoherentStructureFactor, SerialPerQShell):
    pass

class AutoRegressiveDynamicCoherentStructureFactor_parallel(AutoRegressiveDynamicCoherentStructureFactor, ParallelPerQShell):
    pass
