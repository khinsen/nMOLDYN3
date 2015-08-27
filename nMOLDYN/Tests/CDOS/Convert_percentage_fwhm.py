# This gives the correspondance between the old 'percentage of trajectory length' parameter
# and the new 'energy fwhm' parameter.

from Scientific import N

dt = N.array([0.005,0.02,0.01,0.02,0.01,0.01,0.01,0.08,0.04,0.08,0.08,0.08,0.08,0.08])

n_frames = N.array([49,13,11,6,10,10,10,10,19,10,10,10,10,20])

per_traj_length = N.array([20.0, 10.0, 5.0, 20.0, 4.0, 40.0, 20.0, 80.0, 4.0, 20.0, 10.0, 20.0, 20.0, 12.0])

fwhm_e = 2.0*N.sqrt(2.0*N.log(2.0))*100.0/(1.5192669*dt*per_traj_length*(n_frames - 1.0))

print fwhm_e
