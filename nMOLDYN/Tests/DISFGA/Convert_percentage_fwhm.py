# This gives the correspondance between the old 'percentage of trajectory length' parameter
# and the new 'energy fwhm' parameter.

from Scientific import N

dt = N.array([0.005,0.005,0.01,0.005,0.015,0.01,0.005,0.02,0.005,0.015,0.0200005,0.0100002])

n_frames = N.array([49,19,10,9,8,10,19,9,19,9,7,19])

per_traj_length = N.array([80.0,10.0,20.0,50.0,10.0,40.0,80.0,10.0,20.0,25.0,10.0,8.0])

fwhm_e = 2.0*N.sqrt(2.0*N.log(2.0))*100.0/(1.5192669*dt*per_traj_length*(n_frames - 1.0))

print fwhm_e