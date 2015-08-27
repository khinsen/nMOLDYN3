# This gives the correspondance between the old 'percentage of trajectory length' parameter
# and the new 'energy fwhm' parameter.

from Scientific import N

dt = N.array([1, 4, 2, 2, 2, 2, 2, 2, 2, 3])

time_step = N.array([0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01, 0.01])

n_frames = N.array([19, 13, 11, 10, 10, 10, 10, 20, 9, 15])

per_traj_length = N.array([4, 10, 25, 20, 40, 8, 20, 20, 20, 20])

fwhm_e = 2.0*N.sqrt(2.0*N.log(2.0))*100.0/(1.5192669*dt*time_step*per_traj_length*(n_frames - 1.0))

print fwhm_e
