"""
Functions for initializing the simulations for tracking.
"""

import numpy as np
import os
import netCDF4 as netCDF
import pdb
import glob
from datetime import datetime, timedelta
from matplotlib.mlab import *


def parameters():
	'''
	Parameters for running numerical simulation.
	loc 	Path to directory of grid and output files
	nsteps 	Number of steps to do between model outputs (iter in tracmass)
	ndays 	number of days to track the particles from start date
	ff 		ff=1 to go forward in time and ff=-1 for backward in time
	tseas	Time between outputs in seconds
	ah 		Horizontal diffusion in m^2/s. 
			See project values of 350, 100, 0, 2000. For -turb,-diffusion
	av 		Vertical diffusion in m^2/s.
	do3d 	for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
	doturb	turbulence/diffusion flag. 
			doturb=0 means no turb/diffusion,
			doturb=1 means adding parameterized turbulence
			doturb=2 means adding diffusion on a circle
			doturb=3 means adding diffusion on an ellipse (anisodiffusion)
	z0/zpar For 3D drifter movement, turn off twodim flag in makefile.
			Then z0 should be an array of initial drifter depths. 
			The array should be the same size as lon0 and be negative
			for under water. Currently drifter depths need to be above 
			the seabed for every x,y particle location for the script to run.
			To do 3D but start at surface, use z0=zeros(ia.shape) and have
			 either zpar='fromMSL'
			choose fromMSL to have z0 starting depths be for that depth below the base 
			time-independent sea level (or mean sea level).
			choose 'fromZeta' to have z0 starting depths be for that depth below the
			time-dependent sea surface. Haven't quite finished the 'fromZeta' case.
			For 2D drifter movement, turn on twodim flag in makefile.
			Then: 
			set z0 to 's' for 2D along a terrain-following slice
			 and zpar to be the index of s level you want to use (0 to km-1)
			set z0 to 'rho' for 2D along a density surface
			 and zpar to be the density value you want to use
			 Can do the same thing with salinity ('salt') or temperature ('temp')
			 The model output doesn't currently have density though.
			set z0 to 'z' for 2D along a depth slice
			 and zpar to be the constant (negative) depth value you want to use
			To simulate drifters at the surface, set z0 to 's' 
			 and zpar = grid['km']-1 to put them in the upper s level
			 z0='s' is currently not working correctly!!!
			 In the meantime, do surface using the 3d set up option but with 2d flag set
	xp 		x-locations in x,y coordinates for drifters
	yp 		y-locations in x,y coordinates for drifters
	zp 		z-locations (depths from mean sea level) for drifters
	t 		time for drifter tracks
	'''

	# Location of TXLA model output
	# file and then grid. 
	# 0150 file goes from (2009, 11, 19, 12, 0) to (2009, 12, 6, 0, 0)
	# loc = ['http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6/', \
	# 		'http://barataria.tamu.edu:8080//thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc']
	loc = ['http://barataria.tamu.edu:8080/thredds/dodsC/txla_nesting6/ocean_his_0150.nc', \
			'http://barataria.tamu.edu:8080//thredds/dodsC/txla_nesting6_grid/txla_grd_v4_new.nc']
	# # Location of TXLA model output
	# if 'rainier' in os.uname():
	# 	loc = '/Users/kthyng/Documents/research/postdoc/' # for model outputs
	# elif 'hafen.tamu.edu' in os.uname():
	# 	loc = '/home/kthyng/shelf/' # for model outputs

	# Initialize parameters
	nsteps = 5
	ndays = 1 #16
	ff = 1
	# date = datetime(2009,11, 20, 0)

	# Time between outputs
	# Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
	tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
	ah = 5. #100.
	av = 1.e-5 # m^2/s, or try 5e-6
=
	## Choose method for vertical placement of drifters
	# Also update makefile accordingly. Choose the twodim flag for isoslice.
	# See above for more notes, but do the following two lines for an isoslice
	z0 = 's'  #'salt' #'s' #'z' #'salt' #'s' 
	zpar = 29 #30 #29 #-10 #grid['km']-1 # 30 #grid['km']-1
	# Do the following two for a 3d simulation
	# z0 = np.ones(xstart0.shape)*-40 #  below the surface
	# zpar = 'fromMSL' 
	# pdb.set_trace()

	# for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
	do3d = 0
	# turbulence/diffusion flag. doturb=0 means no turb/diffusion,
	# doturb=1 means adding parameterized turbulence
	# doturb=2 means adding diffusion on a circle
	# doturb=3 means adding diffusion on an ellipse (anisodiffusion)
	doturb = 0

	return loc,nsteps,ndays,ff,tseas,ah,av,z0,zpar,do3d,doturb


def locations(test):
	'''
	Contains the locations and name information for the simulations.

	Inputs:
		test 	Index of the test case we want to run

	Outputs:
		lon0 	Drifter starting locations in x/zonal direction for test
		lat0 	Drifter starting locations in y/meridional direction for test
		name	Name of simulation to be used for netcdf file containing final tracks for test
	'''

	# Name of locations packages were found
	name = list((
			'Matagorda Island, TX',
			'Padre Island National Seashore, TX',
			'Matagorda Island, TX',
			'In Gulf of Mexico',
			'Mustang Island, TX',
			'Padre Island National Seashore, TX',
			'South Padre Island, TX',
			'Padre Island National Seashore, TX',
			'Padre Island National Seashore, TX',
			'Mustang Island, TX',		
			'Padre Island National Seashore, TX',
			'Padre Island National Seashore, TX',
			'Ship Island, MS',
			'Gilchrist Beach, TX',
			'Galveston Beach, Galveston, TX',
			'Galveston Beach, Galveston, TX',
			'Crystal Beach, TX',
			'Mud Lake, TX',
			'Galveston Beach, Galveston, TX',
			'High Island, TX',
			'Galveston Seawall, Galveston, TX',
			'Crystal Beach, TX',
			'McFadden Beach, TX',
			'Sabine Pass Jetties, TX',
			'Padre Island, TX',
			'McFaddin Beach, Sabine Pass, TX',	
			'McFaddin Beach, Sabine Pass, TX',		
			'Galveston Beach, Galveston, TX',
			'Texas Point, Sabine Pass, TX',
			'Texas Point, Sabine Pass, TX',
			'Sabine Pass, TX',
			'Sabine Pass, TX',
			'Sea Rim State Park, Sabine Pass, TX',
			'Sea Rim State Park, Sabine Pass, TX',
			'Pirates Beach, Galveston, TX',
			'Bolivar Peninsula, TX',
			'GOM, 3 miles south of Destin, FL',
			'Cameron, LA',
			'Cameron, LA',
			'Cameron, LA'))

	# Latitudes of where packages were found
	lat = np.array([28 + 10/60. + 8.03/3600.,
					27 + 25/60. + 26.96/3600.,
					28 + 10/60. + 8.03/3600.,
					np.nan,
					27 + 44/60. + 21.98/3600.,
					27 + 25/60. + 26.96/3600.,
					26 + 6/60. + 42.62/3600.,
					27 + 25/60. + 26.96/3600.,
					27 + 25/60. + 26.96/3600.,
					27 + 44/60. + 21.98/3600.,
					27 + 25/60. + 26.96/3600.,
					27 + 25/60. + 26.96/3600.,
					30 + 12/60. + 40.09/3600.,
					np.nan,
					29 + 16/60. + 16.79/3600.,
					29 + 16/60. + 16.79/3600.,
					29 + 27/60. + 25.82/3600.,
					29 + 32/60. + 33.32/3600.,
					29 + 16/60. + 16.79/3600.,
					29 + 34/60. + 0.83/3600.,
					29 + 18/60. + 8.85/3600.,
					29 + 27/60. + 25.82/3600.,
					29 + 39/60. + 24.04/3600.,
					29 + 41/60. + 40.00/3600.,
					27 + 15/60. + 49.24/3600.,
					29 + 39/60. + 24.04/3600.,
					29 + 39/60. + 24.04/3600.,
					29 + 16/60. + 16.79/3600.,
					29 + 41/60. + 3.88/3600.,
					29 + 41/60. + 3.88/3600.,
					29 + 41/60. + 3.88/3600.,
					29 + 41/60. + 3.88/3600.,
					29 + 41/60. + 12.84/3600.,
					29 + 41/60. + 12.84/3600.,
					29 + 12/60. + 15.35/3600.,
					29 + 28/60. + 41.84/3600.,
					np.nan,
					29 + 47/60. + 51.80/3600.,
					29 + 47/60. + 51.80/3600.,
					29 + 47/60. + 51.80/3600.])

	# Longitudes of where packages were found
	lon = np.array([96 + 44/60. + 19.30/3600.,
					97 + 17/60. + 57.06/3600.,
					96 + 44/60. + 19.30/3600.,
					np.nan,
					97 +  7/60. + 54.01/3600.,
					97 + 17/60. + 57.06/3600.,
					97 + 10/60. + 5.25/3600.,
					97 + 17/60. + 57.06/3600.,
					97 + 17/60. + 57.06/3600.,
					97 +  7/60. + 54.01/3600.,
					97 + 17/60. + 57.06/3600.,
					97 + 17/60. + 57.06/3600.,
					88 + 57/60. + 50.85/3600.,
					np.nan,
					94 + 49/60. + 25.87/3600.,
					94 + 49/60. + 25.87/3600.,
					94 + 38/60. + 22.71/3600.,
					95 +  0/60. + 56.82/3600.,
					94 + 49/60. + 25.87/3600.,
					94 + 23/60. + 36.68/3600.,
					94 + 46/60. + 27.70/3600.,
					94 + 38/60. + 22.71/3600.,
					94 +  6/60. + 20.05/3600.,
					93 + 51/60. + 9.07/3600.,
					97 + 21/60. + 49.13/3600.,
					94 +  6/60. + 20.05/3600.,
					94 +  6/60. + 20.05/3600.,
					94 + 49/60. + 25.87/3600.,
					93 + 51/60. + 2.66/3600.,
					93 + 51/60. + 2.66/3600.,
					93 + 51/60. + 2.66/3600.,
					93 + 51/60. + 2.66/3600.,
					94 +  2/60. + 38.40/3600.,
					94 +  2/60. + 38.40/3600.,
					94 + 56/60. + 10.69/3600.,
					94 + 34/60. + 47.69/3600.,
					np.nan,
					93 + 19/60. + 30.16/3600.,
					93 + 19/60. + 30.16/3600.,
					93 + 19/60. + 30.16/3600.])

	# Select out the lon/lat for test
	lon0,lat0 = np.meshgrid(np.linspace(lon[test]-0.1, lon[test]+0.1,30), \
							np.linspace(lat[test]+0.1, lat[test]+0.1,30))

	# Eliminate points that are outside domain or in masked areas
	lon0,lat0 = tools.check_points(lon0,lat0,grid)

	return lon0, lat0, name[test]


def start_times(test):
	'''
	Contains the times for finding the locations.

	Inputs:
		test 	Index of the test case we want to run

	Outputs:
		date 	Start date for test in datetime object
	'''

	dates = np.array([datetime(2013, 1, 11, 0),
					datetime(2013, 1, 16, 0),
					datetime(2013, 1, 18, 0),
					datetime(2013, 1, 25, 0),
					datetime(2013, 1, 29, 0),
					datetime(2013, 2, 10, 0),
					datetime(2013, 4, 10, 0),
					datetime(2013, 4, 13, 0),
					datetime(2013, 4, 23, 0),
					datetime(2013, 4, 30, 0),
					datetime(2013, 5, 7, 0),
					datetime(2013, 5, 7, 0),
					datetime(2013, 5, 12, 0),
					datetime(2013, 5, 19, 0),
					datetime(2013, 5, 24, 0),
					datetime(2013, 5, 26, 0),
					datetime(2013, 5, 26, 0),
					datetime(2013, 5, 28, 0),
					datetime(2013, 5, 29, 0),
					datetime(2013, 5, 29, 0),
					datetime(2013, 5, 30, 0),
					datetime(2013, 6, 2, 0),
					datetime(2013, 6, 2, 0),
					datetime(2013, 6, 3, 0),
					datetime(2013, 6, 3, 0),
					datetime(2013, 6, 7, 0),
					datetime(2013, 6, 7, 0),
					datetime(2013, 6, 7, 0),
					datetime(2013, 6, 9, 0),
					datetime(2013, 6, 9, 0),
					datetime(2013, 6, 10, 0),
					datetime(2013, 6, 10, 0),
					datetime(2013, 6, 10, 0),
					datetime(2013, 6, 10, 0),
					datetime(2013, 6, 11, 0),
					datetime(2013, 6, 11, 0),
					datetime(2013, 6, 12, 0),
					datetime(2013, 6, 13, 0),
					datetime(2013, 6, 13, 0),
					datetime(2013, 6, 13, 0)])

	return dates[test]