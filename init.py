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
import inout
import tools

def parameters(loc,nsteps,ndays,ff,tseas,ah,av,z0,zpar,do3d,doturb):
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

	# grid = netCDF.Dataset(loc+'grid.nc')
	# lonr = grid.variables['lon_rho'][:]
	# latr = grid.variables['lat_rho'][:]
	grid = inout.readgrid(loc)


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

	# simulation name, used for saving results into netcdf file
	name = 'temp_mapcoord_order1' #'5_5_D5_F'

	# # Save these settings to a file
	# np.savez('tracks/' + name + 'header.in',loc=loc,nsteps=nsteps,ndays=ndays, 
	# 		ff=ff,date=date,tseas=tseas,ah=ah,av=av,lon0=lon0,lat0=lat0,
	# 		z0=z0,zpar=zpar,do3d=do3d,doturb=doturb,name=name)

	return loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name


def locations():
	'''
	lon0 	Drifter starting locations in x/zonal direction.
	lat0 	Drifter starting locations in y/meridional direction.
	name	Name of simulation to be used for netcdf file containing final tracks
	'''

	## Input starting locations as real space lon,lat locations
	# lon0,lat0 = np.meshgrid(-95.498218005315309,23.142258627126882) # [0,0] (SE) corner
	# lon0,lat0 = np.meshgrid(-97.748582291691989,23.000027311710628) # [-1,0] (SW) corner
	# lon0,lat0 = np.meshgrid(-87.757124031927574,29.235771320764623) # [0,-1] (NE) corner
	# lon0,lat0 = np.meshgrid(-88.3634073986196,30.388542615201313) # [-1,-1] (NW) corner
	# lon0,lat0 = np.meshgrid(np.linspace(-94,-93,10),np.linspace(28,29,10)) # grid outside Galveston Bay
	# lon0,lat0 = np.meshgrid(np.linspace(-95,-91,100),np.linspace(28,29,50)) # rectangle outside Galveston

	# lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,1100),np.linspace(22.5,31,980)) # whole domain, 1 km
	# lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,220),np.linspace(22.5,31,196)) # whole domain, 5 km
	# FOR TEST1:
	lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,110),np.linspace(22.5,31,98)) # whole domain, 10 km
	# lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,21),np.linspace(22.5,31,20)) # whole domain, 50 km

	# Eliminate points that are outside domain or in masked areas
	lon0,lat0 = tools.check_points(lon0,lat0,grid)

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

	lat =  28°10'8.03"N
	27°25'26.96"N
	28°10'8.03"N
	np.nan
	 27°44'21.98"N
	27°25'26.96"N
	 26° 6'42.62"N
	27°25'26.96"N
	27°25'26.96"N
	 27°44'21.98"N
	27°25'26.96"N
	27°25'26.96"N
	 30°12'40.09"N
	 np.nan
	  29°16'16.79"N
	  29°16'16.79"N
	  29°27'25.82"N
	   29°32'33.32"N
	  29°16'16.79"N
	  29°34'0.83"N
	  29°18'8.85"N
	  29°27'25.82"N
	  29°39'24.04"N
	  29°41'40.00"N
	   27°15'49.24"N
	  29°39'24.04"N
	  29°39'24.04"N
	  29°16'16.79"N
	  29°41'3.88"N
	  29°41'3.88"N
	  29°41'3.88"N
	  29°41'3.88"N
	  29°41'12.84"N
	  29°41'12.84"N
	  29°12'15.35"N
	  29°28'41.84"N
	  np.nan
	  29°47'51.80"N
	  29°47'51.80"N
	  29°47'51.80"N

	lon =  96°44'19.30"W
	 97°17'57.06"W
	 96°44'19.30"W
	 np.nan
	  97° 7'54.01"W
	 97°17'57.06"W
	  97°10'5.25"W
	 97°17'57.06"W
	 97°17'57.06"W
	  97° 7'54.01"W
	 97°17'57.06"W
	 97°17'57.06"W
	  88°57'50.85"W
	  np.nan
	  94°49'25.87"W
	  94°49'25.87"W
	   94°38'22.71"W
	    95° 0'56.82"W
	  94°49'25.87"W
	  94°23'36.68"W
	   94°46'27.70"W
	   94°38'22.71"W
	   94° 6'20.05"W
	   93°51'9.07"W
	   97°21'49.13"W
	   94° 6'20.05"W
	   94° 6'20.05"W
	   94°49'25.87"W
	   93°51'2.66"W
	   93°51'2.66"W
	   93°51'2.66"W
	   93°51'2.66"W
	   94° 2'38.40"W
	   94° 2'38.40"W
	   94°56'10.69"W
	   94°34'47.69"W
	   np.nan
	   93°19'30.16"W
	   93°19'30.16"W
	   93°19'30.16"W

def dates():
	'''
	date 	Start date in datetime object
	'''

	# Start date
	date = datetime(2009,11, 25, 0)

