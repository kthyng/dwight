import numpy as np
import sys
import os
import op
import tracmass
import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap
import pdb
from matplotlib import delaunay
from matplotlib.pyplot import *
import glob
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import time
import tracpy
import init
from scipy import ndimage

ntests = 40 # number of packages found

# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

# Read in parameter initialization
loc, nsteps, ndays, ff, tseas, ah, av, z0, zpar, do3d, doturb, tout = init.parameters()

# Get the first nc object to initialize grid
# Read in time initializations
date = datetime(2009,12,2,0)
# date = init.start_times(0)
# Convert date to number
date = netCDF.date2num(date,units)
# Figure out what files will be used for this tracking
nc, _ = tracpy.inout.setupROMSfiles(loc,date,ff,tout)
# Read in grid
grid = tracpy.inout.readgrid(loc,nc)

lonpsave = np.ones(40*900)*np.nan
latpsave = np.ones(40*900)*np.nan

# Loop through start locations and times for running different simulations
for test in xrange(ntests):

	# Read in location initializations
	lon0, lat0, name = init.locations(test,grid)

	# Read in time initializations
	# date = init.start_times(test)

	date = datetime(2009,12,2,0)

	# pdb.set_trace()

	# Add information to name
	name = str(test) + '-' + str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '-' + name

	# If the particle trajectories have not been run, run them
	if not os.path.exists('tracks/' + name + '.nc'):
		# TODO: Try to put each simulation on a different core of the current machine, except 1 or 2
		lonp, latp, zp, t, grid = tracpy.run.run(loc, nsteps, ndays, ff, date, \
										tseas, ah, av, lon0, lat0, \
										z0, zpar, do3d, doturb, name)

	else: # if the files already exist, just read them in for plotting
		d = netCDF.Dataset('tracks/' + name + '.nc')
		lonp = d.variables['lonp'][:]
		latp = d.variables['latp'][:]

	# pdb.set_trace()
	ln = lonp.shape[1]
	lonpsave[ln*test:ln*test+ln] = lonp[-1,:]
	latpsave[ln*test:ln*test+ln] = latp[-1,:]

	# Plot tracks
	# pdb.set_trace()
	# tracpy.plotting.tracks(lonp,latp,name,grid=grid)

	# Plot final location (by time index) histogram
	# tracpy.plotting.hist(lonp,latp,name,grid=grid,which='contour')
	# xmin, ymin = grid['basemap'](lonp.min()-.1, latp.min()+.1)
	# xmax, ymax = grid['basemap'](lonp.max()+.1, latp.max()+.1)
	# tracpy.plotting.hist(lonp,latp,name,grid=grid, \
	# 						which='pcolor',bins=(80,80), \
	# 						xlims=[xmin, xmax], \
	# 						ylims=[ymin, ymax])	

	# pdb.set_trace()
	# ADD ABILITY TO JUST READ IN TRACKS IF ALREADY DONE

# pdb.set_trace()
# Make histogram of all final locations
tracpy.plotting.hist(lonpsave,latpsave,name,grid=grid,tind='vector', \
							which='pcolor',bins=(80,80))

# Compile tex document with figures in it