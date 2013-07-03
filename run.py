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

	# # Convert date to number
	# date = netCDF.date2num(date,units)

	# # Figure out what files will be used for this tracking
	# nc,tinds = tracpy.inout.setupROMSfiles(loc,date,ff,tout)

	# TODO: Try to put each simulation on a different core of the current machine, except 1 or 2
	lonp, latp, zp, t, grid = tracpy.run.run(loc, nsteps, ndays, ff, date, \
									tseas, ah, av, lon0, lat0, \
									z0, zpar, do3d, doturb, name)

	pdb.set_trace()

	# Plot tracks
	tracpy.plotting.tracks(lonp,latp,name,grid=grid)

	# Plot final location (by time index) histogram
	tracpy.plotting.hist(lonp,latp,name,grid=grid,which='contour')
	tracpy.plotting.hist(lonp,latp,name,grid=grid,which='pcolor')	

# Compile tex document with figures in it