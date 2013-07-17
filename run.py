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
import time
import tracpy
import init
from scipy import ndimage

ntests = 39 # number of packages found

# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

# Read in parameter initialization
loc, nsteps, ndays, ff, tseas, ah, av, z0, zpar, do3d, doturb, tout = init.parameters()

# Get the first nc object to initialize grid
# Read in time initializations
# date = datetime(2009,12,2,0)
date = init.start_times(0)
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
    lon0, lat0, testname = init.locations(test,grid)

    # Read in time initializations
    date = init.start_times(test)
    # date = datetime(2009,12,2,0)

    datevec = date - timedelta(days=2)
    for hour in range(0,48,4):
        dat = datevec + timedelta(hours=hour)

        # Add information to name
        name = str(test) + '-' + str(dat.year) + '-' + str(dat.month).zfill(2) \
             + '-' + str(dat.day).zfill(2) + '-' + str(dat.hour).zfill(2) + '-' + testname


        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc'):
            # TODO: Try to put each simulation on a different core of the current machine, except 1 or 2
            lonp, latp, zp, t, grid = tracpy.run.run(loc, nsteps, ndays, ff, dat, \
                                            tseas, ah, av, lon0, lat0, \
                                            z0, zpar, do3d, doturb, name)

        else: # if the files already exist, just read them in for plotting
            d = netCDF.Dataset('tracks/' + name + '.nc')
            lonp = d.variables['lonp'][:]
            latp = d.variables['latp'][:]

        # If the particle trajectories have not been plotted, plot them
        if not os.path.exists('figures/' + name + 'tracks.png'):
            tracpy.plotting.tracks(lonp, latp, name, grid=grid)
        if not os.path.exists('figures/' + name + 'histhexbin.png'):
            tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')

# pdb.set_trace()
# Make histogram of all final locations
d = netCDF.MFDataset('tracks/*',aggdim='ntrac')
name = 'overall'
lonp = d.variables['lonp'][:]; latp = d.variables['latp'][:]
tracpy.plotting.hist(lonp,latp, name, tind='vector', \
                    grid=grid, which='hexbin')
tracpy.plotting.tracks(lonp, latp, name, grid=grid)

# Compile tex document with figures in it
# !pdflatex dwight.tex