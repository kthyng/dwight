import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import sys
import os
import op
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

packs = 39 # number of packages found

# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

loc = ['/pong/raid/kthyng/forecast/txla_oof_his_jan_jul_2013.nc', \
        'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc']

# Read in grid
grid = tracpy.inout.readgrid(loc)

# Loop through packages found
for pack in xrange(packs):

    # Start simulations from the date the package was found
    date = init.start_times(pack)

    name = date.isoformat()[0:13]

    # Initialize counter for number of hours to increment through simulation by
    nh = 0

    # Start simulations for a span of 2 days before the package was found
    while nh < 48:

        # Read in parameter initialization
        nsteps, ndays, ff, tseas, ah, av, lon0, lat0, z0, zpar, \
            do3d, doturb, dostream, T0, U, V, name = init.parameters(loc, grid, date, pack)

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc'):

            # Run tracpy
            lonp, latp, zp, t, grid, T0, U, V \
                = tracpy.run.run(loc, nsteps, ndays, ff, date, tseas, ah, av, \
                                    lon0, lat0, z0, zpar, do3d, doturb, name, \
                                    grid=grid, dostream=dostream, T0=T0, U=U, V=V)

        # If basic figures don't exist, make them
        if not os.path.exists('figures/' + name + '*.png'):

            # Read in and plot tracks
            d = netCDF.Dataset('tracks/' + name + '.nc')
            lonp = d.variables['lonp'][:]
            latp = d.variables['latp'][:]
            tracpy.plotting.tracks(lonp, latp, name, grid=grid)
            tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')

        # Increment by 4 hours for next loop
        nh = nh + 4
        date = date + timedelta(hours=nh)

# Compile tex document with figures in it
# !pdflatex dwight.tex