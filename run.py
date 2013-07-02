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

# Read in parameter initialization
loc, nsteps, ndays, ff, tseas, ah, av, z0, zpar, do3d, doturb = init.parameters()

# Read in grid
grid = tracpy.inout.readgrid(loc,nc)

# Loop through start locations and times for running different simulations
for test in xrange(ntests):

	# Read in location initializations
	lon0, lat0, name = init.locations(test)

	# Read in time initializations
	date = init.start_times(test)

	# TODO: Try to put each simulation on a different core of the current machine, except 1 or 2
	lonp, latp, zp, t, grid = run(loc, nsteps, ndays, ff, date, \
									tseas, ah, av, lon0, lat0, \
									z0, zpar, do3d, doturb, name)

# Compile tex document with figures in it