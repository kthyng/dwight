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
import tracpy


def parameters():
    '''
    Parameters for running numerical simulation.
    loc     Path to directory of grid and output files
    nsteps  Number of steps to do between model outputs (iter in tracmass)
    ndays   number of days to track the particles from start date
    ff      ff=1 to go forward in time and ff=-1 for backward in time
    tseas   Time between outputs in seconds
    ah      Horizontal diffusion in m^2/s. 
            See project values of 350, 100, 0, 2000. For -turb,-diffusion
    av      Vertical diffusion in m^2/s.
    do3d    for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    doturb  turbulence/diffusion flag. 
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
    xp      x-locations in x,y coordinates for drifters
    yp      y-locations in x,y coordinates for drifters
    zp      z-locations (depths from mean sea level) for drifters
    t       time for drifter tracks
    '''

    # Location of TXLA forecast model output file and then grid. 
    loc = ['/pong/raid/kthyng/forecast/txla_oof_his_jan_jul_2013.nc', \
            'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc']

    # Initialize parameters
    nsteps = 5
    ndays = 5 #16
    ff = -1
    # date = datetime(2009,11, 20, 0)


    # Time between outputs
    # Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 5. #100.
    av = 1.e-5 # m^2/s, or try 5e-6

    # Number of model outputs to use
    tout = np.int((ndays*(24*3600))/tseas)

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

    return loc,nsteps,ndays,ff,tseas,ah,av,z0,zpar,do3d,doturb,tout


def locations(test,grid):
    '''
    Contains the locations and name information for the simulations.

    Inputs:
        test    Index of the test case we want to run

    Outputs:
        lon0    Drifter starting locations in x/zonal direction for test
        lat0    Drifter starting locations in y/meridional direction for test
        name    Name of simulation to be used for netcdf file containing final tracks for test
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
    lat = np.array([28 + 9/60. + 22/3600.,
                    27 + 3/60. + 35/3600.,
                    28 + 15/60. + 16/3600.,
                    28 + 3.7/60.,
                    27 + 44/60. + 21.98/3600.,
                    27 + 26/60. + 14/3600.,
                    26 + 10/60. + 35/3600.,
                    27 + 28/60. + 13/3600.,
                    26 + 42/60. + 5/3600.,
                    27 + 47/60. + 9/3600.,
                    26 + 54/60. + 57/3600.,
                    26 + 59/60. + 13/3600.,
                    30 + 12/60. + 54/3600.,
                    29 + 29/60. + 18/3600.,
                    29 + 14/60. + 59/3600.,
                    29 + 14/60. + 59/3600.,
                    29 + 25/60. + 44/3600.,
                    29 + 35/60. + 53/3600.,
                    29 + 15/60. + 1/3600.,
                    29 + 33/60. + 21/3600.,
                    29 + 17/60. + 19/3600.,
                    29 + 25/60. + 44/3600.,
                    29 + 39/60. + 20/3600.,
                    29 + 42/60. + 13/3600.,
                    27 + 34/60. + 55/3600.,
                    29 + 39/60. + 20/3600.,
                    29 + 39/60. + 20/3600.,
                    29 + 10/60. + 48/3600.,
                    29 + 40/60. + 45/3600.,
                    29 + 40/60. + 45/3600.,
                    29 + 40.97/60.,
                    29 + 40.88/60.,
                    29 + 40.94/60.,
                    29 + 40.9/60.,
                    29 + 12/60. + 9/3600.,
                    29 + 29/60. + 31/3600.,
                    30 + 20/60. + 0/3600.,
                    29 + 44/60. + 48/3600.,
                    29 + 44/60. + 42/3600.,
                    29 + 45/60. + 49/3600.])

    # Longitudes of where packages were found
    lon = np.array([96 + 43/60. + 52/3600.,
                    97 + 22/60. + 44/3600.,
                    96 + 35/60. + 39/3600.,
                    96 + 29.7/60.,
                    97 +  7/60. + 54.01/3600.,
                    97 + 17/60. + 23/3600.,
                    97 + 10/60. + 28/3600.,
                    97 + 13/60. + 36/3600.,
                    97 + 19/60. + 22/3600.,
                    97 +  5/60. + 13/3600.,
                    97 + 22/60. + 18/3600.,
                    97 + 22/60. + 39/3600.,
                    88 + 56/60. + 42/3600.,
                    94 + 31/60. + 18/3600.,
                    94 + 51/60. + 17/3600.,
                    94 + 51/60. + 17/3600.,
                    94 + 40/60. + 28/3600.,
                    94 + 15/60. + 41/3600.,
                    94 + 51/60. + 15/3600.,
                    94 + 22/60. + 19/3600.,
                    94 + 47/60. + 19/3600.,
                    94 + 40/60. + 28/3600.,
                    94 +  6/60. + 38/3600.,
                    93 + 49/60. + 4/3600.,
                    97 + 13/60. + 4/3600.,
                    94 +  6/60. + 38/3600.,
                    94 +  6/60. + 38/3600.,
                    94 + 58/60. + 26/3600.,
                    93 + 51/60. + 12/3600.,
                    93 + 51/60. + 12/3600.,
                    93 + 52.64/60.,
                    93 + 53.96/60.,
                    93 + 57.48/60.,
                    93 + 59.29/60.,
                    94 + 56/60. + 8/3600.,
                    94 + 32/60. + 19/3600.,
                    86 + 30/60. + 0/3600.,
                    93 + 40/60. + 50/3600.,
                    93 + 41/60. + 31/3600.,
                    93 + 21/60. + 38/3600.])
    lon = -lon # add in negative sign!

    # Select out the lon/lat for test
    dlon = 0.5; dlat = 0.5 # delta degree distances for starting particles
    lon0,lat0 = np.meshgrid(np.linspace(lon[test]-dlon, lon[test]+dlon,30), \
                            np.linspace(lat[test]-dlat, lat[test]+dlat,30))

    # pdb.set_trace()

    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tracpy.tools.check_points(lon0,lat0,grid)

    return lon0, lat0, name[test]


def start_times(test):
    '''
    Contains the times for finding the locations.

    Inputs:
        test    Index of the test case we want to run

    Outputs:
        date    Start date for test in datetime object
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