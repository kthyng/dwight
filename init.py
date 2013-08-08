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

# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

def parameters(loc, grid, date, pack):
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
    
    # Initialize parameters
    nsteps = 5
    ndays = 7
    ff = -1

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 20.
    av = 0

    # Number of drifters
    N = 1000

    # Read in location initializations
    lon0, lat0, packname = locations(pack, grid)
    lon0 = np.ones(N, np.dtype(int))*lon0
    lat0 = np.ones(N, np.dtype(int))*lat0

    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')

    # Initialize seed locations 
    ia = np.ceil(xstart0).astype(int) #[253]#,525]
    ja = np.ceil(ystart0).astype(int) #[57]#,40]

    # Choose method for vertical placement of drifters
    z0 = 's'
    zpar = 29

    do3d = 0
    doturb = 1

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    dostream = 1

    # convert date to number
    datenum = netCDF.date2num(date, units)

    # Number of model outputs to use
    tout = np.int((ndays*(24*3600))/tseas)

    # Figure out what files will be used for this tracking - to get tinds for
    # the following calculation
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)

    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0],grid,nc,z0,zpar)
    nc.close()

    # Initial total volume transport as a scalar quantity to be conserved, I think
    T0 = (abs(uf[ia, ja, 0]) + abs(vf[ia, ja, 0]))/N

    # Initialize the arrays to save the transports on the grid in the loop.
    # These arrays aggregate volume transport when a drifter enters or exits a grid cell
    # These should start at zero since we don't know which way things will travel yet
    U = np.ma.zeros(grid['xu'].shape,order='F')
    V = np.ma.zeros(grid['xv'].shape,order='F')

    # Add information to name
    name = str(pack) + '-' + date.isoformat()[0:13] + '-' + packname

    return nsteps, ndays, ff, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, dostream, T0, U, V, name


def seed(lon, lat, dlon=.5, dlat=.5, N=30):
    '''
    Chose array of starting locations based on the location
    the package was found (from locations()). A Gaussian 
    distribution is used to distribute points around the 
    find location.
    Inputs:
        lon, lat    Location where package was found
        dlon, dlat  Distance in degrees in which to seed drifters.
                    Default is 0.5 degrees.
        N           Number of drifters in x and y. Default is 30.

    Returns:
        lon0, lat0  Points in lon/lat at which to seed drifters
    '''

    # Find 2D distribution of points around the package location
    # the center is indicated using (lon, lat)
    # The variance is given by [[.25,0],[0,.25]] indicates
    #  a standard deviation away from the center of .5 degrees
    #  or .5**2=.25
    # There are N points in both x and y
    dist = np.random.multivariate_normal((lon, lat), \
                                    [[dlon**2,0],[0,dlat**2]], \
                                    [N,N])
    return dist[:,:,0], dist[:,:,1]

def locations(pack, grid):
    '''
    Contains the locations and name information for the simulations.

    Inputs:
        pack    Index of the pack case we want to run
        grid    grid dictionary as read in by tracpy.inout()

    Outputs:
        lon0    Drifter starting locations in x/zonal direction for pack
        lat0    Drifter starting locations in y/meridional direction for pack
        name    Name of simulation to be used for netcdf file containing final tracks for pack
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
            # 'GOM, 3 miles south of Destin, FL',
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
                    # 30 + 20/60. + 0/3600.,
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
                    # 86 + 30/60. + 0/3600.,
                    93 + 40/60. + 50/3600.,
                    93 + 41/60. + 31/3600.,
                    93 + 21/60. + 38/3600.])
    lon = -lon # add in negative sign!

    # Select location for this package
    lon0 = lon[pack]; lat0 = lat[pack];

    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')

    # Move all in the negative y direction, to avoid masked areas
    ystart0 = ystart0 - 2

    # Recover lon,lat
    lon0, lat0, _ = tracpy.tools.interpolate2d(xstart0, ystart0, \
                                             grid, 'm_ij2ll')

    # # masked is 1 if that starting location is masked out
    # masked = grid['mask'][ia.astype(int), ja.astype(int)]

    # # Try starting drifters at a single location instead of spread around.
    # # But, the found location might be masked out in the domain.
    # # If it is, then move in the y direction toward the ocean
    # for i in xrange(len(lat)):
    #     if masked[i] == 1.:
    #         # if masked, move 2 in y direction
    #         ja[i] = ja[i] - 2.

    # # Recover corresponding lat/lons
    # lon[masked.astype(int)], lat[masked.astype(int)], _ \
    #     = tracpy.tools.interpolate2d(ia[masked.astype(int)] - .5, \
    #                                  ja[masked.astype(int)] - .5, \
    #                                  grid, 'm_ij2ll')


    # # Select out the lon/lat for pack
    # dlon = 0.5; dlat = 0.5 # delta degree distances for starting particles

    # # Time Gaussian to set number of drifters used in (x,y)
    # H = np.arange(48) # hours in 2 days
    # mu = 24 # 1 day into the 2 days of simulation starts is the mean
    # sigma = 16. # Standard deviation
    # # pdb.set_trace()
    # N = 30*np.exp(-(H-mu)**2/(2*sigma**2))
    # # Choose N value for hour
    # Nh = np.floor(N[H==hour])
    # # N = 1/(sigma*sqrt(2*pi))*exp(-(H-mu)**2/(2*sigma**2))

    # lon0, lat0 = seed(lon[pack], lat[pack], dlon=dlon, dlat=dlat, N=Nh)
    # # lon0,lat0 = np.meshgrid(np.linspace(lon[pack]-dlon, lon[pack]+dlon,30), \
    # #                         np.linspace(lat[pack]-dlat, lat[pack]+dlat,30))

    # # pdb.set_trace()

    # # Eliminate points that are outside domain or in masked areas
    # lon0,lat0 = tracpy.tools.check_points(lon0,lat0,grid)

    return lon0, lat0, name[pack]

def start_times(pack):
    '''
    Contains the times for finding the locations.

    Inputs:
        pack    Index of the pack case we want to run

    Outputs:
        date    Start date for pack in datetime object
    '''

    dates = np.array([datetime(2013, 1, 11, 0, 1),
                    datetime(2013, 1, 16, 0, 1),
                    datetime(2013, 1, 18, 0, 1),
                    datetime(2013, 1, 25, 0, 1),
                    datetime(2013, 1, 29, 0, 1),
                    datetime(2013, 2, 10, 0, 1),
                    datetime(2013, 4, 10, 0, 1),
                    datetime(2013, 4, 13, 0, 1),
                    datetime(2013, 4, 23, 0, 1),
                    datetime(2013, 4, 30, 0, 1),
                    datetime(2013, 5, 7, 0, 1),
                    datetime(2013, 5, 7, 0, 1),
                    datetime(2013, 5, 12, 0, 1),
                    datetime(2013, 5, 19, 0, 1),
                    datetime(2013, 5, 24, 0, 1),
                    datetime(2013, 5, 26, 0, 1),
                    datetime(2013, 5, 26, 0, 1),
                    datetime(2013, 5, 28, 0, 1),
                    datetime(2013, 5, 29, 0, 1),
                    datetime(2013, 5, 29, 0, 1),
                    datetime(2013, 5, 30, 0, 1),
                    datetime(2013, 6, 2, 0, 1),
                    datetime(2013, 6, 2, 0, 1),
                    datetime(2013, 6, 3, 0, 1),
                    datetime(2013, 6, 3, 0, 1),
                    datetime(2013, 6, 7, 0, 1),
                    datetime(2013, 6, 7, 0, 1),
                    datetime(2013, 6, 7, 0, 1),
                    datetime(2013, 6, 9, 0, 1),
                    datetime(2013, 6, 9, 0, 1),
                    datetime(2013, 6, 10, 0, 1),
                    datetime(2013, 6, 10, 0, 1),
                    datetime(2013, 6, 10, 0, 1),
                    datetime(2013, 6, 10, 0, 1),
                    datetime(2013, 6, 11, 0, 1),
                    datetime(2013, 6, 11, 0, 1),
                    datetime(2013, 6, 12, 0, 1),
                    datetime(2013, 6, 13, 0, 1),
                    datetime(2013, 6, 13, 0, 1),
                    datetime(2013, 6, 13, 0, 1)])

    return dates[pack]