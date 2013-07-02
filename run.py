

# Read in parameter initialization

# Read in location initializations

# Read in time initializations

loc,nsteps,ndays,ff,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name = init.parameters()

# Loop through start locations and times for running different simulations
for test in xrange(tests):
	# TODO: Try to put each simulation on a different core of the current machine, except 1 or 2
	lonp,latp,zp,t,grid = run(loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name)

# Compile tex document with figures in it