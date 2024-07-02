#####################################################################################################################
#
# COMPUTING CONVEXITY AND CURVATURE
#
# OPTIONS: I created a log file to track per frame timings and then discovered you can get a
# performance report that does the same thing only significantly better. I left them both in. To get rid of
# the log, just comment out the creation of hte log file, starting with "Lognum" and ending with
# "Log = open(..." and then every line with "Log." after that.
# To get rid of the performance report, comment out the "with performance_report..." line and then
# remove the indent from all lines under that. 
#
# Here we compute the convexity estimator using the 3D extension of the area based method
# I got the idea from the reference listed (although I believe it originated from a previous work).
# In addition, we are computing the mean and Gaussian curvatures of the surface without abridgement.
# See the corresponding jupyter notebook for potential ways to reduce the output file size.
#
# Zunic, J.; Rosin, P.L.; A Convexity Measurement for Polygons. In Proceedings of the British Machine
# Vision Conference, Cardiff, UK, Sept 2-5, 2002; Rosin, P.L.; Marshall, A.D., Eds.; British Machine
# Vision Association, 173-182. DOI: 10.5244/C.16.15
#
#####################################################################################################################

#################################################################################
# Here's the code to read in the files created by this code.
# It's not exactly pretty, but it works fine and I refuse to comment it.
# I wrote it and I'm still fairly certain it's black magic. 
# It works. Just use it and don't ask questions. 
#
#def convexity_reader(IO):
#    output = np.empty((0,2))
#    with open(IO,newline='') as csvfile:
#        file = csv.reader(csvfile,delimiter=',')
#        for row in file:
#            output = np.append(output,np.array([[row[0],row[1]]]),axis=0)
#    return output.astype(np.float64)
#
#def curvature_reader(IO):
#    output = {}
#    with open(IO) as file:
#        lines = file.readlines()
#        for line in lines:
#            line = line.strip('\n')
#            if ',' in line:    # times and start of values are only lines with comma
#                if not line.split(',')[0] == '0.0':
#                    output.update({time:array})
#                    time = line.split(',')[0]
#                    array = []
#                    parsed = line.split('[')[1].split(' ')
#                    while('' in parsed):
#                        parsed.remove('')
#                    [array.append(float(parsed[i])) for i in range(len(parsed))]
#                else:
#                    time = line.split(',')[0]
#                    array = []
#                    parsed = line.split('[')[1].split(' ')
#                    while('' in parsed):
#                        parsed.remove('')
#                    [array.append(float(parsed[i])) for i in range(len(parsed))]
#            else:
#                line = line.strip(']')
#                parsed = line.split(' ')
#                while('' in parsed):
#                        parsed.remove('')
#                [array.append(float(parsed[i+1])) for i in range(len(parsed)-1)]
#    return output
#
#alternate curvature reader which only finds a specific time point's curvature
#def curvature_single_time(IO,time):
#    flipper = 0
#    with open(IO) as file:
#        while True:
#            line = file.readline().strip('\n')
#            if ',' in line:    # times and start of values are only lines with comma
#                if line.split(',')[0] == time:
#                    flipper = 1
#                    array = []
#                    parsed = line.split('[')[1].split(' ')
#                    while('' in parsed):
#                        parsed.remove('')
#                    [array.append(float(parsed[i])) for i in range(len(parsed))]
#                elif flipper == 1:
#                    output = np.asarray(array)
#                    return output
#            elif flipper == 1:
#                line = line.strip(']')
#                parsed = line.split(' ')
#                while('' in parsed):
#                    parsed.remove('')
#                [array.append(float(parsed[i+1])) for i in range(len(parsed)-1)]
# 
#################################################################################

#####################################
# Preamble and set up
#####################################

# Loading packages
import numpy as np
import MDAnalysis as mda
from scipy import spatial
import pytim
import pyvista as pv
from pytim.datafiles import *
import time
import math
import csv
import sys
from os import path
import multiprocessing as mp
import dask
import dask.multiprocessing
from dask.distributed import Client, LocalCluster, performance_report

# force print to write the entire array out
np.set_printoptions(threshold=sys.maxsize)
# sets up dask to run multiprocessing rather than multithreading
dask.config.set(scheduler='processes')

#####################################
# CHANGE-ABLES
# 
# It was convenient to put everything
# you'd change for analyzing a different
# trajectory up here at the top.
#####################################

# Trajectory
Top = 'Your_Topology.gro'
Traj = 'Your_Trajectory.xtc'

# Settings

# how many frames to analyze, e.g. "analyze every [step]th frame"
step = 4
# how many total frames to analyze. It was convenient for parallel version to specify exact number
# and for that number to be a nice, round, divisible number. It's less important for Serial
nframes = 51000

# specifying a suffix for your log file
Log_name = '_fullSA_Dask'
# I'm also using a Dask profiler to log resource management much better
# for multiprocessing. This sets the name.
Profile_name = 'My_Profile_Name-Dask-profile.html'

# Core selection
CoreSelection = '(resname AOT and (type O or type S or name C1 or name H1 or name C2 or name H2 or name H3 or name C3 or name C12)) or resname SOL'

# radii dictionary for creating Willard-Chandler surface
radiidict = pytim_data.vdwradii(CHARMM27_TOP)

#####################################
# Set Up
#
# Setting up input and output paths
# Creating a log file with a unique number
# Log file is mostly for tracking time per frame
# and total time for the analysis.
# Comment out all Log lines if you don't need this
#
# Side note: if you run this on slurm and it fails for any reason,
# slurm may continue to attempt to re-run the code and you'll end up with
# hundreds of log files. 
#####################################

# Loading in files
Path = '/projects/cgale@colostate.edu/GMX/Production/'
Outputs_Path = '/projects/cgale@colostate.edu/GMX/Python_Outputs/'

# Recording start time
initialization = time.time()

# Creating a log file with a new name to avoid overwrite issues and closing issues
Lognum = 1
TrajName = Top.split('.')[0]
while path.exists(Outputs_Path+TrajName+Log_name+r'{}.log'.format(Lognum)):
	Lognum += 1
Log = open(Outputs_Path+TrajName+Log_name+r'{}.log'.format(Lognum),'w')


# Log file header
Log.write('\n##########################################################################\n# Performing Shape Analysis of the {} Trajectory\n##########################################################################\n\n'.format(TrajName))
Log.flush()		# DO WHAT I SAY AND WRITE

#####################################
# Defining convexity and curvature
# calculating function
#####################################

def Conv_and_Curve(frame,u,Group):
	# setting frame
	Group.universe.trajectory[frame]

	# getting trajectory time
	utime = u.trajectory.time

	# Getting Surface
	WC = pytim.WillardChandler(u,group=Group,alpha=3.0,mesh=1.1,fast=False,radii_dict=radiidict)
	
	# converting PyTim to PyVista surface
	verts = WC.triangulated_surface[0]
	faces = WC.triangulated_surface[1]
	threes = 3*np.ones((faces.shape[0],1),dtype=int)
	faces = np.concatenate((threes,faces),axis=1)
	Poly = pv.PolyData(verts,faces)
	
	# Getting actual surface's volume
	volume = Poly.volume
	
	# Getting convex hull's volume
	CH = spatial.ConvexHull(verts)
	CH_volume = CH.volume
	
	# computing convexity
	convexity = volume/CH_volume
	
	# Getting curvature data
	mean_curv = Poly.curvature(curv_type='mean')
	G_curv = Poly.curvature(curv_type='Gaussian')
	
	# Computing convexity estimator
	return np.array([utime,convexity]),np.array([mean_curv]),np.array([G_curv])
	
#####################################
# Defining a Batch Operator
#####################################

def Batch_Computer(range,u,Group):
	# defining output array for convexity
	out_conv = np.empty((len(range),2))
	# defining output list for curvature (variable length, can't use array)
	out_mcurv = []
	out_gcurv = []
	# setting a number to count iterations
	call = 0
	for i in range:
		# we're sharing convexities time with all metrics
		out_conv[call,:],mcurv,gcurv = Conv_and_Curve(i,u,Group)
		out_mcurv.append(mcurv)
		out_gcurv.append(gcurv)
		call += 1
	return out_conv,out_mcurv,out_gcurv
	
#####################################
# Functions to write to file
#####################################

# there exists some limits in how many characters and columns you can throw in a csv with 
# csv.writer. So, we are bypassing that with some (simple) custom csv writing for the curvatures.
def CurvWriter(file,val_list,Conv_array):
	with open(file,'w') as fh:
		for i in range(len(val_list)):
			fh.write(str(Conv_array[i,0])+','+','.join([str(i) for i in val_list[i]])+'\n')
			
def ConvWriter(file,array):
	with open(file,'w') as fh:
		for i in range(array.shape[0]):
			fh.write(str(array[i,0])+','+str(array[i,1])+'\n')
			
#####################################
# Setting up I/O
#####################################

# generating trajectories
Log.write('Reading in trajectories and selecting micelle core\n')
Log.flush()		# DO WHAT I SAY AND WRITE

# initiating MDA Universe
u = mda.Universe(Path+Top,Path+Traj)
Core = u.select_atoms(CoreSelection)
	
Log.write('Finished reading in trajectories\n')
Log.flush()		# DO WHAT I SAY AND WRITE

#####################################
# Computing over Trajectory
#####################################

if __name__ == '__main__':
	# getting local architecture for Dask-y boy
	client = Client()
	# we're using the Dask performance report to track timing and resource management
	# This runs the code inside of the performance report so that we actually track those resources
	with performance_report(filename=Outputs_Path+Profile_name):
		Log.write('There are {} cores available\n'.format(len(client.ncores().values())))
		Log.flush()
		
		# batching
		partitions = len(client.ncores().values())
		batch_size = nframes*step//partitions

		# creating start/stop ranges
		start = np.empty(partitions)
		stop = np.empty(partitions)
		for i in range(partitions):
			start[i] = i*batch_size
			stop[i] = (i+1)*batch_size
		# adjusting last batch to make sure we get all frames (accounts for remainders)
		stop[-1] = nframes

		# Creating list of dask.delayed jobs: lazy functions that can be evaluated later
		job_list = []
		for i in range(partitions):
			job_list.append(dask.delayed(Batch_Computer)(range(int(start[i]),int(stop[i]),int(step)),
														u,
														Core))

		'''
		Notes on running and output
		
		dask.compute(*job_list) will evaluate the delayed objects directly and return a Collections object
		client.compute(job_list) will create a list of futures which have not been evaluated yet
		The syntactic difference with the asterisk IS important
		
		To evaluate the futures and return a list, we need to (naively) call the "Future.result()" method
		But this would create a list of results that need compiling. We call the "Client.gather()" method,
		which is faster than calling "result()" individually. Note that the result() method does block the
		code until computations are complete. However, it returns the results in the same list format as the
		Futures. So you still need to compile the results back together.
		'''
		# Still not actually evaluating, just creating a set of futures
		futures = client.compute(job_list)
		
		Log.write('Created futures\nBlocking here for evaluation of futures\n')
		Log.flush()		# DO WHAT I SAY AND WRITE
		
		# actual computation takes place here. Code is blocked here until complete
		result = client.gather(futures)
		
		Log.write('Finished evaluating futures\n')
		Log.flush()		# DO WHAT I SAY AND WRITE
		
		# repacking results into singular list
		# Conv_and_Curve creates 3 outputs, results will then create a list with as many entries
		# as there are cores being utilized. In this case, there are 6, so it will be a list with the
		# approximate shape (6,3) (although those entries are themselves lists and arrays of varying shape)
		# this code just combines those entries regardless of number of cores used and splits the function
		# output into individual outputs. With the batching I've done above, the results should always be
		# in order. Dask appears to handle that much, at least. 
		master_conv = np.empty((0,2))
		master_mcurv = []
		master_gcurv = []
		for i in range(len(result)):
			master_conv = np.concatenate((master_conv,result[i][0]),axis=0)
			master_mcurv += result[i][1]
			master_gcurv += result[i][2]

		Log.write('Writing results to csv files\n')
		Log.flush()		# DO WHAT I SAY AND WRITE

		#####################################
		# Writing Results to File
		#####################################

		ConvWriter(Outputs_Path+TrajName+'_Convexity.csv',master_conv)
		CurvWriter(Outputs_Path+TrajName+'_mCurvature.csv',master_mcurv,master_conv)
		CurvWriter(Outputs_Path+TrajName+'_gCurvature.csv',master_gcurv,master_conv)

#####################################
# Finishing Log file and wrapping up
#####################################

# Log wrap up,
# getting the final time
enditilization = time.time()
Total_time = (enditilization-initialization)/(60)

Log.write('Finished writing to file\nComputation complete\nbeginning Robot Revolution\n>>>END LOG<<<\n\n\n... ...\n... ...\n... ...\n...\nHave a nice day!\n\n#####################################\n# Computation took {:.5f} minutes\n#####################################'.format(Total_time))
Log.flush()		# DO WHAT I SAY AND WRITE

Log.close()

