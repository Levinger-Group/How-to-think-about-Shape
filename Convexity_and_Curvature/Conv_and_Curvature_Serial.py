#####################################################################################################################
#
# COMPUTING CONVEXITY AND CURVATURE
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

# force print to write the entire array out
np.set_printoptions(threshold=sys.maxsize)

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
while path.exists(Outputs_Path+TrajName+r'_Shape_Analysis{}.log'.format(Lognum)):
	Lognum += 1
Log = open(Outputs_Path+TrajName+r'_Shape_Analysis{}.log'.format(Lognum),'w')

# Log file header
Log.write('\n##########################################################################\n# Performing Shape Analysis of the {} Trajectory\n##########################################################################\n\n'.format(TrajName))
Log.flush()

#####################################
# Defining convexity and curvature
# calculating function
#####################################

def Conv_and_Curve(u,Group):
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
# Defining main iterator loop
#####################################

# Main loop
def Creator_Loop(u,CoreSelect,nframes,step):
	# setting a start time for timing information in log file
	starttime = time.time()
	# creating empty dictionaries for output
	conv_out = {}
	mcurve_out = {}
	gcurve_out = {}
	
	# counter for number of times through for loop, just for logging purposes
	p = 0
	selection = u.select_atoms(CoreSelect)
	
	# iterating over trajectory to max of nframes and only doing every [step]th frame
	for ts in u.trajectory[:nframes:step]:
		# getting current time of trajectory
		utime = u.trajectory.time
		# this is meant to ensure we have a complete log even if we hit the time limit. Gives just under 24 hours
		if (starttime - initialization) > 85000:
			  break
		# computing shape metrics at time point
		Conv,mcurve,gcurve = Conv_and_Curve(u,selection)
		# adding shape metrics to output dictionary
		conv_out.update({utime:Conv})
		mcurve_out.update({utime:mcurve})
		gcurve_out.update({utime:gcurve})
		
		# periodically updating log file to check progress
		if p%5000 == 0:
			Log.write('Computed time {} ps\n'.format(utime))
			Log.flush()
		p += 1
			
	return conv_out,mcurve_out,gcurve_out
	
#####################################
# Functions to write to file
#####################################

# there exists some limits in how many characters and columns you can throw in a csv with 
# csv.writer. So, we are bypassing that with some (simple) custom csv writing for the curvatures.
def CurvWriter(file,dict):
	with open(file,'w') as fh:
		for key in dict.keys():
			fh.write(str(key)+','+','.join([str(i) for i in dict[key]])+'\n')
			
def ConvWriter(file,dict):
	with open(file,'w') as csv_file:
		writer = csv.writer(csv_file,lineterminator='\n')
		for key, value in dict.items():
			writer.writerow([key,value])
	
#####################################
# Computing over Trajectory
#####################################

# initiating MDA Universe 
u = mda.Universe(Path+Top,Path+Traj)
	
Log.write('Finished reading in trajectory\nStarting main computation\n')
Log.flush()		# DO WHAT I SAY AND WRITE

# running actual computation of shape metrics
master_conv,master_mcur,master_gcur = Creator_Loop(u,CoreSelection,nframes,step)
	
Log.write('Finished computing shape analyses\nWriting results to csv files\n')
Log.flush()		# DO WHAT I SAY AND WRITE

#####################################
# Writing Results to File
#####################################

ConvWriter(Outputs_Path+TrajName+'_Convexity.csv',master_conv)
CurvWriter(Outputs_Path+TrajName+'_mCurvature.csv',master_mcur)
CurvWriter(Outputs_Path+TrajName+'_gCurvature.csv',master_gcur)

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

