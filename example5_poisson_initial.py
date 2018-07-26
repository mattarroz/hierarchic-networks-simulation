# coding: utf-8
import stochpy
import numpy as np
import math
import pickle
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages

def poisson(mu,k):
	return (mu**k/math.factorial(k))*math.exp(-mu);

def simulate(mu,initmRNA):
	smod = stochpy.SSA()
	smod.Model('Transc-Transl.psc')

 	smod.ChangeInitialSpeciesCopyNumber("S1",initmRNA)
	traj = max(10,int(round(1e6*(poisson(mu,initmRNA)))))
	print "Simulating " + str(traj) + " trajectories"
	smod.DoStochSim(end = 5,mode = 'time',trajectories = traj)

	smod.GetRegularGrid()

	fp = open("poisssim"+str(mu)+"_"+str(initmRNA)+".pkl", "wb")
	pickle.dump(smod.data_stochsim_grid, fp)
	fp.close()
	del smod
	return;

mu = 1.2
histlist = []
maxmRNA = 10
initval = range(0,maxmRNA-1)
maxlen = 0

#Parallel(n_jobs=maxmRNA)(delayed(simulate)(mu,initmRNA) for initmRNA in initval)

for initmRNA in initval:
	simulate(mu,initmRNA)

simtime = np.zeros(shape=(maxmRNA,100))

for initmRNA in initval:
	print "Processing simulation, initmRNA="+str(initmRNA)
	fp = open ("poisssim"+str(mu)+"_"+str(initmRNA)+".pkl", "rb")
	data_stochsim_grid = pickle.load(fp)
	fp.close()
	# species[species_index][trajectory_index][time_step_index]
	protarray=np.array(data_stochsim_grid.species[1])
	# we transpose to get the time step index as the first index
	protarray=protarray.transpose()
	
	histlistlist = []
	np.resize(simtime,(initmRNA,len(protarray)))
	for ts in range(0,len(protarray)):
		simtime[initmRNA][ts] = data_stochsim_grid.getTime()[ts]
		if max(protarray[ts] != 0):
			hist = np.histogram(protarray[ts],range=[0,max(protarray[ts])],bins=max(protarray[ts]))
			histlistlist.append(hist[0])
			if maxlen < len(hist[0]):
				maxlen = len(hist[0])
	histlist.append(histlistlist)

histlist = np.array(histlist)
histlist = histlist.transpose()
histpoisson = np.zeros(shape=(len(protarray),maxlen))
for initRNA in initval:
	for ts in range(0,len(protarray)-1):
		for k in range(0,len(histlist[ts][initmRNA])-1):
			histpoisson[ts][k] += poisson(mu,initmRNA)*histlist[ts][initmRNA][k]
		np.savetxt("transc-trans-poiss-"+str(simtime[initmRNA][ts])+".csv",histpoisson[ts],delimiter=",")

