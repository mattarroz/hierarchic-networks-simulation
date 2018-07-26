# coding: utf-8
import matplotlib
matplotlib.use('Agg')


import stochpy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#plt=matplotlib.pyplot

smod = stochpy.SSA()
smod.model_dir='./'
smod.Model('example5_deterministic_initial_mRNA_steady.psc')

smod.DoStochSim(end = 15.0,mode = 'time',trajectories = 100000)
with PdfPages('Transc-Transl-mRNA-steady-det-init-trajectories.pdf') as pdf:
	fig=smod.PlotAverageSpeciesTimeSeries()
	stochpy.plt.title("")
	stochpy.plt.ylabel("Cells")
	stochpy.plt.ylabel("")
	stochpy.plt.xlabel("Time [minutes]")
	stochpy.plt.legend(["mRNA","Protein"],numpoints=1,loc='center right')
	pdf.savefig(fig)

smod.GetRegularGrid()
# species[species_index][trajectory_index][time_step_index]
protarray=np.array(smod.data_stochsim_grid.species[1])
# we transpose to get the time step index as the first index
protarray=protarray.transpose()
# select last time step
#ts=len(protarray)-1

for ts in range(0,len(protarray)):
	simtime = smod.data_stochsim_grid.getTime()[ts]
	if max(protarray[ts] != 0):
		
		fig=plt.figure()
		plt.ioff()
		with PdfPages('Transc-Transl-mRNA-steady-det-init-'+str(simtime)+'.pdf') as pdf:
			plt.hist(protarray[ts],range=[0,max(protarray[ts])],bins=max(protarray[ts]))
			plt.title("")
			plt.ylabel("Frequency")
			plt.xlabel("Number of molecules")
#			plt.show()
#+			plt.pause(0.001)
			pdf.savefig(fig)
	
		hist = np.histogram(protarray[ts],range=[0,max(protarray[ts])],bins=max(protarray[ts]))
		np.savetxt("transc-trans-mRNA-steady-det-init-"+str(simtime)+".csv",hist[1],delimiter=",")
