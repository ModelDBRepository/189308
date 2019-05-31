# A script that runs mcell.hoc and plots the results (dashed) against the experimental data (solid).
# The figure fit.eps produced by the script corresponds to Figure 7B in the paper.
# Tuomo Maki-Marttunen, 2013-2017 (CC-BY 4.0)

import scipy.io
from pylab import *
from neuron import h
import mytools

h("""load_file("mcell.hoc")""")

experimental_data=scipy.io.loadmat('experimental_data.mat')
close("all")

f,axs = subplots(1,2)

stimAmps = [10,30,50,70,90,190,170,150,140,130,110,100,20,-20,-50,-70]
inds = [15,14,13,0,12,1,2,3,4,11]
cols = mytools.colorsredtolila(10,0.75)
myhandles = []
legends = []
for ii in range(0,len(inds)):
      i = inds[ii]
      axs[0].plot(experimental_data['ts'].T[0],experimental_data['medianVs'].T[i],'k-',linewidth=2)
      myhandles.append(axs[0].plot(h.timeList[i],h.vrecList[i],'r--', dashes=(2,1),color=cols[ii],linewidth=2,label=str(stimAmps[i])+' nA')[0])
      legends.append(str(stimAmps[i])+' nA')
axs[0].legend(myhandles,legends,fontsize=6)

axs[1].plot(experimental_data['ts'].T[0],experimental_data['medianVs'].T[9],'k-',linewidth=2)
myhandle=axs[1].plot(h.timeList[9],h.vrecList[9],'r--', dashes=(2,1),color=cols[0],linewidth=2,label=str(stimAmps[9])+' nA')
axs[1].legend(myhandle,[str(stimAmps[9])+' nA'],fontsize=6)

axs[0].set_xlim([4,12])
axs[0].set_xticks([5,7,9,11])
axs[0].set_xticklabels(['0','2','4','6'])
axs[0].set_title('Square 5ms',fontsize=10)

axs[1].set_xlim([4,7.5])
axs[1].set_xticks([5,6,7])
axs[1].set_xticklabels(['0','1','2'])
axs[1].set_title('Square 5ms',fontsize=10)

for i in range(0,2):
  axs[i].set_xlabel('$t$ (ms)',fontsize=10)
  axs[i].set_ylim([-92,-45])
  for tick in axs[i].xaxis.get_major_ticks()+axs[i].yaxis.get_major_ticks():
    tick.label.set_fontsize(6)

axs[0].set_ylabel('$V_m$ (mV)',fontsize=10)
axs[1].set_ylabel('$V_m$ (mV)',fontsize=10)
f.savefig("fit.eps")

