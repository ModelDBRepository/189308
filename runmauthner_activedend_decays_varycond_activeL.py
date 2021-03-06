# A script that runs the simulations with the dendritic stimulations (active dendrites) using
# various values for active conductances, fits the maximal membrane potential along the dendrites
# to an exponential decay and analyzes the results. Both dendrites are active and have a constant
# value for the active conductances along the dendrites.
# The script plots two figures:
# decays_robustness_activeL.eps:
#   Plots the obtained exponential decay constants for each ventral and lateral conductance value.
#   These data were not plotted in the paper but may help to understand the results.
# decay_orders_robustness_activeL.eps
#   Plots the order of decay constants between the two cases (orthodromic vs antidromic).
#   Corresponds to Figure 9B.
# Tuomo Maki-Marttunen, 2013-2017 (CC-BY 4.0)

import mcell_activedend_varycoeffs_activeL as mcell
from pylab import *
from neuron import h
import numpy.matlib
import pickle
import time
from os.path import exists

params = [ 0.008692502978958,  9.033350623275275,  4.380690156620405, -83.399741825852217,  0.000204081169084, -55.832300314690642, -67.436819173866425,  8.100021281440229,  9.559735981936562,  0.020391367006446,  1.412486203424452, -63.999689527874885,  6.046642114148621,  0.170597469041809 ] #(THE BEST ONE FOR ACTIVE VENTRAL DENDRITE)

stimbranches = ['lateral','ventral']
close("all")
f,axarr = subplots(1,2)
cols = ['#800000','#000080'][::-1]
linecols = ['#FF0000','#0000FF'][::-1]

# The variables 'cond_Ls' and 'cond_Vs' define the parameter space where the effect of the values of active conductances along the lateral and ventral dendrites is analyzed.
# The variable 'coeffs' will determine the active conductances along the dendrites. The first compartment ("dend[20]", i.e. previously passive AIS between soma and active AIS) is set
# the conductances such that Na+ channel conductance is coeffs[0]*(Na+ conductance at active AIS) and K+ channel conductance is coeffs[0]*(K+ conductance at active AIS). In a similar
# manner, the following compartments along the main branch of ventral dendrite until the end of the ventral dendrite ("soma", "dend[21]", "dend[25]", "dend[27]", "dend[29]") receive
# values that are scaled by coeffs[1-5] from the values at active AIS. Likewise, the compartments along the main branch of lateral dendrite ("dend[19]", "dend[17]", "dend[0]",
# "dend[4]", "dend[8]","dend[10]", "dend[16]") receive values that are scaled by coeffs[6-12] from the values at active AIS. Here, the values of 'coeffs' are constant within ventral
# dendrite (except the last compartment which is zero) and another constant within lateral dendrite (except the last compartment which is zero).
cond_Ls = [0.0+0.00125*i for i in range(0,41)] 
cond_Vs = [0.0+0.0025*i for i in range(0,51)]
amps = [0.5,1,1.5,2,2.5,3,4,5]

if exists('activedend_decays_varyconds_activeL.sav'): # If the result file already exists, use it.
 unpicklefile = open('activedend_decays_varyconds_activeL.sav', 'r')
 unpickledlist = pickle.load(unpicklefile)
 unpicklefile.close()
 cond_Ls = unpickledlist[0]
 cond_Vs = unpickledlist[1]
 cMs_all = unpickledlist[2]
 maxVs_all = unpickledlist[3]
 cMsA_all = unpickledlist[4]
 maxVsA_all = unpickledlist[5]
else: # Simulate the model for all combinations of the cond_Ls and cond_Vs and save the data in the result file
 cMs_all = []
 maxVs_all = []
 for icond_L in range(0,len(cond_Ls)):
  cMs_thiscond_L = []
  maxVs_thiscond_L = []
  print "ORTHO: icond_L = "+str(icond_L)+"/"+str(len(cond_Ls))
  for icond_V in range(0,len(cond_Vs)):
   cMs_thiscond_V = []
   maxVs_thiscond_V = []
   for ibranch in range(0,len(stimbranches)):
     branch = stimbranches[ibranch]
     coeffs = [cond_Vs[icond_V] for i in range(0,5)]+[0.0]+[cond_Ls[icond_L] for i in range(0,6)]+[0.0]
     data = mcell.run_model_dendritic_stims(params, branch, 5, 50.0, amps, [16, 29], [1.0, 1.0], 15, ["dend[20]", "soma", "dend[21]", "dend[25]", "dend[27]", "dend[29]", "dend[19]", "dend[17]", "dend[0]", "dend[4]", "dend[8]","dend[10]", "dend[16]"], coeffs)

     times = data[0]
     Vrecs = data[1]
     dists = data[2]
     branches = data[3]

     cMs = []
     print "coeffs = "+str(coeffs)
     for iamp in range(0,len(amps)):
       XL = array([dists[iloc] for iloc in range(0,len(dists)) if dists[iloc] < 400 and branches[iloc]==ibranch or branches[iloc]==-1])
       XM = matrix(c_[np.matlib.repmat(XL,1,1).T, kron([1],ones([len(XL),])).T])
       YL = [array([max(Vrecs[iamp][iloc]) for iloc in range(0,len(dists)) if dists[iloc] < 400 and branches[iloc]==ibranch or branches[iloc]==-1])]
       YM = matrix(concatenate(YL)).T
       thiserr = inf
  
       Vappr = params[3]

       cM = inv(XM.T*XM)*XM.T*(log(YM-Vappr))
       thiserr = sum(sum(abs(YM-Vappr-exp(XM*cM))))
       cMs.append(cM[:])
     cMs_thiscond_V.append(cMs[:])
     maxVs_thiscond_V.append([max(Vrecs[i][30]) for i in range(0,len(amps))])
   cMs_thiscond_L.append(cMs_thiscond_V[:])
   maxVs_thiscond_L.append(maxVs_thiscond_V[:])
  cMs_all.append(cMs_thiscond_L[:])
  maxVs_all.append(maxVs_thiscond_L[:])

 cMsA_all = []
 maxVsA_all = []
 for icond_L in range(0,len(cond_Ls)):
  cMsA_thiscond_L = []
  maxVsA_thiscond_L = []
  print "ANTI: icond_L = "+str(icond_L)+"/"+str(len(cond_Ls))
  for icond_V in range(0,len(cond_Vs)):
   cMsA_thiscond_V = []
   maxVsA_thiscond_V = []

   stim = []
   for i in range(0,len(amps)):
     stim.append([5,55,amps[i]])
   coeffs = [cond_Vs[icond_V] for i in range(0,5)]+[0.0]+[cond_Ls[icond_L] for i in range(0,5)]+[0.0]
   data = mcell.run_model_somatic_stims(params, stim,[5.0, 10.0, 0], True, ["dend[20]", "soma", "dend[21]", "dend[25]", "dend[27]", "dend[29]", "dend[19]", "dend[17]", "dend[0]", "dend[4]", "dend[10]", "dend[16]"], coeffs)

   times = data[0]
   Vrecs = data[1]
   VrecsDend = data[2]
   dists = data[3]
   branches = data[4]

   for ibranch in range(0,2):
     cMs = []
     cMolds = []
     for iamp in range(0,len(amps)):
       XL = array([dists[iloc] for iloc in range(0,len(dists)) if dists[iloc] < 400 and branches[iloc]==ibranch or branches[iloc]==-1])
       XM = matrix(c_[np.matlib.repmat(XL,1,1).T, kron([1],ones([len(XL),])).T])
       YL = [array([max(VrecsDend[iamp][iloc]) for iloc in range(0,len(dists)) if dists[iloc] < 400 and branches[iloc]==ibranch or branches[iloc]==-1])]
       YM = matrix(concatenate(YL)).T
       thiserr = inf
  
       Vappr = params[3]

       cM = inv(XM.T*XM)*XM.T*(log(YM-Vappr))
       thiserr = sum(sum(abs(YM-Vappr-exp(XM*cM))))
       cM = cM
       cMs.append(cM[:])
     cMsA_thiscond_V.append(cMs[:])
   maxVsA_thiscond_L.append([max(VrecsDend[i][30]) for i in range(0,len(amps))])
   cMsA_thiscond_L.append(cMsA_thiscond_V[:])
  cMsA_all.append(cMsA_thiscond_L[:])
  maxVsA_all.append(maxVsA_thiscond_L[:])

 picklelist = [cond_Ls, cond_Vs, cMs_all, maxVs_all, cMsA_all, maxVsA_all]
 file = open('activedend_decays_varyconds_activeL.sav', 'w')
 pickle.dump(picklelist,file)
 file.close()

cMsL = [[[log(3.0)/cMs_all[i][j][0][k][0][0,0] for j in range(0,len(cMs_all[i]))] for i in range(0,len(cMs_all))] for k in range(0,len(amps))]
cMsV = [[[log(3.0)/cMs_all[i][j][1][k][0][0,0] for j in range(0,len(cMs_all[i]))] for i in range(0,len(cMs_all))] for k in range(0,len(amps))]
mcMsL = mean(array(cMsL),axis=0)
mcMsV = mean(array(cMsV),axis=0)

cMsAL = [[[-log(3.0)/cMsA_all[i][j][0][k][0][0,0] for j in range(0,len(cMsA_all[i]))] for i in range(0,len(cMsA_all))] for k in range(0,len(amps))]
cMsAV = [[[-log(3.0)/cMsA_all[i][j][1][k][0][0,0] for j in range(0,len(cMsA_all[i]))] for i in range(0,len(cMsA_all))] for k in range(0,len(amps))]
mcMsAL = mean(array(cMsAL),axis=0)
mcMsAV = mean(array(cMsAV),axis=0)

f,axarr = subplots(2,2)
a1 = axarr[0,0].imshow(mcMsL,extent=[min(cond_Vs),max(cond_Vs),min(cond_Ls),max(cond_Ls)],interpolation='nearest')
colorbar(a1,ax = axarr[0,0], orientation='horizontal')
axarr[0,0].set_title('Ortho, L')

a2 = axarr[0,1].imshow(mcMsV,extent=[min(cond_Vs),max(cond_Vs),min(cond_Ls),max(cond_Ls)],interpolation='nearest')
colorbar(a2,ax = axarr[0,1], orientation='horizontal')
axarr[0,1].set_title('Ortho, C')

a3 = axarr[1,0].imshow(mcMsAL,extent=[min(cond_Vs),max(cond_Vs),min(cond_Ls),max(cond_Ls)],interpolation='nearest')
colorbar(a3,ax = axarr[1,0], orientation='horizontal')
axarr[1,0].set_title('Anti, L')

a4 = axarr[1,1].imshow(mcMsAV,extent=[min(cond_Vs),max(cond_Vs),min(cond_Ls),max(cond_Ls)],interpolation='nearest')
colorbar(a4,ax = axarr[1,1], orientation='horizontal')
axarr[1,1].set_title('Anti, V')

f.savefig("decays_robustness_activeL.eps")


mask_chosen = zeros(array(cMsL[0]).shape)
mask_chosen[find(array(cond_Ls)==0.0)[0],find(array(cond_Vs)==0.035)[0]] = True
f,axarr = subplots(1,1)
data = (12*logical_and(greater(numpy.min(array(cMsL[1:]),axis=0),numpy.max(array(cMsV[1:]),axis=0)),greater(numpy.min(array(cMsAV[1:-2]),axis=0),numpy.max(array(cMsAL[1:-2]),axis=0))) + \
        7.5*logical_and(greater(numpy.min(array(cMsL[1:]),axis=0),numpy.max(array(cMsV[1:]),axis=0)),logical_not(greater(numpy.min(array(cMsAV[1:-2]),axis=0),numpy.max(array(cMsAL[1:-2]),axis=0)))) + \
        8*logical_and(logical_not(greater(numpy.min(array(cMsL[1:]),axis=0),numpy.max(array(cMsV[1:]),axis=0))),greater(numpy.min(array(cMsAV[1:-2]),axis=0),numpy.max(array(cMsAL[1:-2]),axis=0))) + \
        6*logical_and(logical_not(greater(numpy.min(array(cMsL[1:]),axis=0),numpy.max(array(cMsV[1:]),axis=0))),logical_not(greater(numpy.min(array(cMsAV[1:-2]),axis=0),numpy.max(array(cMsAL[1:-2]),axis=0)))) + \
        -0.5*logical_and(greater(numpy.min(array(cMsL[1:]),axis=0),numpy.max(array(cMsV[1:]),axis=0)),less(numpy.max(array(cMsAV[1:-2]),axis=0),numpy.min(array(cMsAL[1:-2]),axis=0))) + \
        0.5*logical_and(less(numpy.max(array(cMsL[1:]),axis=0),numpy.min(array(cMsV[1:]),axis=0)),greater(numpy.min(array(cMsAV[1:-2]),axis=0),numpy.max(array(cMsAL[1:-2]),axis=0))) \
        )*(1-mask_chosen)

data2 = zeros([2*x for x in data.shape])
for i in range(0,data.shape[0]):
  for j in range(0,data.shape[1]):
    if abs(data[i,j]-12) < 1e-7 or abs(data[i,j]-7.0) < 1e-7 or abs(data[i,j]-0) < 1e-7 or abs(data[i,j]-8.5) < 1e-7:
      data2[2*i:2*i+2,2*j:2*j+2] = data[i,j]
    elif abs(data[i,j]-7.5) < 1e-7:
      data2[2*i:2*i+2,2*j:2*j+2] = array([[7,12],[12,7]])
    elif abs(data[i,j]-8.0) < 1e-7:
      data2[2*i:2*i+2,2*j:2*j+2] = array([[12,8],[8,12]])
    elif abs(data[i,j]-6.5) < 1e-7 or abs(data[i,j]-6.0) < 1e-7:
      data2[2*i:2*i+2,2*j:2*j+2] = array([[7,8],[8,7]])
a3 = axarr.imshow(data2, extent=[min(cond_Vs),max(cond_Vs),max(cond_Ls),min(cond_Ls)],interpolation='nearest', cmap="hot")

axarr.set_ylim([0.0,0.05])
axarr.set_aspect('auto')
axarr.set_position([0.2,0.33,0.6,0.5])
axarr.set_xlabel('$c_{VD}$')
myyl = axarr.set_ylabel('$c_{LD}$')
myyl.set_rotation(0)
axnew = [0,0,0,0,0]
axnew[0] = f.add_axes([0.18, 0.19, 0.03, 0.03],axisbg='w')
axnew[1] = f.add_axes([0.18, 0.11, 0.03, 0.03],axisbg='w')
axnew[2] = f.add_axes([0.48, 0.19, 0.03, 0.03],axisbg='w')
axnew[3] = f.add_axes([0.48, 0.11, 0.03, 0.03],axisbg='w')
axnew[4] = f.add_axes([0.48, 0.03, 0.03, 0.03],axisbg='w')
a4 = axnew[0].imshow(array([[12,12],[12,12]]), interpolation='nearest', cmap="hot", vmin=0, vmax=12)
a4 = axnew[1].imshow(array([[7,7],[7,7]]), interpolation='nearest', cmap="hot", vmin=0, vmax=12)
a4 = axnew[2].imshow(array([[7,12],[12,7]]), interpolation='nearest', cmap="hot", vmin=0, vmax=12)
a4 = axnew[3].imshow(array([[12,8],[8,12]]), interpolation='nearest', cmap="hot", vmin=0, vmax=12)
a4 = axnew[4].imshow(array([[7,8],[8,7]]), interpolation='nearest', cmap="hot", vmin=0, vmax=12)
for ax in axnew:
  ax.set_xticks([])
  ax.set_yticks([])
axnew[0].text(4,2,'L(O) $>$ V(O), L(A) $<$ V(A)\nfor each amplitude',fontsize=10)
axnew[1].text(4,2,'L(O) $>$ V(O), L(A) $>$ L(A)\nfor each amplitude',fontsize=10)
axnew[2].text(4,2,'L(O) $>$ V(O) for each amplitude,\nL(A) $>$ V(A) for some amplitudes',fontsize=10)
axnew[3].text(4,2,'L(O) $<$ V(O) for some amplitudes,\nL(A) $<$ V(A) for each amplitude',fontsize=10)
axnew[4].text(4,2,'L(O) $<$ V(O) for some amplitudes,\nL(A) $>$ V(A) for some amplitudes',fontsize=10)
f.savefig("decay_orders_robustness_activeL.eps")

