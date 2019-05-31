# An example of the dimension-by-dimension parameter fitting.
# Tuomo Maki-Marttunen, 2013-2017 (CC-BY 4.0)

import scipy.io
from pylab import *
import minimizedimbydim
import mcell
import mytools


stims = [[5.0, 10.0, 10.], [5.0, 10.0, 30.], [5.0, 10.0, 50.], [5.0, 10.0, 70.], [5.0, 10.0, 90.], [5.0, 10.0, 190.], [5.0, 10.0, 170.], [5.0, 10.0, 150.], [5.0, 10.0, 140.], [5.0, 10.0, 130.], [5.0, 10.0, 110.], [5.0, 10.0, 100.], [5.0, 10.0, 20.], [5.0, 10.0, -20.], [5.0, 10.0, -50.], [5.0, 10.0, -70.]]
stimsR = [52.0, 72.0, 200.]
experimental_data = scipy.io.loadmat('experimental_data.mat')
experimental_spt = [mytools.spike_times(experimental_data['ts'].T[0], x, -62, inf) for x in experimental_data['medianVs'].T]+[mytools.spike_times(experimental_data['tsR'].T[0], experimental_data['medianVsR'].T[0], -62, inf)]

weight_ramp = 3
coeff_mV = 1./200
coeff_ms = 1./20

def objective_function(params_this,savefig=""):
  myparams = [0.0]*28
  myparams[0] = params_this[0] #0.008700000039526
  myparams[1] = params_this[1] #20.999999990461987
  myparams[2] = params_this[2] #15.289301400499159
  myparams[3] = params_this[3] #-83.400007934423883
  myparams[4] = params_this[4] #0.000300000004592
  myparams[5] = params_this[5] #-56.700000003615479
  myparams[6] = params_this[6] #-67.499999903453016
  myparams[7] = params_this[7] #8.100000020602771
  myparams[8] = params_this[8] #9.570002271582146
  myparams[9] = params_this[9] #0.017999999846640
  myparams[10] = params_this[10] #1.399997381068154
  myparams[11] = params_this[11] #-64.000000481147609
  myparams[12] = params_this[12] #6.060000000757244
  myparams[13] = params_this[13] #0.209992252219524
  return calc_error(myparams,savefig)

def calc_error(params_this,savefig=""):
  global stims, stimsR, experimental_data

  data = mcell.run_model_somatic_stims(params_this,stims,stimsR)
  times = data[0]
  Vrecs = data[1]

  errs = []
  for irun in range(0,len(times)):
    if irun == len(times)-1:
      tref = experimental_data['tsR'].T[0]
      vref = experimental_data['medianVsR'].T[0]
    else:
      tref = experimental_data['ts'].T[0]
      vref = experimental_data['medianVs'].T[irun]
    sptref = experimental_spt[irun]

    tthis = times[irun]
    vthis = Vrecs[irun]
    sptthis = mytools.spike_times(tthis, vthis, -62, inf)
    vthis = mytools.interpolate_extrapolate_constant(tthis,vthis,tref)

    meantracediff_this = 1.0*mean([abs(x-y) for x,y in zip(vthis, vref)])
    sp_N_err_this = abs((not len(sptthis))-(not len(sptref)))
    sp_t_err_this = 0
    if len(sptref) > 0:
      for ispike in range(0,min(len(sptthis),len(sptref))):
        sp_t_err_this = sp_t_err_this + min([abs(sptthis[ispike] - x) for x in sptref])

    if irun == len(times)-1:
      errs.append( weight_ramp * (coeff_mV * meantracediff_this + coeff_ms * sp_t_err_this + sp_N_err_this) )
    else:
      errs.append( coeff_mV * meantracediff_this + coeff_ms * sp_t_err_this + sp_N_err_this )

  if len(savefig) > 0:
    close("all")
    f,axs = subplots(1,2)
    for i in [0,1,2,3,4,11,12,13,14,15]:
      axs[0].plot(times[i],Vrecs[i],'b-')
      axs[0].plot(experimental_data['ts'].T[0],experimental_data['medianVs'].T[i],'r--')
    axs[1].plot(times[9],Vrecs[9],'b-')
    axs[1].plot(experimental_data['ts'].T[0],experimental_data['medianVs'].T[9],'r--')

    for i in range(0,2):
      axs[i].set_xlim([4,12])
      axs[i].set_xticks([5,7,9,11])
      axs[i].set_xticklabels(['0','2','4','6'])

    f.savefig(savefig)
  return sum(errs)

thrs = [ [0.002, 0.02],
         [0, 30],
         [0, 20],
         [-80, -90],
         [0.000, 0.01], #gleak_A
         [-60, -30],    #Voffa_Na
         [-70, -30],    #Voffa_K
         [6.0, 10.0],
         [7.0, 12.0],
         [0.015, 0.025],
         [1.0, 2.0],
         [-75, -25],    #Voffi_Na
         [4.0, 10.0],
         [0.1, 0.4]]

init_params = array([0.0087, 21., 15.2893, -83.4, 0.0003, -56.700000003615479, -67.499999939430012, 8.1, 9.57, 0.018, 1.399997369505575, -64, 6.06, 0.209992280912807])
init_error = objective_function(init_params,"init.eps")

params = minimizedimbydim.minimizedimbydim(lambda x: objective_function(x),array(thrs).T, init_params)
fitted_error = objective_function(params[0][0],"fitted.eps")
scipy.io.savemat('fitted.mat',{'params': params, 'init_params': init_params})

