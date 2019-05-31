# Python file for running the Mauthner cell simulation
#
# Function run_model_somatic_stims runs by default the same simulation
# as mcell.hoc. Function run_model_dendritic stims runs
# another simulation, where the end of one of the dendrites
# is stimulated, and response is measured along the dendrites.
#
# HH formalism according to Buhry et al. 2013: "Global
# parameter estimation of an Hodgkin-Huxley formalism
# using membrane voltage recordings: Application to
# neuro-mimetic analog integrated circuits", Neurocomputing
# 81 (2012) 75-85
#
# Parameters obtained by hand-fitting and dimension-by-
# dimension local optimization, see an example in runfit.py
#
# Tuomo Maki-Marttunen, 2013-2017 (CC-BY 4.0)
#



import numpy as np
from neuron import h

def run_model_somatic_stims(params = [], stims = [[5.0, 10.0, 10.], [5.0, 10.0, 30.], [5.0, 10.0, 50.], [5.0, 10.0, 70.], [5.0, 10.0, 90.], [5.0, 10.0, 190.], [5.0, 10.0, 170.], [5.0, 10.0, 150.], [5.0, 10.0, 140.], [5.0, 10.0, 130.], [5.0, 10.0, 110.], [5.0, 10.0, 100.], [5.0, 10.0, 20.], [5.0, 10.0, -20.], [5.0, 10.0, -50.], [5.0, 10.0, -70.]], stimsR = [52.0, 72.0, 200.], recordDend = False, activesecs = ["dend[20]", "soma", "dend[21]", "dend[25]", "dend[27]", "dend[29]"], activecoeffs = [0.035, 0.035, 0.035, 0.035, 0.035, 0.0]): 

  if len(params)==0:
    params = [0.008699989010953,  4.454066447429189, 10.606380438297881, -83.389488181026280,  0.000955040259787, -52.802714891048296, -58.684517509388293,  8.103232645890719,  8.897735502857380,  0.016351981798972,  1.355604123083193, -63.065631347739128,  5.403637679068333,  0.114211915503005]

  gl=params[0]
  g1=params[1]
  g2=params[2]
  el=params[3]
  e1=55
  e2=-90
  Cm=2.5
  Ra=120
  glA=params[4]
  VoffaNa=params[5]   
  VoffaK=params[6]    
  VsloaNa=params[7]   
  VsloaK=params[8]    
  tauaNa=params[9]    
  tauaK=params[10]   
  VoffiNa=params[11]
  VsloiNa=params[12]
  tauiNa=params[13]

  t_sim = 15

  dt = 0.01

  alllengths = [28.263, 65.172, 137.513, 128.311, 39.547, 48.292, 49.275, 27.859, 36.188, 47.769, 49.665, 60.565, 69.584, 34.191, 125.280, 15.586, 72.951, 51.087, 59.876, 46.099, 32.208, 45.457, 75.944]
  allcumlengths = [25+x for x in [0.0, 28.26286, 93.43471, 230.94766, 230.94766, 270.49515, 270.49515, 93.43471, 28.26286]]+[25+x for x in [0.00000, 47.76941, 97.43476, 157.99953, 227.58324, 261.77471, 261.77471, 277.36042, 277.36042, 157.99953, 97.43476, 143.53365, 143.53365]]+[0.00]

  dendxs = [allcumlengths[i]+alllengths[i]/2 for i in range(len(alllengths))]
  for i in range(0,7):
    dendxs.append(100)
  axondiams = [54,54,54]
  h.load_file("neurmorph_activedend.hoc")
  h("access soma")
  h("soma.cm="+str(Cm))
  h("soma.Ra="+str(Ra))
  h("soma.e_pas="+str(el))
  h("soma.g_pas="+str(gl))
  for i in range(0,30):
    h("dend["+str(i)+"].cm="+str(Cm))
    h("dend["+str(i)+"].Ra="+str(Ra))
    h("dend["+str(i)+"].e_pas="+str(el))
    h("dend["+str(i)+"].g_pas="+str(gl))
  h("axonhillock.cm="+str(Cm))
  h("axonhillock.Ra="+str(Ra))
  h("axonhillock.e_pas="+str(el))
  h("axonhillock.g_pas="+str(glA))
  h("axonhillock.g_I1="+str(g1))
  h("axonhillock.g_I2="+str(g2))
  h("axonhillock.E_I1="+str(e1))
  h("axonhillock.E_I2="+str(e2))
  h("axonhillock.Voffa_I1="+str(VoffaNa))
  h("axonhillock.Voffa_I2="+str(VoffaK))
  h("axonhillock.Vsloa_I1="+str(VsloaNa))
  h("axonhillock.Vsloa_I2="+str(VsloaK))
  h("axonhillock.taua_I1="+str(tauaNa))
  h("axonhillock.taua_I2="+str(tauaK))
  h("axonhillock.Voffi_I1="+str(VoffiNa))
  h("axonhillock.Vsloi_I1="+str(VsloiNa))
  h("axonhillock.taui_I1="+str(tauiNa))
  for i in range(0,3):
    h("axon["+str(i)+"].cm="+str(Cm))
    h("axon["+str(i)+"].Ra="+str(Ra))
    h("axon["+str(i)+"].diam="+str(axondiams[i]))  
    h("axon["+str(i)+"].e_pas="+str(el))
    h("axon["+str(i)+"].g_pas="+str(glA))
  for i in range(0,len(activesecs)):
    h(activesecs[i]+".g_I1="+str(activecoeffs[i]*g1))
    h(activesecs[i]+".g_I2="+str(activecoeffs[i]*g2))
    h(activesecs[i]+".E_I1="+str(e1))
    h(activesecs[i]+".E_I2="+str(e2))
    h(activesecs[i]+".Voffa_I1="+str(VoffaNa))
    h(activesecs[i]+".Voffa_I2="+str(VoffaK))
    h(activesecs[i]+".Vsloa_I1="+str(VsloaNa))
    h(activesecs[i]+".Vsloa_I2="+str(VsloaK))
    h(activesecs[i]+".taua_I1="+str(tauaNa))
    h(activesecs[i]+".taua_I2="+str(tauaK))
    h(activesecs[i]+".Voffi_I1="+str(VoffiNa))
    h(activesecs[i]+".Vsloi_I1="+str(VsloiNa))
    h(activesecs[i]+".taui_I1="+str(tauiNa))
  h("forall nseg=20")

  h("objref stims[1]")
  h("soma stims[0] = new IClamp(0.5)")

  h("""
	v_init = """ + str(el) + """
	tstop = """ + str(t_sim) + """
        dt = """ + str(dt) + """

        cvode_active(1)
        cvode.atol(0.00005)
	objref time, vrec

	time = new Vector()
        time.record(&t)
        vrec = new Vector()
        vrec.record(&dend[1].v(0.5))
  """)
  dists_rec = []
  reclocs_branch = []
  if recordDend:
    h("distance()")
    reclocs_seg =    [-1, 17, 0,  4  ,10, 16,   16,  16,   16, 16,   16,  16,   16,21,  21, 21,  21,25,  25, 25, 25, 25, 27,  27, 27,  27,29,  29, 29,  29, 0]
    reclocs_x =      [0.5,0.5,0.5,0.5,0.5,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1, 0.25,0.5,0.75,1, 0.25,0.5,0.5,0.75,1, 0.25,0.5,0.75,1, 0.25,0.5,0.75,1, 0.5]
    reclocs_branch = [-1, 0,  0,  0,  0,  0,    0,   0,    0,  0,    0,   0,    0, 1,   1,  1,   1, 1,   1,  1,  1,  1,  1,   1,  1,   1, 1,   1,  1,   1, -2]
    h("""
objref recs
recs = new List()
""")
    for irec in range(0,len(reclocs_seg)):
      h("{recs.append(new Vector())}")
      if reclocs_seg[irec]==-1:
        h("{recs.o["+str(irec)+"].record(&soma.v("+str(reclocs_x[irec])+"))}")
        dists_rec.append(h.distance(reclocs_x[irec],sec=h.soma))
      elif reclocs_seg[irec]==-2:
        h("{recs.o["+str(irec)+"].record(&axonhillock.v("+str(reclocs_x[irec])+"))}")
        dists_rec.append(h.distance(reclocs_x[irec],sec=h.axonhillock))
      else:
        h("{recs.o["+str(irec)+"].record(&dend["+str(reclocs_seg[irec])+"].v("+str(reclocs_x[irec])+"))}")
        dists_rec.append(h.distance(reclocs_x[irec],sec=h.dend[reclocs_seg[irec]]))


  Vrecs = np.empty((np.shape(stims)[0]+1,), dtype=np.object)
  times = np.empty((np.shape(stims)[0]+1,), dtype=np.object)
  VrecsDend = []
  for istim in range(0,np.shape(stims)[0]):
    h("stims[0].del = "+str(stims[istim][0]))
    h("stims[0].dur = "+str(stims[istim][1]-stims[istim][0]))
    h("stims[0].amp = "+str(stims[istim][2]))
    h.init()
    h.run()

    Vrecs[istim] = np.array(h.vrec)
    times[istim] = np.array(h.time)
    if recordDend:
      VrecsDend.append(np.array(h.recs))

  h("stims[0].amp = 0")
  h("objref stimsR[50]")
  for i in range(0,50):
    h("soma stimsR["+str(i)+"] = new IClamp(0.5)")
  for i in range(0,50):
    h("stimsR["+str(i)+"].del = "+str(stimsR[0] + (stimsR[1]-stimsR[0])/50.0*i))
    h("stimsR["+str(i)+"].dur = "+str((stimsR[1]-stimsR[0])/50.0))
    h("stimsR["+str(i)+"].amp = "+str(stimsR[2]/50.0*i))
  t_sim = 72
  h("tstop = " + str(t_sim))
  h.init()
  h.run()

  Vrecs[np.shape(stims)[0]] = np.array(h.vrec)
  times[np.shape(stims)[0]] = np.array(h.time)

  return [times, Vrecs, VrecsDend, dists_rec, reclocs_branch]



def run_model_dendritic_stims(params = [], stimloc = 'lateral', stim_onset = 5, stim_dur=0.1, stim_amps = [10,20,40,60,80,100,120,150,200], idendstims = [16, 29], xdendstims = [1.0, 1.0], t_sim=15, activesecs = ["dend[20]", "soma", "dend[21]", "dend[25]", "dend[27]", "dend[29]"], activecoeffs = [0.035, 0.035, 0.035, 0.035, 0.035, 0.0]): 

  if len(params)==0:
    params = [0.008700000039526154, 20.999999990461987, 15.28930140049916, -83.40000793442388, 0.0003000000045918104, -56.70000000361548, -67.49999990345302, 8.10000002060277, 9.570002271582146,
              0.017999999846640063, 1.399997381068154, -64.00000048114761, 6.060000000757244, 0.20999225221952442]

  gl=params[0]
  g1=params[1]
  g2=params[2]
  el=params[3]
  e1=55
  e2=-90
  Cm=2.5
  Ra=120
  glA=params[4]
  VoffaNa=params[5]   
  VoffaK=params[6]    
  VsloaNa=params[7]   
  VsloaK=params[8]    
  tauaNa=params[9]    
  tauaK=params[10]   
  VoffiNa=params[11]
  VsloiNa=params[12]
  tauiNa=params[13]

  axondiams = [54,54,54]

  dt = 0.01

  h.load_file("neurmorph_activedend.hoc")
  h("access soma")
  h("distance()")
  h("soma.cm="+str(Cm))
  h("soma.Ra="+str(Ra))
  h("soma.e_pas="+str(el))
  h("soma.g_pas="+str(gl))
  for i in range(0,30):
    h("dend["+str(i)+"].cm="+str(Cm))
    h("dend["+str(i)+"].Ra="+str(Ra))
    h("dend["+str(i)+"].e_pas="+str(el))
    h("dend["+str(i)+"].g_pas="+str(gl))
  h("axonhillock.cm="+str(Cm))
  h("axonhillock.Ra="+str(Ra))
  h("axonhillock.e_pas="+str(el))
  h("axonhillock.g_pas="+str(glA))
  h("axonhillock.g_I1="+str(g1))
  h("axonhillock.g_I2="+str(g2))
  h("axonhillock.E_I1="+str(e1))
  h("axonhillock.E_I2="+str(e2))
  h("axonhillock.Voffa_I1="+str(VoffaNa))
  h("axonhillock.Voffa_I2="+str(VoffaK))
  h("axonhillock.Vsloa_I1="+str(VsloaNa))
  h("axonhillock.Vsloa_I2="+str(VsloaK))
  h("axonhillock.taua_I1="+str(tauaNa))
  h("axonhillock.taua_I2="+str(tauaK))
  h("axonhillock.Voffi_I1="+str(VoffiNa))
  h("axonhillock.Vsloi_I1="+str(VsloiNa))
  h("axonhillock.taui_I1="+str(tauiNa))
  for i in range(0,3):
    h("axon["+str(i)+"].cm="+str(Cm))
    h("axon["+str(i)+"].Ra="+str(Ra))
    h("axon["+str(i)+"].diam="+str(axondiams[i]))  
    h("axon["+str(i)+"].e_pas="+str(el))
    h("axon["+str(i)+"].g_pas="+str(glA))
  for i in range(0,len(activesecs)):
    h(activesecs[i]+".g_I1="+str(activecoeffs[i]*g1))
    h(activesecs[i]+".g_I2="+str(activecoeffs[i]*g2))
    h(activesecs[i]+".E_I1="+str(e1))
    h(activesecs[i]+".E_I2="+str(e2))
    h(activesecs[i]+".Voffa_I1="+str(VoffaNa))
    h(activesecs[i]+".Voffa_I2="+str(VoffaK))
    h(activesecs[i]+".Vsloa_I1="+str(VsloaNa))
    h(activesecs[i]+".Vsloa_I2="+str(VsloaK))
    h(activesecs[i]+".taua_I1="+str(tauaNa))
    h(activesecs[i]+".taua_I2="+str(tauaK))
    h(activesecs[i]+".Voffi_I1="+str(VoffiNa))
    h(activesecs[i]+".Vsloi_I1="+str(VsloiNa))
    h(activesecs[i]+".taui_I1="+str(tauiNa))

  h("forall nseg=20")

  h("objref stims[2]")
  h("dend["+str(idendstims[0])+"] stims[0] = new IClamp("+str(xdendstims[0])+")")
  h("dend["+str(idendstims[1])+"] stims[1] = new IClamp("+str(xdendstims[1])+")")

  reclocs_seg =    [-1, 17, 0,  4  ,10, 16,   16,  16,   16, 16,   16,  16,   16,21,  21, 21,  21,25,  25, 25, 25, 25, 27,  27, 27,  27,29,  29, 29,  29, 0]
  reclocs_x =      [0.5,0.5,0.5,0.5,0.5,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1, 0.25,0.5,0.75,1, 0.25,0.5,0.5,0.75,1, 0.25,0.5,0.75,1, 0.25,0.5,0.75,1, 0.5]
  reclocs_branch = [-1, 0,  0,  0,  0,  0,    0,   0,    0,  0,    0,   0,    0, 1,   1,  1,   1, 1,   1,  1,  1,  1,  1,   1,  1,   1, 1,   1,  1,   1, -2]

  h("""
	v_init = """ + str(el) + """
	tstop = """ + str(t_sim) + """
        dt = """ + str(dt) + """

        cvode_active(1)
        cvode.atol(0.00005)
	objref time, recs

	time = new Vector()
        time.record(&t)
        recs = new List()
""")

  dists_rec = []
  for irec in range(0,len(reclocs_seg)):
    h("{recs.append(new Vector())}")
    if reclocs_seg[irec]==-1:
      h("{recs.o["+str(irec)+"].record(&soma.v("+str(reclocs_x[irec])+"))}")
      dists_rec.append(h.distance(reclocs_x[irec],sec=h.soma))
    elif reclocs_seg[irec]==-2:
      h("{recs.o["+str(irec)+"].record(&axonhillock.v("+str(reclocs_x[irec])+"))}")
      dists_rec.append(h.distance(reclocs_x[irec],sec=h.axonhillock))
    else:
      h("{recs.o["+str(irec)+"].record(&dend["+str(reclocs_seg[irec])+"].v("+str(reclocs_x[irec])+"))}")
      dists_rec.append(h.distance(reclocs_x[irec],sec=h.dend[reclocs_seg[irec]]))
  dists_stim = []
  for istim in range(0,2):
    dists_stim.append(h.distance(xdendstims[istim],sec=h.dend[idendstims[istim]]))

  Vrecs = np.empty((len(stim_amps),), dtype=np.object)
  times = np.empty((len(stim_amps),), dtype=np.object)
  if stimloc=="lateral":
    istimloc = 0
  elif stimloc=="ventral":
    istimloc = 1
  else:
    print "Unknown stimulus location!"

  for istim in range(0,len(stim_amps)):
    h("stims["+str(istimloc)+"].del = "+str(stim_onset))
    h("stims["+str(istimloc)+"].dur = "+str(stim_dur))
    h("stims["+str(istimloc)+"].amp = "+str(stim_amps[istim]))
    h.init()
    h.run()

    Vrecs[istim] = np.array(h.recs)
    times[istim] = np.array(h.time)

  return [times, Vrecs, dists_rec, reclocs_branch, dists_stim]


