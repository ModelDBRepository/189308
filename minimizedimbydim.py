# Dimension-by-dimension local optimization algorithm
# Tuomo Maki-Marttunen, 2013-2017 (CC-BY 4.0)

from pylab import *
import time

def minimizedimbydim(fcn=lambda x: x[0]**2+x[1]**2,thrs=array([[-1.,-1.],[1.,1.]]),initx=[],powinterval=1./3,Ninterval=40,Niter=10,powintervalstd=[],nhoursmax=inf):

  if len(initx) == 0:
    initx = thrs[0] + (thrs[1]-thrs[0])*rand(1)

  if len(powintervalstd) == 0:
    powintervalstd = powinterval/5

  if len(initx.shape) == 1:
    initx = [initx]
  #if type(initx[0]) is not list:
  #  initx = [initx]
  myx = initx[:]
    
  Npart = len(myx)
  errs = []
  for j in range(0, Npart):
    errs.append(fcn(myx[j]))

  errsall = zeros([Niter,Npart])

  timenow = time.time()
  dimord = range(0,len(initx[0]))
  shuffle(dimord)
  #print str(powintervalstd)

  for iter in range(0,Niter):
    if time.time()-timenow > nhoursmax*3600:
      break
    dimord2 = range(0,len(initx[0])-1)
    shuffle(dimord2)
    dimord = [dimord[0]] + [dimord[1+i] for i in dimord2]

    for iidim in range(0,len(dimord)):
      idim = dimord[iidim]
      thispowinterval = powinterval + powintervalstd*randn(1)
      k=0
      while thispowinterval <= 0 or thispowinterval >= 1:
        thispowinterval = powinterval + powintervalstd*randn(1)
        k = k+1
        if k > 100000:
          print('try smaller powintervalstd!')
          return [nan,nan]
        
      parvals = []
      for j in range(0,Npart):
        parvals.append((myx[j][idim] - (myx[j][idim]-thrs[0][idim])*thispowinterval**array(range(1,Ninterval+1))).tolist() +
                       (myx[j][idim] + (thrs[1][idim]-myx[j][idim])*thispowinterval**array(range(1,Ninterval+1))).tolist())

      #print "        parvals[0]="+str(parvals[j][0])+", [end-1]="+str(parvals[j][Ninterval-1])+", val="+str(myx[j][idim])+", [2*end-1]="+str(parvals[j][2*Ninterval-1])+", [end]="+str(parvals[j][Ninterval])
      for j in range(0,Npart):
        these_errs = zeros([2*Ninterval,1])
        for k in range(0,len(parvals[0])):
          if parvals[j][k] == myx[j][idim]:
            these_errs[k] = errs[j]
          else:
            thisx = myx[j]
            thisx[idim] = parvals[j][k]
            these_errs[k] = fcn(thisx)

        if any(these_errs < errs[j]):
          ind = argmin(these_errs)
          myx[j][idim] = parvals[j][ind]
          errs[j] = min(these_errs)
    errsall[iter,:] = errs;
    print("minimizedimbydim: iter = "+str(iter)+", err = "+str(min(errs))) #+", x = "+str(myx))
  return [myx[:],errsall[:]]
