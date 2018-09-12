import numpy as np
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
from matplotlib import mlab
from matplotlib import rc
import argparse

rc('text', usetex=True)
np.random.seed(20)

def defdis(sig,l,nx,ny):
  # sig and l are parameters in the cov kernel
  # nx and ny are the number of cells along x and y
  max_def = 100
  x = np.linspace(0,1,nx)
  y = np.linspace(0,1,ny)
  xx, yy = np.meshgrid(x,y)
  S = (sig**2)*np.exp(-(xx-yy)**2/(2*l**2))
  D = np.zeros((S.shape[0],S.shape[1]))
  c = 0 # number of defects
  xD, yD = np.zeros((max_def,1)),np.zeros((max_def,1))

  for j in range(S.shape[0]):
    for i in range(S.shape[1]):
      sam = np.random.normal(0,S[j,i],1)
      D[j,i] = np.exp(sam)  
       
  for j in range(S.shape[0]):
    for i in range(S.shape[1]):
      if (D[j,i] > 0.8*np.amax(D)):
        xD[c], yD[c] = xx[j,i], yy[j,i] 
        c += 1
#        if (c == max_def):
#          break
#    else:
#      continue
#    break

# Plot the random field
  fig = plt.figure()
  ax = fig.add_subplot(111)
  plt.contourf(x,y,D)
  plt.colorbar()
  plt.hold(True)
  for i in range(c):
    ax.plot(xD[i],yD[i],'r*',markersize=12) 

  plt.title('$\mathrm{l = %3.2f}$' %l)

# Construct a filename automatically and replace '.' in l with a p 
  filename, ui = 'loc_defect_l', str(l)
  for i in range(len(ui)):
    if (ui[i] == '.'):
      p = i
      break
  filename = filename + ui[:p] + 'p' + ui[p+1:] + '.pdf'
  fig.savefig(filename)

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('corr_len',help='correlation length')
  args = parser.parse_args()
  l = float(args.corr_len)
  defdis(1.0,l,100,100)
