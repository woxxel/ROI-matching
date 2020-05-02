import numpy as np
import cv2,sys
from scipy.io import loadmat
import matplotlib.pyplot as plt
from scipy import signal


sys.path.append('/home/wollex/Data/Science/PhD/Programs/PC_analysis')
from utils import pathcat, get_shift_and_flow, normalize_sparse_array

def get_session_perspective(basePath,mouse,s_in):
  
  pathMouse = pathcat([basePath,mouse])
  
  dims = (512,512)
  Cn = np.zeros((2,512,512))
  for i,s in enumerate(s_in):
    pathSession = pathcat([pathMouse,'Session%02d'%s])
    pathResults = pathcat([pathSession,'results_redetect.mat'])
    ld = loadmat(pathResults,variable_names=['A'])
    A = ld['A']
    #A=normalize_sparse_array(A)
    Cn[i,...] = A.sum(1).reshape(512,512).astype('float32')
  
  shift, flow, (x_grid,y_grid) = get_shift_and_flow(Cn[0,...],Cn[1,...],dims,projection=None,plot_bool=False)
  #out = cv2.getPerspectiveTransform(Cn[0,...],Cn[1,...])
  #grid_ref = np.stack([x_grid,y_grid,np.zeros(dims)],2).reshape(-1,3)
  #grid2 = np.stack([x_grid+flow[...,0],y_grid+flow[...,1],np.zeros(dims)],2).reshape(-1,3)
  #print(shift)
  #grid_ref = np.stack([x_grid,y_grid],2).reshape(-1,2)
  #grid2 = np.stack([x_grid-shift[0]+flow[...,0],y_grid-shift[1]+flow[...,1]],2).reshape(-1,2)
  #return grid_ref
  #print(grid_ref.shape)
  #out = cv2.estimateAffine2D(grid_ref,grid2)
  #out = cv2.estimateAffine3D(grid_ref,grid2)
  
  #out = cv2.calibrateCamera(grid_ref,grid2,dims,None,None)
  print(Cn[0,...].max())
  print(Cn[1,...].max())
  plt.figure(figsize=(10,6))
  plt.subplot(231)
  plt.imshow(Cn[0],origin='lower')
  plt.title('Session 1')
  plt.subplot(234)
  plt.imshow(Cn[1],origin='lower')
  plt.title('Session 2')
  ax = plt.subplot(232)
  C = signal.convolve(Cn[0,...]-Cn[0,...].mean(),Cn[1,::-1,::-1]-Cn[1,...].mean(),mode='same')/(np.prod(dims)*Cn[0,...].std()*Cn[1,...].std())
  ax.imshow(C,origin='lower')
  ax.set_xticks(range(240,270,5))
  ax.set_xticklabels(range(-15,15,5))
  ax.set_yticks(range(240,270,5))
  ax.set_yticklabels(range(-15,15,5))
  
  #plt.set('yticks',range(240,270,7),'yticklabels',range(-15,15,7))
  ax.set_xlim([240,271])
  ax.set_ylim([240,271])
  #plt.colorbar()
  plt.title('image correlation (= %5.3g)'%C.max())
  
  plt.subplot(235)
  idxes=15
  plt.quiver(x_grid[::idxes,::idxes], y_grid[::idxes,::idxes], flow[::idxes,::idxes,0], flow[::idxes,::idxes,1], angles='xy', scale_units='xy', scale=0.5, headwidth=4,headlength=4, width=0.002, units='width')
  plt.xlim([0,512])
  plt.ylim([0,512])
  
  plt.title('optical flow')
  
  plt.subplot(233)
  
  Cn_plot = np.zeros(dims+(3,))
  Cn_plot[...,0] = Cn[0]
  Cn_plot[...,1] = Cn[1]
  Cn_plot /= Cn_plot.max()
  plt.imshow(Cn_plot,origin='lower')
  plt.title('1+2 (unaligned)')
  
  plt.subplot(236)
  x_remap = (x_grid - shift[0] + flow[:,:,0])
  y_remap = (y_grid - shift[1] + flow[:,:,1])
  Cn_new = cv2.remap(Cn[1], x_remap, y_remap, cv2.INTER_CUBIC)
  
  Cn_plot = np.zeros(dims+(3,))
  Cn_plot[...,0] = Cn[0]
  Cn_plot[...,1] = Cn_new
  Cn_plot /= Cn_plot.max()
  plt.imshow(Cn_plot,origin='lower')
  plt.title('1+2 (aligned)')
  plt.tight_layout()
  plt.show(block=False)
  svPath = pathcat([pathMouse,'alignExample.png'])
  plt.savefig(svPath,format='png',dpi=300)
  
  #plt.figure()
  #plt.imshow(Cn[0])
  #plt.show(block=False)
  
  #plt.figure()
  #plt.imshow(Cn_new)
  #plt.show(block=False)

  #grid_new = np.zeros(grid2.shape)
  
  #im_out = cv2.warpAffine(Cn[1],out[0],(Cn.shape[1:]))
  #print(Cn.shape[1:])
  #plt.figure()
  #plt.subplot(121)
  #plt.imshow(im_out)
  #plt.subplot(122)
  #plt.imshow(im_out)
  #plt.show(block=False)
  
  #shift, flow, (x_grid,y_grid) = get_shift_and_flow(im_out,Cn[1],dims,projection=None,plot_bool=False)
  
  return out#,im_out
