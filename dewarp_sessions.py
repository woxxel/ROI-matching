from scipy.io import loadmat,savemat
import numpy as np
import scipy, time, cv2

import matplotlib.pyplot as plt

from utils import pathcat, pickleData

def dewarp_sessions(basePath,mouse,sessions,use_opt_flow=True,plot_bool=False):
  
  max_thr = 0.01
  idxes = 15
  
  t_start = time.time()
  pathMouse = pathcat([basePath,mouse])
  
  nSes = sessions[-1]-sessions[0]+1
  
  pathSession = pathcat([pathMouse,'Session%02d'%sessions[0]])
  pathData = pathcat([pathSession,'results_OnACID.mat'])
  
  ld = loadmat(pathData)
  template1 = ld['Cn']
  template1 -= template1.min()
  template1 /= template1.max()
  dims = template1.shape
  
  A1 = ld['A']
  A1 = scipy.sparse.vstack([a.multiply(a>max_thr*a.max())/a.sum() for a in A1.T]).T
  if 'csc_matrix' not in str(type(A1)):
    A1 = scipy.sparse.csc_matrix(A1)
  
  pathSave = pathcat([pathSession,'footprints_dewarped'])
  #pickleData(A1,'%s.pkl'%pathSave,'save')
  savemat('%s.mat'%pathSave, {'A':A1,'Cn':template1})
  for s in range(sessions[0],nSes):
    
    pathSession = pathcat([pathMouse,'Session%02d'%(s+1)])
    pathData = pathcat([pathSession,'results_OnACID.mat'])
    
    print('Now processing %s'%pathData)
    ld = loadmat(pathData)
    template2 = ld['Cn']
    template2 -= template2.min()
    template2 /= template2.max()
    
    A2 = ld['A']
    if 'ndarray' not in str(type(A2)):
        A2 = A2.toarray()
    
    x_grid, y_grid = np.meshgrid(np.arange(0., dims[1]).astype(
                np.float32), np.arange(0., dims[0]).astype(np.float32))
    
    ## align ROIs from session 2 to the template from session 1
    print('adjust session position')
    C = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.fft2(template1) * np.fft.fft2(np.rot90(template2,2)))))
    max_pos = np.where(C==np.max(C))
    x_shift = (max_pos[1] - (dims[1]/2-1)).astype(int)
    y_shift = (max_pos[0] - (dims[0]/2-1)).astype(int)
    
    print('shift by x,y: %5.3f,%5.3f'%(x_shift,y_shift))
    x_remap = (x_grid - x_shift).astype(np.float32)
    y_remap = (y_grid - y_shift).astype(np.float32)
    
    A_2t = np.reshape(A2, dims + (-1,), order='F').transpose(2, 0, 1) # cast A2 to (x,y,n) -> (n,x,y) representation
    A2_shift = np.stack([cv2.remap(img.astype(np.float32), x_remap,
                            y_remap, cv2.INTER_NEAREST) for img in A_2t], axis=0)
    template2 = cv2.remap(template2, x_remap, y_remap, cv2.INTER_NEAREST)
    
    if use_opt_flow:    ## for each pixel, find according position in other map
      
      template1_norm = np.uint8(template1*(template1 > 0)*255)
      template2_norm = np.uint8(template2*(template2 > 0)*255)
      flow = cv2.calcOpticalFlowFarneback(np.uint8(template1_norm*255),
                                          np.uint8(template2_norm*255),
                                          None,0.5,3,128,3,7,1.5,0)
      x_remap = (flow[:,:,0] + x_grid).astype(np.float32) 
      y_remap = (flow[:,:,1] + y_grid).astype(np.float32)
      
      A2 = np.stack([cv2.remap(img.astype(np.float32), x_remap,
                              y_remap, cv2.INTER_NEAREST) for img in A2_shift], axis=0)
      template2 = cv2.remap(template2, x_remap, y_remap, cv2.INTER_NEAREST)
    else:
      A2 = A2_shift
    A2 = np.reshape(A2.transpose(1, 2, 0),                      # cast A2 to (x,y,n) again
                    (A1.shape[0], A_2t.shape[0]), order='F')    # and 
      
    print("Aligning done, t = %5.3fs"%(time.time()-t_start))
    
    #if max_thr > 0:
    A2 = scipy.sparse.csc_matrix(A2)
    A2 = scipy.sparse.vstack([a.multiply(a>max_thr*a.max())/a.sum() for a in A2.T]).T    ## normalizing
    
    pathSave = pathcat([pathSession,'footprints_dewarped'])
    #pickleData(A2,'%s.pkl'%pathSave,'save')
    savemat('%s.mat'%pathSave, {'A':A2,'Cn':template2})
    
    if plot_bool:
      plt.figure()
      plt.quiver(x_grid[::idxes,::idxes], y_grid[::idxes,::idxes], flow[::idxes,::idxes,0], flow[::idxes,::idxes,1], angles='xy', scale_units='xy', scale=0.25, headwidth=4,headlength=4, width=0.002, units='width')
      plt.show(block=False)