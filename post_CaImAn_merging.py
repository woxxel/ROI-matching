
import scipy
import os
import caiman
import numpy as np
from scipy.io import loadmat
from scipy.io import savemat
import matplotlib.pyplot as plt


def post_CaImAn_merging(pathMouse,sessions=None,rm_list_in=[]):
  
  ### check for pairs of componentents, that have spatial & temporal correlation > 0.5 and remove the one with less SNR
  
  for s in range(sessions[0],sessions[-1]+1):
    # load data from file
    pathData = '%sSession%02d/results_OnACID.mat'%(pathMouse,s)
    if not os.path.exists(pathData):
      print("File %s not found. skipping this session..."%pathData)
      continue
    else:
      print("Now processing file %s..."%pathData)
    
    f_ld = loadmat(pathData)
    
    #calculate center for all
    A = f_ld['A']
    C = f_ld['C']
    S = f_ld['S']
    
    if not rm_list_in:
      dims = f_ld['Cn'].shape
      cmA = com(A, dims[0], dims[1])
      
      print('ROIs found: %d'%A.shape[1])
      #find neurons with distance < 10px
      d_cmA = scipy.spatial.distance.pdist(cmA)
      d_cmA_sq = scipy.spatial.distance.squareform(d_cmA)
      d_cmA_bool = scipy.spatial.distance.squareform(d_cmA < 10)
      idxes = np.where(d_cmA_bool)
      
      # calculate spatial correlation for those
      normA = scipy.linalg.norm(A.todense(),axis=0)
      rm_list = []
      for ID1,ID2 in zip(idxes[0],idxes[1]):
        if ID1 > ID2:
          sp_corr = (A[:,ID1].T.dot(A[:,ID2])/(normA[ID1]*normA[ID2])).todense()
          
          if sp_corr > 0.5:           # for those with sp_corr > 0.5, calculate temporal correlation (or vice versa)
            tmp_corr = np.corrcoef(C[:,ID1],C[:,ID2])[0,1]
            
            if tmp_corr > 0.5:        # if tmp_corr & sp_corr > 0.5, obtain SNR of both
              print('(%04d,%04d): distance: %5.3g, \t sp_corr=%5.3g, \t tmp_corr=%5.3g'%(ID1,ID2,d_cmA_sq[ID1,ID2],sp_corr,tmp_corr))
              # remove the one with lower SNR
              fitness_raw,erc,sd_r,md = caiman.components_evaluation.compute_event_exceptionality(C[:,[ID1,ID2]].T)
              print("fitness:")
              print(fitness_raw)
              if fitness_raw[0] < fitness_raw[1]:
                rm_list.append(ID2)
              else:
                rm_list.append(ID1)
    else:
      rm_list = rm_list_in
      
      #print(rm_list)
    if len(rm_list)>0:
      print('removing %d ROIs'%len(rm_list))
      mask = np.ones(A.shape[1], dtype=bool)
      mask[rm_list] = False
      A = A[:,mask]
      if C.shape[0] > 4000:
        C = np.delete(C, rm_list, axis=1)
        S = np.delete(S, rm_list, axis=1)
      else:
        C = np.delete(C, rm_list, axis=0)
        S = np.delete(S, rm_list, axis=0)
      
      #save results under old name
      results = dict(A=A,
                    C=C,
                    S=S,
                    Cn=f_ld['Cn'],
                    b=f_ld['b'],
                    f=f_ld['f'])
      savemat(pathData,results)
        
        
    #if C.shape[1] == 8988:
      
      
      #C = np.insert(C,rm_list[0]+1,C[:,rm_list[0]],axis=1)
      #S = np.insert(S,rm_list[0]+1,S[:,rm_list[0]],axis=1)
      
      #results = dict(A=A,
                    #C=C,
                    #S=S,
                    #Cn=f_ld['Cn'],
                    #b=f_ld['b'],
                    #f=f_ld['f'])
      #savemat(pathData,results)
      
    #print(C.shape)
    #print(S.shape)
  
  
  
  
  
  
  
  
def com(A,d1,d2):
  ## from "rois.py" of CaImAn
  
  Coor = np.matrix([np.outer(np.ones(d2), np.arange(d1)).ravel(),
                          np.outer(np.arange(d2), np.ones(d1)).ravel()], dtype=A.dtype)
  cm = (Coor * A / A.sum(axis=0)).T
  return np.array(cm)