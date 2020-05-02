import cv2, sys, os, time, scipy
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
from tqdm import *
from scipy.io import loadmat, savemat
from utils import pathcat, find_modes, pickleData, com, get_shift_and_flow

sys.path.append('/home/wollex/Data/Science/PhD/Programs/CaImAn')

import caiman as cm
from caiman.source_extraction import cnmf as cnmf
from caiman.source_extraction.cnmf.initialization import hals

from caiman.paths import caiman_datadir
from caiman.motion_correction import MotionCorrect
from past.utils import old_div

from build_PC_cluster import *

### import packages from CNMF, etc

class test_silence_in_sessions:
  
  def __init__(self,basePath,mouse,nSes,session_classification=None,dataSet='OnACID'):
    
    #pathMouse = pathcat([basePath,mouse])
    self.dataSet = dataSet
    self.cluster = cluster(basePath,mouse,nSes,dataSet=dataSet)
    if (not os.path.exists(self.cluster.meta['svIDs'])) & (not os.path.exists(self.cluster.meta['svActivity'])):
      self.cluster.allocate_cluster()
      self.cluster.extend_dicts()
      self.cluster.get_matching()
      self.cluster.session_classification(sessions=session_classification)
      self.cluster.cluster_classification()
      self.cluster.get_PC_fields()
      self.cluster.find_PCs()
    else:
      self.cluster.load([True,False,True,False,False])
    
    self.dims = self.cluster.meta['dims']
    
  def run_all(self,s_process,n_processes,plt_bool=False):
    
    for s in s_process:
      
      if self.cluster.boolSessions[s]:
        t_start = time.time()
        self.obtain_footprints(s)
        t_end = time.time()
        print('### --------- footprints constructed - time took: %5.3fs ---------- ###'%(t_end-t_start))
        
        self.prepare_CNMF(n_processes)
        self.run_detection(n_processes,cleanup=True)
        t_end = time.time()
        print('### --------- rerun completed - time took: %5.3fs ---------- ###'%(t_end-t_start))
        
        self.analyze_traces(redo=True)
        self.save_results('mat')
        t_end = time.time()
        print('### --------- analysis and saving completed - time took: %5.3fs ---------- ###'%(t_end-t_start))
        
        if plt_bool:
          self.analyze_pre_plot(redo=True)
          self.plot_analysis()
          plt.close('all')
      
    
  def obtain_footprints(self,s,max_diff=None):
    
    print('Run redetection of "silent" cells on session %d'%(s+1))
    
    self.s = s
    self.idxes = {}
    self.dataIn = {}
    self.dataOut = {}
    self.data = {}
    
    ## find silent neurons in session s
    self.idxes['silent']= ~self.cluster.activity['status'][:,self.s,1]
    
    if max_diff is None:
      max_diff = self.cluster.meta['nSes']
    
    dims = self.dims
    T = 8989
    self.dataIn = {}
    self.dataIn['A'] = scipy.sparse.csc_matrix((np.prod(self.dims),self.cluster.meta['nC']))
    self.dataIn['C'] = np.random.rand(self.cluster.meta['nC'],T)
    
    ## load footprints of active cells
    self.pathSession = pathcat([self.cluster.meta['pathMouse'],'Session%02d'%(s+1)])
    pathData = pathcat([self.pathSession,'results_OnACID.mat'])
    ld = loadmat(pathData,variable_names=['A','C','b','f'])
    
    c_idx = np.where(~self.idxes['silent'])[0]
    n_idx = self.cluster.IDs['neuronID'][c_idx,s,1].astype('int')
    
    self.dataIn['A'][:,c_idx] = scipy.sparse.vstack([a/a.sum() for a in ld['A'][:,n_idx].T]).T
    self.dataIn['C'][c_idx,:] = ld['C'][n_idx,:]
    self.dataIn['b'] = ld['b'].T
    self.dataIn['f'] = ld['f'].T
    ## construct footprint of silent cells (interpolation between adjacent sessions?) and adjust for shift & rotation
    s_ref = np.zeros((self.cluster.meta['nC'],2))*np.NaN
    n_ref = np.zeros((self.cluster.meta['nC'],2))*np.NaN
    
    for c in np.where(self.idxes['silent'])[0]:
      s_pre = np.where(self.cluster.activity['status'][c,:s,1])[0]
      s_post = s+1+np.where(self.cluster.activity['status'][c,s+1:,1])[0]
      
      if len(s_pre)>0:
        s_ref[c,0] = s_pre[-1] if (s-s_pre[-1]) <= max_diff else np.NaN
        if not np.isnan(s_ref[c,0]):
          n_ref[c,0] = self.cluster.IDs['neuronID'][c,int(s_ref[c,0]),1]
      
      if len(s_post)>0:
        s_ref[c,1] = s_post[0] if (s_post[0] - s) <= max_diff else np.NaN
        if not np.isnan(s_ref[c,1]):
          n_ref[c,1] = self.cluster.IDs['neuronID'][c,int(s_ref[c,1]),1]
    
    s_load = np.unique(s_ref[~np.isnan(s_ref)])
    
    print('loading footprints...')
    
    progress = tqdm(zip(s_load.astype('int'),range(len(s_load))),total=len(s_load))
    for (s_ld,i) in progress:
      
      pathSession = pathcat([self.cluster.meta['pathMouse'],'Session%02d'%(s_ld+1)])
      pathData = pathcat([pathSession,'results_OnACID.mat'])
      ld = loadmat(pathData,variable_names=['A'])
      
      A_tmp = ld['A'].tocsc()
      
      (x_shift,y_shift),flow,(x_grid,y_grid),_ = get_shift_and_flow(self.dataIn['A'],A_tmp,dims,projection=1,plot_bool=False)
        
      x_remap = (x_grid - x_shift + flow[:,:,0])
      y_remap = (y_grid - y_shift + flow[:,:,1])
      
      A_tmp = sp.sparse.hstack([sp.sparse.csc_matrix(cv2.remap(img.reshape(self.dims), x_remap,y_remap, cv2.INTER_CUBIC).reshape(-1,1)) for img in A_tmp[:,n_ref[s_ref==s_ld].astype('int')].toarray().T])
      
      self.dataIn['A'][:,np.where(s_ref==s_ld)[0]] += scipy.sparse.vstack([a/a.sum() for a in A_tmp.T]).T
      #print('%d neuron footprints taken from session %d'%(A_tmp.shape[1],s_ld+1))
      
    
    max_thr = 0.001
    self.dataIn['A'] = scipy.sparse.vstack([a.multiply(a>(max_thr*a.max()))/a.sum() for a in self.dataIn['A'].T]).T
    
  def prepare_CNMF(self,n_processes):
    
    fname = None
    for f in os.listdir(self.pathSession):
      if f.startswith("thy") or f.startswith("shank"):
        fname = pathcat([self.pathSession,f])
        if f.endswith('.h5'):
          break
    
    if not fname or not os.path.exists(fname):
      print("No file here to process :(")
      return
    
    #sv_dir = "/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/tmp"
    sv_dir = "/home/wollex/Data/Science/PhD/Data/tmp"
    #svname = self.pathSession + "results_OnACID.mat"
    #if os.path.exists(svname):
      #print("Processed file already present - skipping")
      #return
    
    svname_h5 = self.pathSession + "results_OnACID.hdf5"
    
    t_start = time.time()
    
    # set up CNMF parameters
    params_dict ={
            
            #general data
            'fnames': [fname],
            'fr': 15,
            'decay_time': 0.47,
            
            # init
            'ssub': 2,                          # spatial subsampling during initialization
            'tsub': 5,                          # temporal subsampling during initialization
            
            #online
            'init_batch': 300,                  # number of frames for initialization
            
            #preprocess
            #'compute_g': True,
            
            #quality
            'min_SNR': 2,                     # minimum SNR for accepting candidate components
            'SNR_lowest': 0.5,                  # minimum accepted SNR value
            'rval_thr': 0.85,                   # space correlation threshold for accepting a component
            'rval_lowest': 0,                   
            #'test_both': True,                 # use CNN and correlation to test for new components
            'use_cnn': True
    }
    
    self.opts = cnmf.params.CNMFParams(params_dict=params_dict)
    
    self.fname_memmap = pathcat([sv_dir,'memmap%s_%d__d1_512_d2_512_d3_1_order_C_frames_8989_.mmap'%(self.cluster.meta['mouse'],(self.s+1))])
    print(self.fname_memmap)
    if not os.path.exists(self.fname_memmap):
      print("Start writing memmapped file @t = " +  time.ctime())
      if n_processes > 1:
        c, self.dview, n_processes = cm.cluster.setup_cluster(backend='local', n_processes=n_processes, single_thread=False)
      else:
        self.dview=None
      
      self.fname_memmap = cm.save_memmap([fname], base_name='memmap%s_%d_'%(self.cluster.meta['mouse'],(self.s+1)), save_dir=sv_dir, n_chunks=100, order='C', dview=self.dview)  # exclude borders
      if n_processes > 1:
          cm.stop_server(dview=self.dview)      ## restart server to clean up memory
    #self.preprocessed = False
      
  
  def run_detection(self,n_processes,cleanup=True):
    
    #if not self.preprocessed:
    self.cnm = cnmf.CNMF(n_processes)
    
    Yr, dims, T = cm.load_memmap(self.fname_memmap)
    images = np.reshape(Yr.T, [T] + list(dims), order='F')
    
    Y = np.reshape(Yr.T, [T] + list(self.cluster.meta['dims']), order='F')
    self.cnm.estimates.dims = self.cnm.dims = dims
    self.idxes['in'] = np.array(self.dataIn['A'].sum(0)>0).ravel()
    self.cnm.estimates.A = self.dataIn['A'][:,self.idxes['in']]
    
    self.cnm.estimates.C = self.dataIn['C'][self.idxes['in'],:]
    self.cnm.estimates.b = self.dataIn['b'].T
    self.cnm.estimates.f = self.dataIn['f'].T
    
    self.cnm.params.change_params({'p':1})
    
    self.cnm.params.set('online', {'init_batch': 2000})
    Yr = np.transpose(np.reshape(images, (T, -1), order='F'))
    if np.isfortran(Yr):
        raise Exception('The file is in F order, it should be in C order (see save_memmap function')
    
    # Make sure filename is pointed correctly (numpy sets it to None sometimes)
    try:
        Y.filename = images.filename
        Yr.filename = images.filename
    except AttributeError:  # if no memmapping cause working with small data
        pass
    
    if self.cnm.estimates.sn is None:
      print('preprocessing...')
      Yr = self.cnm.preprocess(Yr)
      self.preprocessed = True
      print('done!')
    
    print('temporal and spatial update...')
    self.cnm.update_temporal(Yr)
    self.cnm.update_spatial(Yr)
    self.cnm.update_temporal(Yr)
    print('done!')
    
    self.cnm.deconvolve(optimize_g=5,s_min=0)
    self.cnm.estimates.evaluate_components(Y,self.opts)
    
    if cleanup:
      os.remove(self.fname_memmap)
    
    
  ## call cnmf to run 
  ###   1. temporal trace updates on ROIs and background
  ###   2. spatial update on silent neurons (?) - rather not!
  
  ## evaluate results: 
  ### significant calcium-activity? / spikes?
  ### correlation with overlapping/nearby neurons?
  
  def run_evaluation(self):
    Yr, dims, T = cm.load_memmap(self.fname_memmap)
    images = np.reshape(Yr.T, [T] + list(dims), order='F')
    
    Y = np.reshape(Yr.T, [T] + list(self.cluster.meta['dims']), order='F')
    Yr = np.transpose(np.reshape(images, (T, -1), order='F'))
    #print(Y.shape)
    #print(Yr.shape)
    
    #self.Y_proj = np.sum(Y,0)
    #self.Yr_proj = np.sum(Yr,1).reshape(dims)
    
    #plt.figure()
    #plt.subplot(121)
    #plt.imshow(self.Yr_proj,origin='lower')
    
    #plt.subplot(122)
    #plt.imshow(self.Yr_proj,origin='lower')
    #plt.show(block=False)
    
    
    #plt.figure()
    #plt.subplot(121)
    #plt.imshow(self.dataOut['A'].sum(1).reshape(dims),origin='lower')
    
    #plt.subplot(122)
    #plt.imshow(self.dataOut['A'].sum(1).reshape(dims),origin='lower')
    #plt.show(block=False)
    
    #plt.pause(1)
    
    #self.cnm.estimates.dims = self.cnm.dims = dims
    self.cnm.deconvolve(optimize_g=5,s_min=0)
    self.cnm.estimates.evaluate_components(Y,self.opts)
    
    self.analyze_traces(redo=True)
    #self.analyze_pre_plot(redo=True)
    #self.plot_analysis()
    
    
  
  def analyze_traces(self,redo=False):
    
    dims = self.cluster.meta['dims']
    nC = self.cluster.meta['nC']
    
    
    self.idxes['in'] = np.array(self.dataIn['A'].sum(0)>0).ravel()
    idxes_Ain = np.where(self.idxes['in'])[0]
    
    ### find proper IDs --- some are removed, some are weirdly remapped to completely other places
    print('Finding proper IDs...')
    nC_in = self.idxes['in'].sum()
    nC_out = self.cnm.estimates.A.shape[1]
    if nC_out < nC_in:
      print('\t %d neurons were removed ...'%(nC_in-nC_out))
    Aout = self.cnm.estimates.A
    
    ## find next footprints with low distance
    cm_in = com(self.dataIn['A'][:,self.idxes['in']],dims[0],dims[1])
    cm_out = com(Aout,dims[0],dims[1])
    
    idx_out_bool = np.zeros(nC_out,'bool')
    idx_in_bool = np.zeros(nC_in,'bool')
    
    idx_new = []
    dc = 0              ## difference in index
    for c in range(nC_out):
      d = np.linalg.norm(cm_in[c+dc,:]-cm_out[c,:])
      if d > 10:
        dist = np.linalg.norm(cm_in-cm_out[c,:],axis=1)
        dc_tmp = np.argmin(dist)-c
        if (dc_tmp>=dc) & (dc_tmp<=(nC_in-nC_out)):
          dc = dc_tmp
          idx_in_bool[c+dc] = True
          idx_out_bool[c] = True
          
          print('new match has been found: %d to %d with dc=%d'%(c,c+dc,dc))
        else:
          idx_new.append(c)
          print('no match could be found for neuron %d, minimum distance: %5.3g, dc = %d'%(c,dist.min(),dc_tmp))
      else:
        idx_in_bool[c+dc] = True
        idx_out_bool[c] = True
    
    nC_new = len(idx_new)
    
    print('%d new neurons detected'%nC_new)
    self.idxes['new'] = np.hstack([np.zeros(nC,'bool'),np.ones(nC_new,'bool')])
    self.idxes['previous'] = np.hstack([~self.idxes['silent'],np.zeros(nC_new,'bool')])
    self.idxes['in'] = np.zeros(nC+nC_new,'bool')
    self.idxes['in'][idxes_Ain[idx_in_bool]] = True
    
    print('done!')
    
    idxes_Ain = np.where(self.idxes['in'])[0]
    
    idx_out_evaluate = np.zeros(nC_out,'bool')
    idx_out_evaluate[self.cnm.estimates.idx_components] = True
    idx_out_evaluate_new = idx_out_evaluate[~idx_out_bool]
    idx_out_evaluate_old = idx_out_evaluate[idx_out_bool]
    
    
    idx_evaluate = np.zeros(nC+nC_new,'bool')
    idx_evaluate[idxes_Ain[idx_out_evaluate_old]] = True
    idx_evaluate[nC:] = idx_out_evaluate_new
    
    self.dataOut = {}
    self.idxes['evaluate'] = idx_evaluate
    
    idx_bad = ~idx_evaluate
    idx_good = idx_evaluate
    
    self.idxes['active_good'] = self.idxes['previous'] & idx_good
    self.idxes['active_bad'] = self.idxes['previous'] & idx_bad
    
    self.idxes['silent_good'] = ~self.idxes['previous'] & idx_good
    self.idxes['silent_bad'] = ~self.idxes['previous'] & idx_bad
    
    print('active good: %d, \t active bad: %d'%(self.idxes['active_good'].sum(),self.idxes['active_bad'].sum()))
    print('silent good: %d, \t silent bad: %d'%(self.idxes['silent_good'].sum(),self.idxes['silent_bad'].sum()))
    
    
    self.dataOut['A'] = scipy.sparse.csc_matrix((np.prod(dims),nC+nC_new))
    self.dataOut['C'] = np.zeros((nC+nC_new,8989))
    self.dataOut['S'] = np.zeros((nC+nC_new,8989))
    
    self.dataOut['SNR'] = np.zeros(nC+nC_new)
    self.dataOut['r_values'] = np.zeros(nC+nC_new)
    self.dataOut['CNN'] = np.zeros(nC+nC_new)
    
    
    self.dataOut['A'][:,self.idxes['in']] = self.cnm.estimates.A.tocsc()[:,idx_out_bool]
    self.dataOut['A'][:,nC:] = self.cnm.estimates.A.tocsc()[:,~idx_out_bool]
    
    self.dataOut['C'][self.idxes['in'],:] = self.cnm.estimates.C[idx_out_bool,:]
    self.dataOut['C'][nC:,:] = self.cnm.estimates.C[~idx_out_bool,:]
    
    self.dataOut['S'][self.idxes['in'],:] = self.cnm.estimates.S[idx_out_bool,:]
    self.dataOut['S'][nC:,:] = self.cnm.estimates.S[~idx_out_bool,:]
    
    self.dataOut['SNR'][self.idxes['in']] = self.cnm.estimates.SNR_comp[idx_out_bool]
    self.dataOut['SNR'][nC:] = self.cnm.estimates.SNR_comp[~idx_out_bool]
    
    self.dataOut['r_values'][self.idxes['in']] = self.cnm.estimates.r_values[idx_out_bool]
    self.dataOut['r_values'][nC:] = self.cnm.estimates.r_values[~idx_out_bool]
    
    self.dataOut['CNN'][self.idxes['in']] = self.cnm.estimates.cnn_preds[idx_out_bool]
    self.dataOut['CNN'][nC:] = self.cnm.estimates.cnn_preds[~idx_out_bool]
    
    ### analyze active neurons: did they maintain same behavior? (correlation > 0.9)
    #self.Cout[self.idxes['in'],:] = np.vstack([Ca/Ca.max() for Ca in self.cnm.estimates.C])
    
    #t_start = time.time()
    #self.dataOut['fitness'] = np.zeros(nC+nC_new)
    #self.dataOut['fitness'][self.idxes['in']], erfc, sd_r, md = compute_event_exceptionality(self.dataOut['C'][self.idxes['in'],:])
    #if nC_new > 0:
      #self.dataOut['fitness'][nC:], erfc, sd_r, md = compute_event_exceptionality(self.dataOut['C'][nC:,:])
    #self.fitness[self.fitness<-10**4] = -10**4
    
    #t_end = time.time()
    #print('fitness computed - time took: %5.3g'%(t_end-t_start))
    
    
    
    
  def analyze_pre_plot(self,redo=False):
    
    idx_good = self.idxes['evaluate']
    idx_bad = ~self.idxes['evaluate']
    
    self.idxes['active_good'] = self.idxes['previous'] & idx_good
    self.idxes['active_bad'] = self.idxes['previous'] & idx_bad
    
    self.idxes['silent_good'] = ~self.idxes['previous'] & idx_good
    self.idxes['silent_bad'] = ~self.idxes['previous'] & idx_bad
    
    print('active good: %d, \t active bad: %d'%(self.idxes['active_good'].sum(),self.idxes['active_bad'].sum()))
    print('silent good: %d, \t silent bad: %d'%(self.idxes['silent_good'].sum(),self.idxes['silent_bad'].sum()))
    
    nC = self.idxes['previous'].shape[0]
    dims = self.cluster.meta['dims']
    
    self.dataOut['C_std'] = self.dataOut['C'].std(1)
    self.dataIn['C_std'] = self.dataIn['C'].std(1)
    
    self.data['nu'] = np.zeros(nC)*np.NaN
    for c in tqdm(range(nC),leave=False):
      if (self.dataOut['S'][c,:]>0).sum():
        spike_nr, md, sd_r = get_spikeNr(self.dataOut['S'][c,self.dataOut['S'][c,:]>0])
      else:
        spike_nr = 0
      self.data['nu'][c] = spike_nr/(8989/15)
    
    
    #if redo | (not hasattr(self,'A_norm_out')):
    #t_start = time.time()
    
    self.data['cm'] = com(self.dataOut['A'],dims[0],dims[1])
    self.data['cm'][np.all(self.data['cm']==0,1),:] = np.NaN
    self.data['D_ROIs'] = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(self.data['cm']))
    np.fill_diagonal(self.data['D_ROIs'],np.NaN)
    
    self.dataOut['A_norm'] = np.zeros(nC)*np.NaN
    self.dataIn['A_norm'] = np.zeros(nC)*np.NaN
    chunks = 500
    idx_Ain = np.where(self.idxes['in'])[0]
    nC_in = len(idx_Ain)
    steps = int(nC_in/chunks)
    for i in tqdm(range(steps+1),leave=False):
      c_start = i*chunks
      c_end = min(nC_in,(i+1)*chunks)
      self.dataOut['A_norm'][idx_Ain[c_start:c_end]] = np.linalg.norm(self.dataOut['A'][:,idx_Ain[c_start:c_end]].toarray(),axis=0)
      self.dataIn['A_norm'][idx_Ain[c_start:c_end]] = np.linalg.norm(self.dataIn['A'][:,idx_Ain[c_start:c_end]].toarray(),axis=0)
    #t_end = time.time()
    #print('Anorm computed - time took: %5.3g'%(t_end-t_start))
    
    #if redo | (not hasattr(self,'fp_corr')):
    self.data['corr'] = scipy.sparse.csc_matrix((nC,nC))#np.zeros(nC)*np.NaN
    self.data['fp_corr'] = scipy.sparse.csc_matrix((nC,nC))
    for c in tqdm(range(nC),leave=False):
      
      if self.idxes['previous'][c]:
        self.data['fp_corr'][c,c] = self.dataOut['A'][:,c].multiply(self.dataIn['A'][:,c]).sum()/(self.dataOut['A_norm'][c]*self.dataIn['A_norm'][c])
        self.data['corr'][c,c] = np.cov(self.dataIn['C'][c,:],self.dataOut['C'][c,:])[0,1]/(self.dataIn['C_std'][c]*self.dataOut['C_std'][c])
      
      idx_close = np.where(self.data['D_ROIs'][c,:]<10)[0]
      if len(idx_close)>0:
        for cc in idx_close:
          self.data['fp_corr'][c,cc] = self.dataOut['A'][:,c].multiply(self.dataOut['A'][:,cc]).sum()/(self.dataOut['A_norm'][c]*self.dataOut['A_norm'][cc])
          
          if self.idxes['previous'][cc]:
          #corr_tmp = np.cov(self.Cout[c,:],self.Cout[cc,:])[0,1]/(self.dataOut['C_std'][c]*self.dataOut['C_std'][cc])
            self.data['corr'][c,cc] = np.cov(self.dataOut['C'][c,:],self.dataIn['C'][cc,:])[0,1]/(self.dataOut['C_std'][c]*self.dataIn['C_std'][cc])
          
    #self.data['fp_corr'].tocsc()
    
  
  def save_results(self,ext='mat'):
    
    ### save everything important from running the detection
    #idx_evaluate = np.zeros(self.cluster.meta['nC']).astype('bool')
    #idx_evaluate[self.idxes['in']] = True
    #idx_evaluate[self.cnm.estimates.idx_components] = True
    
    results = {'A':self.dataOut['A'],
               'Ain':self.dataIn['A'],
               'C':self.dataOut['C'],
               'S':self.dataOut['S'],
               'b':self.cnm.estimates.b,
               'f':self.cnm.estimates.f,
               'Cin':self.dataIn['C'],
               'idx_previous':self.idxes['previous'],
               'idx_Ain':self.idxes['in'],
               'idx_evaluate':self.idxes['evaluate'],
               'SNR':self.dataOut['SNR'],
               'r_values':self.dataOut['r_values'],
               'CNN':self.dataOut['CNN']}
               #'fitness':self.dataOut['fitness']}
    
    svPath = pathcat([self.cluster.meta['pathMouse'],'Session%02d'%(self.s+1),'results_redetect.%s'%ext])
    if ext == 'mat':
      savemat(svPath,results)
      print('Data saved in %s'%svPath)
    else:
      pickleData(results,svPath,'save')
    #something from evaluation
    
    ### analyze silent ones: do they show firing behavior at all? are they correlated with nearby active neurons?
  
  def load(self,s,ext='mat'):
    
    self.s = s-1
    pathSession = pathcat([self.cluster.meta['pathMouse'],'Session%02d'%s])
    pathData = pathcat([pathSession,'results_redetect.%s'%ext])
    #pathData = pathcat([pathSession,'results_postSilent.%s'%ext])
    print('loading data from %s'%pathData)
    if ext == 'mat':
      ld = loadmat(pathData,squeeze_me=True)
    else:
      ld = pickleData([],pathData,'load')
    self.dataOut = {'A':ld['A'],
                    'C':ld['C'],
                    'S':ld['S'],
                    'SNR':ld['SNR'],
                    'r_values':ld['r_values'],
                    'CNN':ld['CNN']}
                    #'fitness':ld['fitness']}
    
    self.dataIn = {'A':ld['Ain'],
                   #'Cn':ld['Cn'],
                   'C':ld['Cin']}
    
    self.idxes = {'previous':ld['idx_previous'].astype('bool'),
                  'in':ld['idx_Ain'].astype('bool'),
                  'evaluate':ld['idx_evaluate'].astype('bool')}
    
    self.data = {}
    
    
  def plot_analysis(self,SNR_thr=2,r_val_thr=0,CNN_thr=0.8):
    
    
    idx_good = (self.dataOut['SNR']>SNR_thr) & (self.dataOut['r_values']>r_val_thr) & (self.dataOut['CNN']>CNN_thr)
    idx_bad = ~idx_good
    
    self.idxes['active_good'] = self.idxes['previous'] & idx_good
    self.idxes['active_bad'] = self.idxes['previous'] & idx_bad
    
    self.idxes['silent_good'] = ~self.idxes['previous'] & idx_good
    self.idxes['silent_bad'] = ~self.idxes['previous'] & idx_bad
    
    print('active good: %d, \t active bad: %d'%(self.idxes['active_good'].sum(),self.idxes['active_bad'].sum()))
    print('silent good: %d, \t silent bad: %d'%(self.idxes['silent_good'].sum(),self.idxes['silent_bad'].sum()))
    
    
    idxs = [self.idxes['active_good'],self.idxes['active_bad'],self.idxes['silent_good'],self.idxes['silent_bad']]
    col_arr = ['k','b','r','g']
    
    ### plot changes from input to output
    acorr = np.diag(self.data['corr'].toarray())
    fp_acorr = np.diag(self.data['fp_corr'].toarray())
    
    corr = np.copy(self.data['corr'].toarray())
    fp_corr = np.copy(self.data['fp_corr'].toarray())
    
    np.fill_diagonal(corr,np.NaN)
    np.fill_diagonal(fp_corr,np.NaN)
    
    max_fp_corr = np.nanmax(fp_corr,1)
    max_Ca_corr = np.nanmax(corr,1)
    
    al=0.5
    
    plt.figure(figsize=(10,4))
    ax = plt.axes([0.39,0.55,0.1,0.19])
    for i in range(2):
      ax.hist(acorr[idxs[i]],np.linspace(0.5,1,21),facecolor=col_arr[i],alpha=al,orientation='horizontal')
    #ax.xlabel('$C_{Ca^{2+}}^{in-out}$',fontsize=14)
    ax.invert_xaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax = plt.axes([0.5,0.75,0.1,0.19])
    for i in range(4):
      #if i==1:
        #continue
      plt.hist(fp_acorr[idxs[i]],np.linspace(0.5,1,21),facecolor=col_arr[i],alpha=al)
    #plt.xlabel('$C_{fp}^{in-out}$',fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax = plt.axes([0.5,0.55,0.1,0.19])
    for i in range(4):
      ax.scatter(fp_acorr[idxs[i]],acorr[idxs[i]],s=2,c=col_arr[i])  
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_xlabel('$C_{fp}$',fontsize=14)
    ax.set_ylabel('$C_{Ca^{2+}}$',fontsize=14)
    #plt.tight_layout()
    #plt.show(block=False)
    #ext = 'png'
    #pathSession = pathcat([self.cluster.meta['pathMouse'],'Session%02d'%(self.s+1)])
    #pathFigure = pathcat([pathSession,'find_silent1.%s'%(ext)]);
    #plt.savefig(pathFigure,format=ext,dpi=300)
    #print('Figure saved as %s'%pathFigure)
    
    ### plot goodness of detected ROIs
    #plt.figure(figsize=(10,5))
    
    plt.subplot(231)
    plt.hist(self.data['D_ROIs'].flatten(),np.linspace(0,700,101),facecolor='k')
    
    ax_zoom = plt.axes([0.2,0.8,0.15,0.15])
    ax_zoom.hist(self.data['D_ROIs'].flatten(),np.linspace(0,20,41),facecolor='k')
    
    #ax = plt.subplot(232)
    #for i in range(4):
      #ax.hist(-self.dataOut['fitness'][idxs[i]],np.logspace(0.5,4,21),facecolor=col_arr[i],alpha=al)
    #ax.set_xscale('log')
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #ax.set_xlim([10**0.5,10**4])
    #ax.set_xlabel('fitness',fontsize=14)
    
    ax = plt.subplot(233)
    for i in range(4):
      ax.hist(self.dataOut['r_values'][idxs[i]],np.linspace(-1,1,21),facecolor=col_arr[i],alpha=al)
    #plt.xscale('log')
    ax.set_xlim([-1,1])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('r value',fontsize=14)
    
    ax = plt.subplot(235)
    for i in range(4):
      ax.hist(self.dataOut['SNR'][idxs[i]],np.linspace(0,30,21),facecolor=col_arr[i],alpha=al,orientation='horizontal')
    ax.invert_xaxis()
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('SNR',fontsize=14)
    ax.set_ylim([0,30])
    ax.set_xlabel('occurence',fontsize=14)
    
    
    #ax = plt.subplot(235)
    #for i in range(4):
      #ax.scatter(-self.dataOut['fitness'][idxs[i]],self.dataOut['SNR'][idxs[i]],s=2,c=col_arr[i])
    #ax.set_xscale('log')
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    #ax.set_yticks([])
    #ax.set_xlabel('fitness',fontsize=14)
    #ax.set_xlim([10**(0.5),10**4.1])
    #ax.set_ylim([0,30])
    
    ax =plt.subplot(236)
    for i in range(4):
      ax.scatter(self.dataOut['r_values'][idxs[i]],self.dataOut['SNR'][idxs[i]],s=2,c=col_arr[i])
    ax.set_xlabel('r value',fontsize=14)
    #plt.ylabel('SNR',fontsize=14)
    ax.set_xlim([-1,1])
    ax.set_ylim([0,30])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    
    plt.tight_layout()
    plt.show(block=False)
    ext = 'png'
    pathSession = pathcat([self.cluster.meta['pathMouse'],'Session%02d'%(self.s+1)])
    pathFigure = pathcat([pathSession,'find_silent.%s'%(ext)]);
    plt.savefig(pathFigure,format=ext,dpi=300)
    print('Figure saved as %s'%pathFigure)
    
  def plot_traces(self,idxs,n=3):
    
    ## some nice plots:
    ###   1. active vs silent ROI shapes
    ###   2. correlation vs activity
    
    plt.figure()
    idxs = np.where(idxs)[0]
    for i in range(n):
      
      c = idxs[np.random.randint(len(idxs))]
      print(c)
      plt.subplot(2,n,i+1)
      plt.imshow(self.dataIn['A'][:,c].reshape(512,512).toarray())
      plt.ylabel('$C_{out}$')
      
      plt.subplot(4,n,2*n+i+1)
      plt.plot(self.dataOut['C'][c,:])
      plt.ylabel('$C_{out}$')
      
      plt.subplot(4,n,3*n+i+1)
      plt.plot(self.dataOut['C'][c,:])
      plt.ylabel('$C_{in}$')
    
    #plt.subplot(3,2,6)
    #plt.plot(self.cnm.estimates.S[c,:])
    #plt.ylabel('$S_{out}$')
    
    #plt.subplot(315)
    #plt.plot(Cout[c,:])
    #if not self.idxes['silent'][c]:
      #plt.plot(Cin[c,:])
    plt.suptitle('correlation: %5.3f'%self.data['corr'][c,c])
    plt.show(block=False)
    
    
    
  

def compute_event_exceptionality(traces, robust_std=False, N=5, sigma_factor=3.):
    """
    Define a metric and order components according to the probability of some "exceptional events" (like a spike).

    Such probability is defined as the likelihood of observing the actual trace value over N samples given an estimated noise distribution.
    The function first estimates the noise distribution by considering the dispersion around the mode.
    This is done only using values lower than the mode. The estimation of the noise std is made robust by using the approximation std=iqr/1.349.
    Then, the probability of having N consecutive events is estimated.
    This probability is used to order the components.

    Args:
        traces: ndarray
            Fluorescence traces

        N: int
            N number of consecutive events

        sigma_factor: float
            multiplicative factor for noise estimate (added for backwards compatibility)

    Returns:
        fitness: ndarray
            value estimate of the quality of components (the lesser the better)

        erfc: ndarray
            probability at each time step of observing the N consequtive actual trace values given the distribution of noise

        noise_est: ndarray
            the components ordered according to the fitness
    """

    T = np.shape(traces)[-1]
    
    md = find_modes(traces,axis=1)
    ff1 = traces - md[:,None]

    # only consider values under the mode to determine the noise standard deviation
    ff1 = -ff1 * (ff1 < 0)
    if robust_std:

        # compute 25 percentile
        ff1 = np.sort(ff1, axis=1)
        ff1[ff1 == 0] = np.nan
        Ns = np.round(np.sum(ff1 > 0,1) * .5)
        iqr_h = np.zeros(traces.shape[0])

        for idx, _ in enumerate(ff1):
            iqr_h[idx] = ff1[idx, -Ns[idx]]

        # approximate standard deviation as iqr/1.349
        sd_r = 2 * iqr_h / 1.349

    else:
        Ns = np.sum(ff1 > 0, -1)
        sd_r = np.sqrt(old_div(np.sum(ff1**2, -1), Ns))

    # compute z value
    z = old_div((traces - md[:,None]), (sigma_factor * sd_r[:,None]))

    # probability of observing values larger or equal to z given normal
    # distribution with mean md and std sd_r
    #erf = 1 - norm.cdf(z)

    # use logarithm so that multiplication becomes sum
    #erf = np.log(erf)
    # compute with this numerically stable function
    erf = scipy.special.log_ndtr(-z)

    # moving sum
    erfc = np.cumsum(erf, 1)
    erfc[:, N:] -= erfc[:, :-N]

    # select the maximum value of such probability for each trace
    fitness = np.min(erfc, 1)

    return fitness, erfc, sd_r, md

def get_spikeNr(traces):
    md = find_modes(traces)
    ff1 = traces - md

    # only consider values under the mode to determine the noise standard deviation
    ff1 = -ff1 * (ff1 < 0)

    # compute 25 percentile
    ff1 = np.sort(ff1)
    ff1[ff1 == 0] = np.nan
    Ns = round((ff1>0).sum() * .5).astype('int')
    iqr_h = np.zeros(traces.shape[0])

    #for idx, _ in enumerate(ff1):
    iqr_h = ff1[-Ns]

    # approximate standard deviation as iqr/1.349
    sd_r = 2 * iqr_h / 1.349
    data_thr = md+2*sd_r;
    spikeNr = np.floor(traces/data_thr).sum();
    return spikeNr,md,sd_r
