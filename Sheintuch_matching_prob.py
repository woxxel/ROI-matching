import os, sys, cv2, math, time, warnings, logging
from tqdm import *
from itertools import chain
import numpy as np
from scipy.io import loadmat
import scipy as sp
from scipy.optimize import curve_fit, linear_sum_assignment
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm

from mpl_toolkits.mplot3d import Axes3D


sys.path.append('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Programs/PC_analysis')
from utils import com, pathcat, calculate_img_correlation, get_shift_and_flow, fun_wrapper, pickleData

warnings.filterwarnings("ignore")


class Sheintuch_matching:
  
  def __init__(self,basePath,mouse,sessions,footprints_file='results_redetect.mat',SNR_thr=2.,r_thr=0.0,d_thr=12,nbins=50,qtl=[0.05,0.95],use_kde=True,model='new'):
    ## build histograms:
    # nearest neighbour / others
    # distance & footprint correlation 
    ### loading data
    
    print('shifted vs unshifted version')
    print('matching assessment (how to? from Sheintuch)')
    self.para = {'pathMouse': pathcat([basePath,mouse]),
                 'sessions':  list(chain.from_iterable(sessions)) if (type(sessions[0]) is range) else list(chain(sessions)),
                 'fp_file':   footprints_file,
                 'nbins':     nbins,
                 'd_thr':     d_thr,
                 'qtl':       qtl,
                 'dims':      (512,512),
                 'pxtomu':    530.68/512,
                 'SNR_thr':   SNR_thr,
                 'r_thr':     r_thr,
                 'use_kde':   use_kde,
                 'model':     model         ## can be 'new', 'old', or 'both'
                 }  
    
    self.nS = len(list(self.para['sessions']))
    
    self.data = {'nA':                np.zeros(self.nS,'int'),
                 'cm':                {},
                 'p_same':            {}}
    
    self.session_data = {'D_ROIs':[],
                         'fp_corr':[],
                         'fp_corr_max':[],
                         'nearest_neighbour':[]}
    
    self.model = {'counts':         np.zeros((self.para['nbins'],self.para['nbins'],3)),
                  'counts_old':     np.zeros((self.para['nbins'],self.para['nbins'],3)),
                  'fit_function':   {'distance':       {'NN':   [],
                                                        'nNN':  [],
                                                        'all':  []},
                                     'fp_correlation': {'NN':   [],
                                                        'nNN':  [],
                                                        'all':  []}},
                  'fit_parameter':  {'single':  {'distance':        {'NN':    [],
                                                                     'nNN':   [],
                                                                     'all':   []},
                                                 'fp_correlation':  {'NN':    [],
                                                                     'nNN':   [],
                                                                     'all':   []}},
                                     'joint': {}},
                  'pdf':            {'single':  {'distance':[],
                                                 'fp_correlation':[]},
                                     'joint':   []},
                  'p_same':         {'single':{},
                                     'joint':[]},
                  'kernel':         {'idxes':   {},
                                     'kde':     {}}
                  }
  
  def run_matching(self):
    print('Now analyzing mouse %s'%self.para['pathMouse'])
    print('Building model for matching ...')
    self.run_analysis(save_results=True)
    print('Matching neurons ...')
    self.run_registration(save_results=True)
    print('Done!')
    
  def run_analysis(self,save_results=False,suffix=''):
    self.progress = tqdm(zip(self.para['sessions'],range(self.nS)),total=self.nS)
    for (s0,s) in self.progress:
      
      self.A2 = self.load_footprints(s0)
      
      if s>0:
        self.progress.set_description('Aligning data from Session %d'%s0)
        self.prepare_footprints()
      self.A_idx = np.squeeze(np.array((self.A2>0).sum(0)>50,'bool'))
      self.data['nA'][s] = self.A2.shape[1]
        
      self.progress.set_description('Calculate neuron positions for Session %d'%s0)
      self.data['cm'][s] = com(self.A2,self.para['dims'][0],self.para['dims'][1]) * self.para['pxtomu']
      self.data['cm'][s][~self.A_idx,:] = np.NaN    ## some empty footprints might be in here
      
      if s>0:
        self.progress.set_description('Calculate statistics for Session %d'%s0)
        self.calculate_statistics(self.A_ref,self.A2,s)       # calculating distances and footprint correlations
      
      if self.para['use_kde']:
        self.progress.set_description('Calculate kernel density for Session %d'%s0)
        self.position_kde(self.A2,s,self.para['qtl'])         # build kernel
      
      if s>0:
        self.progress.set_description('Update model with data from Session %d'%s0)
        self.update_joint_model(s,self.para['use_kde'])
        
      self.A_ref = self.A2.copy()
    
    self.fit_model()
    
    if save_results:
      self.save_model(suffix=suffix)
  
  
  def run_registration(self, p_thr=0.05, plot_results=False, save_results=False, save_suffix=''):
    
    self.progress = tqdm(zip(self.para['sessions'][1:],range(1,self.nS)),total=self.nS,leave=False)
    self.A0 = self.load_footprints(self.para['sessions'][0])
    self.A_ref = self.A0.copy()
    self.data['nA'][0] = self.A_ref.shape[1]
    self.data['cm'][0] = com(self.A0,self.para['dims'][0],self.para['dims'][1]) * self.para['pxtomu']
    self.matches = np.zeros((self.data['nA'][0],self.nS))*np.NaN
    self.matches[:,0] = range(self.data['nA'][0])
    
    self.p_matched = np.zeros((self.data['nA'][0],self.nS))
    self.p_matched[:,0] = np.NaN
    
    for (s0,s) in self.progress:
      
      self.nA_ref = self.A_ref.shape[1]
      self.cm_ref = com(self.A_ref,self.para['dims'][0],self.para['dims'][1]) * self.para['pxtomu']
      
      self.A2 = self.load_footprints(s0)
      #tqdm.write('Session %d data contains %d neurons'%(s0))
      self.progress.set_description('A union size: %d, Aligning data from Session %d'%(self.nA_ref,s0))
      self.prepare_footprints(A_ref=self.A0)
      self.A_idx = np.squeeze(np.array((self.A2>0).sum(0)>50,'bool'))
      
      self.data['nA'][s] = self.A2.shape[1]
      self.data['cm'][s] = com(self.A2,self.para['dims'][0],self.para['dims'][1]) * self.para['pxtomu']
      self.data['cm'][s][~self.A_idx,:] = np.NaN    ## some empty footprints might be in here
      
      self.progress.set_description('A union size: %d, Calculate statistics for Session %d'%(self.nA_ref,s0))
      self.calculate_statistics(self.A_ref,self.A2,s,self.cm_ref)
      
      self.progress.set_description('A union size: %d, Obtaining matching probability for Session %d'%(self.nA_ref,s0))
      self.data['p_same'][s] = self.calculate_p(model='joint')
      
      #run hungarian matching with these (1-p_same) as score
      self.progress.set_description('A union size: %d, Perform Hungarian matching on Session %d'%(self.nA_ref,s0))
      matches, p_matched = self.find_matches(s, p_thr=p_thr, plot_results=plot_results)
      
      idx_TP = np.where(p_matched > p_thr)[0] ## thresholding results
      if len(idx_TP) > 0:
          matched_ref = matches[0][idx_TP]    # ground truth
          matched2 = matches[1][idx_TP]   # algorithm - comp
          non_matched_ref = np.setdiff1d(list(range(self.nA_ref)), matches[0][idx_TP])
          non_matched2 = np.setdiff1d(list(range(self.data['nA'][s])), matches[1][idx_TP])
          non_matched2 = non_matched2[self.A_idx[non_matched2]]
          TP = np.sum(p_matched > p_thr).astype('float32')
      
      self.A_ref = self.A_ref.tolil()
      #w = 1-p_matched[idx_TP]
      #print(p_matched[idx_TP])
      #print(p_matched[idx_TP].shape)
      #print(type(p_matched[idx_TP]))
      
      #self.A_ref[:,mat_un] = w*self.A_ref[:,matched_ref] + (1-w)*self.A2[:,matched2]
      self.A_ref[:,matched_ref] = self.A_ref[:,matched_ref].multiply(1-p_matched[idx_TP]) + self.A2[:,matched2].multiply(p_matched[idx_TP])
      self.A_ref = sp.sparse.hstack([self.A_ref, self.A2[:,non_matched2]]).asformat('csc')
      
      N_add = len(non_matched2)
      
      self.matches[matched_ref,s] = matched2
      match_add = np.zeros((N_add,self.nS))*np.NaN
      match_add[:,s] = non_matched2
      self.matches = np.vstack([self.matches,match_add])
      
      self.p_matched[matched_ref,s] = p_matched[idx_TP]
      p_same_add = np.zeros((N_add,self.nS))*np.NaN
      self.p_matched = np.vstack([self.p_matched,p_same_add])
      
    if save_results:
      self.save_registration(suffix=save_suffix)
    
  
  
  def calculate_p(self,model='joint'):
    """
    Returns:
        p_same: score matrix
    
    Raises:
        Exception: 'Nan value produced. Error in inputs'
    """
    
    d_arr = np.linspace(0,self.para['d_thr'],self.para['nbins']+1)[:-1]
    fp_arr = np.linspace(0,1,self.para['nbins']+1)[:-1]
    
    d_w = np.diff(d_arr)[0]
    fp_w = np.diff(fp_arr)[0]
    
    close_neighbours = self.session_data['D_ROIs'] < self.para['d_thr']
    D_idx = (self.session_data['D_ROIs'][close_neighbours] / d_w).astype('int')
    c_idx = (self.session_data['fp_corr'].toarray()[close_neighbours] / fp_w).astype('int')
    
    p_same = sp.sparse.lil_matrix(self.session_data['D_ROIs'].shape)
    p_same[close_neighbours] = self.model['p_same']['joint'][D_idx,c_idx]
    p_same.tocsc()
    return p_same
  
  
  def find_matches(self, s, p_thr=0.05,plot_results=False):
    
    matches = linear_sum_assignment(1 - self.data['p_same'][s].toarray())
    p_matched = self.data['p_same'][s].toarray()[matches]
    
    idx_TP = np.where(np.array(p_matched) > p_thr)[0] ## thresholding results
    if len(idx_TP) > 0:
        matched_ROIs1 = matches[0][idx_TP]    # ground truth
        matched_ROIs2 = matches[1][idx_TP]   # algorithm - comp
        non_matched1 = np.setdiff1d(list(range(self.nA_ref)), matches[0][idx_TP])
        non_matched2 = np.setdiff1d(list(range(self.data['nA'][s])), matches[1][idx_TP])
        TP = np.sum(np.array(p_matched) > p_thr).astype('float32')
    else:
        TP = 0.
        plot_results = False
        matched_ROIs1 = []
        matched_ROIs2 = []
        non_matched1 = list(range(self.nA_ref))
        non_matched2 = list(range(self.data['nA'][s]))
        
    FN = self.nA_ref - TP
    FP = self.data['nA'][s] - TP
    TN = 0
    
    performance = dict()
    performance['recall'] = TP / (TP + FN)
    performance['precision'] = TP / (TP + FP)
    performance['accuracy'] = (TP + TN) / (TP + FP + FN + TN)
    performance['f1_score'] = 2 * TP / (2 * TP + FP + FN)
    
    #print(performance)
    #print('are all empty neurons removed?')
    #print([((self.A_ref!=0).sum(0)<50).sum(),((self.A2!=0).sum(0)<50).sum()])
    
    if plot_results:
      
      print('plotting...')
      t_start = time.time()
      cmap = 'viridis'
      
      Cn = self.A_ref.sum(1).reshape(512,512)
      
#        try : #Plotting function
      level = 0.1
      plt.figure(figsize=(15,12))
      plt.rcParams['pdf.fonttype'] = 42
      font = {'family': 'Myriad Pro',
              'weight': 'regular',
              'size': 10}
      plt.rc('font', **font)
      lp, hp = np.nanpercentile(Cn, [5, 95])
      
      ax_matches = plt.subplot(121)
      ax_nonmatches = plt.subplot(122)
      
      ax_matches.imshow(Cn, vmin=lp, vmax=hp, cmap=cmap)
      ax_nonmatches.imshow(Cn, vmin=lp, vmax=hp, cmap=cmap)
      
      A = np.reshape(self.A_ref.astype('float32').toarray(), self.para['dims'] + (-1,), order='F').transpose(2, 0, 1)
      [ax_matches.contour(a, levels=[level], colors='w', linewidths=1) for a in A[matched_ROIs1,...]]
      [ax_nonmatches.contour(a, levels=[level], colors='w', linewidths=1) for a in A[non_matched1,...]]
      
      print('first half done - %5.3f'%(t_end-t_start))
      A = None
      A = np.reshape(self.A2.astype('float32').toarray(), self.para['dims'] + (-1,), order='F').transpose(2, 0, 1)
      
      #plt.title('Matches')
      #plt.axis('off')
      
      [ax_matches.contour(a, levels=[level], colors='r', linewidths=1) for a in A[matched_ROIs2,...]]
      [ax_nonmatches.contour(a, levels=[level], colors='r', linewidths=1) for a in A[non_matched2,...]]
      A = None
      #plt.title('Mismatches')
      #plt.axis('off')
      
      plt.draw()
      #plt.pause(1)
      t_end = time.time()
      print('done. time taken: %5.3f'%(t_end-t_start))
      plt.show(block=False)
       #except Exception as e:
#            logging.warning("not able to plot precision recall usually because we are on travis")
#            logging.warning(e)
    
    return matches, p_matched
  
  
  def load_footprints(self,s):
    
    pathData = pathcat([self.para['pathMouse'],'Session%02d'%s,self.para['fp_file']])
    try:
      ld = loadmat(pathData,variable_names=['A','idx_evaluate','SNR','r_values'],squeeze_me=True)
      idxes = (ld['SNR']>self.para['SNR_thr']) & (ld['r_values']>self.para['r_thr'])
      return ld['A'][:,idxes]
    except:
      ld = loadmat(pathData,variable_names=['A'],squeeze_me=True)
      return ld['A']
  
  
  def prepare_footprints(self,A_ref=None,align_flag=True,use_opt_flow=True,max_thr=0.001):
    
    if A_ref is None:
      A_ref = self.A_ref
    
    if 'csc_matrix' not in str(type(A_ref)):
        A_ref = sp.sparse.csc_matrix(A_ref)
    if 'csc_matrix' not in str(type(self.A2)):
        self.A2 = sp.sparse.csc_matrix(self.A2)
    
    if align_flag:  # first align ROIs from session 2 to the template from session 1
      t_start = time.time()
      
      (x_shift,y_shift),flow,(x_grid,y_grid) = get_shift_and_flow(A_ref,self.A2,self.para['dims'],projection=1,plot_bool=False)
      
      if use_opt_flow:    ## for each pixel, find according position in other map
        x_remap = (x_grid - x_shift + flow[:,:,0])
        y_remap = (y_grid - y_shift + flow[:,:,1])
      else:
        x_remap = (x_grid - x_shift)
        y_remap = (y_grid - y_shift)
        
      self.A2 = sp.sparse.hstack([sp.sparse.csc_matrix(cv2.remap(img.reshape(self.para['dims']), x_remap,y_remap, cv2.INTER_CUBIC).reshape(-1,1)) for img in self.A2.toarray().T])
      
    self.A_ref = sp.sparse.vstack([a.multiply(a>max_thr*a.max())/a.max() if (a>0).sum()>50 else sp.sparse.csr_matrix(a.shape) for a in self.A_ref.T]).T
    self.A2 = sp.sparse.vstack([a.multiply(a>max_thr*a.max())/a.max() if (a>0).sum()>50 else sp.sparse.csr_matrix(a.shape) for a in self.A2.T]).T
    
    
  def calculate_statistics(self,A1,A2,s,cm_ref=None,binary='half'):
    
    if cm_ref is None:
      cm_ref = self.data['cm'][s-1]
    nA_ref = cm_ref.shape[0]
    
    self.session_data['D_ROIs'] = sp.spatial.distance.cdist(cm_ref,self.data['cm'][s])
    self.session_data['fp_corr'] = sp.sparse.lil_matrix((nA_ref,self.data['nA'][s]))
    self.session_data['fp_corr_max'] = sp.sparse.lil_matrix((nA_ref,self.data['nA'][s]))
    
    self.session_data['nearest_neighbour'] = np.zeros((nA_ref,self.data['nA'][s]),'bool')
    
    idx_good = ~np.isnan(cm_ref[:,0])
    idx_NN = np.nanargmin(self.session_data['D_ROIs'][idx_good,:],axis=1)
    self.session_data['nearest_neighbour'][idx_good,idx_NN] = True
    
    for i in tqdm(range(nA_ref),desc='calculating footprint correlation of %d neurons'%nA_ref,leave=False):
      for j in np.where(self.session_data['D_ROIs'][i,:]<self.para['d_thr'])[0]:
        if self.A_idx[j]:
          #try:
          if (self.para['model']=='old') | (self.para['model']=='both'):
            self.session_data['fp_corr'][i,j], shift = calculate_img_correlation(A1[:,i],A2[:,j],shift=False)
          
          self.session_data['fp_corr_max'][i,j], shift = calculate_img_correlation(A1[:,i],A2[:,j],crop=True,shift=True,binary=binary)
            
          #except:
            #tqdm.write('correlation calculation failed for neurons [%d,%d]'%(i,j),end='')
    self.session_data['fp_corr'].tocsc()
    self.session_data['fp_corr_max'].tocsc()
    
  
  def update_joint_model(self,s,use_kde):
    
    distance_arr = np.linspace(0,self.para['d_thr'],self.para['nbins']+1)
    fpcorr_arr = np.linspace(0,1,self.para['nbins']+1)
    
    if use_kde:
      idxes = self.model['kernel']['idxes'][s-1]
    else:
      idxes = np.ones(self.session_data['D_ROIs'].shape[0],'bool') 
    
    ROI_close = self.session_data['D_ROIs'][idxes,:] < self.para['d_thr']
    
    NN_idx = self.session_data['nearest_neighbour'][idxes,:][ROI_close]
    D_ROIs = self.session_data['D_ROIs'][idxes,:][ROI_close]
    fp_corr = self.session_data['fp_corr'][idxes,:][ROI_close].toarray()
    fp_corr_max = self.session_data['fp_corr_max'][idxes,:][ROI_close].toarray()
    
    for i in tqdm(range(self.para['nbins']),desc='updating joint model',leave=False):
      idx_dist = (D_ROIs >= distance_arr[i]) & (D_ROIs < distance_arr[i+1])
      
      for j in range(self.para['nbins']):
        idx_fp = (fp_corr > fpcorr_arr[j]) & (fpcorr_arr[j+1] > fp_corr)
        idx_vals = idx_dist & idx_fp
        self.model['counts_old'][i,j,0] += np.count_nonzero(idx_vals)
        self.model['counts_old'][i,j,1] += np.count_nonzero(idx_vals & NN_idx)
        self.model['counts_old'][i,j,2] += np.count_nonzero(idx_vals & ~NN_idx)
        
        idx_fp = (fp_corr_max > fpcorr_arr[j]) & (fpcorr_arr[j+1] > fp_corr_max)
        idx_vals = idx_dist & idx_fp
        self.model['counts'][i,j,0] += np.count_nonzero(idx_vals)
        self.model['counts'][i,j,1] += np.count_nonzero(idx_vals & NN_idx)
        self.model['counts'][i,j,2] += np.count_nonzero(idx_vals & ~NN_idx)
        
  def position_kde(self,A,s,qtl=[0.05,0.95],plot_bool=False):
    
    #print('calculating kernel density estimates for session %d'%s)
    x_grid, y_grid = np.meshgrid(np.linspace(0,self.para['dims'][0]*self.para['pxtomu'],self.para['dims'][0]), np.linspace(0,self.para['dims'][1]*self.para['pxtomu'],self.para['dims'][1]))
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    kde = sp.stats.gaussian_kde(self.data['cm'][s][self.A_idx,:].T)
    self.model['kernel']['kde'][s] = np.reshape(kde(positions),x_grid.shape)
    
    cm_px = (self.data['cm'][s][self.A_idx,:]/self.para['pxtomu']).astype('int')
    kde_at_cm = np.zeros(self.data['nA'][s])*np.NaN
    kde_at_cm[self.A_idx] = self.model['kernel']['kde'][s][cm_px[:,1],cm_px[:,0]]
    self.model['kernel']['idxes'][s] = (kde_at_cm > np.quantile(self.model['kernel']['kde'][s],qtl[0])) & (kde_at_cm < np.quantile(self.model['kernel']['kde'][s],qtl[1]))
      
    if plot_bool:
      plt.figure()
      h_kde = plt.imshow(self.model['kernel']['kde'][s],cmap=plt.cm.gist_earth_r,origin='lower',extent=[0,self.para['dims'][0]*self.para['pxtomu'],0,self.para['dims'][1]*self.para['pxtomu']])
      #if s>0:
        #col = self.session_data['D_ROIs'].min(1)
      #else:
        #col = 'w'
      plt.scatter(self.data['cm'][s][:,0],self.data['cm'][s][:,1],c='w',s=5+10*self.model['kernel']['idxes'][s],clim=[0,10],cmap='YlOrRd')
      plt.xlim([0,self.para['dims'][0]*self.para['pxtomu']])
      plt.ylim([0,self.para['dims'][1]*self.para['pxtomu']])
      
      cm_px = (self.data['cm'][s]/self.para['pxtomu']).astype('int')
      kde_at_cm = self.model['kernel']['kde'][s][cm_px[:,1],cm_px[:,0]]
      plt.colorbar(h_kde)
      plt.show(block=False) 
  
  
  def fit_model(self,count_thr=0,model='new'):
    
    if ((self.para['model'] == 'old') | (self.para['model']=='both')) & (model=='old'):
      key_counts = 'counts_old'
    else:
      key_counts = 'counts'
    
    nbins = self.para['nbins']
    
    self.set_functions(model)
    
    d_arr = np.linspace(0,self.para['d_thr'],nbins+1)[1:]
    d_w = np.diff(d_arr)[0]
    fp_arr = np.linspace(0,1,nbins+1)[:-1]
    fp_w = np.diff(fp_arr)[0]
    
    bounds_p = np.array([(0,1)]).T
    bounds_d_NN = np.array([(0,np.inf),(-np.inf,np.inf)]).T
    bounds_d_nNN = np.array([(0,np.inf),(0,np.inf),(0,self.para['d_thr']/2)]).T
    bounds_d = np.hstack([bounds_p,bounds_d_NN,bounds_d_nNN])
    
    bounds_corr_NN = np.array([(0,np.inf),(-np.inf,np.inf)]).T
    if model=='old':
      bounds_corr_nNN = np.array([(-np.inf,np.inf),(-np.inf,np.inf)]).T
    else:
      bounds_corr_nNN = np.array([(0,np.inf),(-np.inf,np.inf)]).T
    bounds_corr = np.hstack([bounds_p,bounds_corr_NN,bounds_corr_nNN])
    
    ## build single models
    ### distance
    distance_NN_dat = self.model[key_counts][...,1].sum(1)/self.model[key_counts][...,1].sum()/d_w
    distance_nNN_dat = self.model[key_counts][...,2].sum(1)/self.model[key_counts][...,2].sum()/d_w
    distance_joint_dat = self.model[key_counts][...,0].sum(1)/self.model[key_counts][...,0].sum()/d_w
    self.model['fit_parameter']['single']['distance']['NN'] = curve_fit(self.model['fit_function']['distance']['NN'],d_arr,distance_NN_dat,bounds=bounds_d_NN)[0]
    self.model['fit_parameter']['single']['distance']['nNN'] = curve_fit(self.model['fit_function']['distance']['nNN'],d_arr,distance_nNN_dat,bounds=bounds_d_nNN)[0]
    p0 = (self.model[key_counts][...,1].sum()/self.model[key_counts][...,0].sum(),)+tuple(self.model['fit_parameter']['single']['distance']['NN'])+tuple(self.model['fit_parameter']['single']['distance']['nNN'])
    self.model['fit_parameter']['single']['distance']['all'] = curve_fit(self.model['fit_function']['distance']['all'],d_arr,distance_joint_dat,bounds=bounds_d,p0=p0)[0]
    
    
    ### to fp-correlation: NN - reverse lognormal, nNN - reverse lognormal
    fp_correlation_NN_dat = self.model[key_counts][...,1].sum(0)/self.model[key_counts][...,1].sum()/fp_w
    fp_correlation_nNN_dat = self.model[key_counts][...,2].sum(0)/self.model[key_counts][...,2].sum()/fp_w
    fp_correlation_joint_dat = self.model[key_counts][...,0].sum(0)/self.model[key_counts][...,0].sum()/fp_w
    self.model['fit_parameter']['single']['fp_correlation']['NN'] = curve_fit(self.model['fit_function']['fp_correlation']['NN'],fp_arr,fp_correlation_NN_dat,bounds=bounds_corr_NN)[0]
    self.model['fit_parameter']['single']['fp_correlation']['nNN'] = curve_fit(self.model['fit_function']['fp_correlation']['nNN'],fp_arr,fp_correlation_nNN_dat,bounds=bounds_corr_nNN)[0]
    p0 = (self.model[key_counts][...,1].sum()/self.model[key_counts][...,0].sum(),)+tuple(self.model['fit_parameter']['single']['fp_correlation']['NN'])+tuple(self.model['fit_parameter']['single']['fp_correlation']['nNN'])
    
    self.model['fit_parameter']['single']['fp_correlation']['all'] = curve_fit(self.model['fit_function']['fp_correlation']['all'],fp_arr,fp_correlation_joint_dat,bounds=bounds_corr,p0=p0)[0]
    
    d_NN = fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['single']['distance']['all'][1:3])*self.model['fit_parameter']['single']['distance']['all'][0]
    d_total = fun_wrapper(self.model['fit_function']['distance']['all'],d_arr,self.model['fit_parameter']['single']['distance']['all'])
    self.model['p_same']['single']['distance'] = d_NN/d_total
    
    corr_NN = fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['all'][1:3])*self.model['fit_parameter']['single']['fp_correlation']['all'][0]
    corr_total = fun_wrapper(self.model['fit_function']['fp_correlation']['all'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['all'])
    self.model['p_same']['single']['fp_correlation'] = corr_NN/corr_total
    
    ## build joint model
    # preallocate
    self.model['fit_parameter']['joint'] = {
      'distance':{'NN':np.zeros((self.para['nbins'],len(self.model['fit_parameter']['single']['distance']['NN'])))*np.NaN,
                  'nNN':np.zeros((self.para['nbins'],2))*np.NaN},
                  #'all':np.zeros((self.para['nbins'],len(self.model['fit_parameter']['single']['distance']['all'])))*np.NaN},
      'fp_correlation':{'NN':np.zeros((self.para['nbins'],len(self.model['fit_parameter']['single']['fp_correlation']['NN'])))*np.NaN,
                        'nNN':np.zeros((self.para['nbins'],2))*np.NaN}
                        #'all':np.zeros((self.para['nbins'],len(self.model['fit_parameter']['single']['fp_correlation']['all'])))*np.NaN}
      }
    
    joint_hist_norm_dist = self.model[key_counts]/self.model[key_counts].sum(0)
    joint_hist_norm_dist[np.isnan(joint_hist_norm_dist)] = 0
    joint_hist_norm_corr = self.model[key_counts]/self.model[key_counts].sum(1)[:,np.newaxis,:]
    joint_hist_norm_corr[np.isnan(joint_hist_norm_corr)] = 0
    
    bounds_d_nNN = np.array([(0,np.inf),(-np.inf,np.inf)]).T
    bounds_corr_nNN = np.array([(0,np.inf),(-np.inf,np.inf)]).T
    counts_thr = 100
    for i in tqdm(range(nbins)):
      ### to distance distribution: NN - lognormal, nNN - large lognormal?!
      if self.model[key_counts][:,i,1].sum() > counts_thr:
        dat = joint_hist_norm_dist[:,i,1]/joint_hist_norm_dist[:,i,1].sum()*nbins/self.para['d_thr']
        self.model['fit_parameter']['joint']['distance']['NN'][i,:] = curve_fit(self.model['fit_function']['distance']['NN'],d_arr,dat,bounds=bounds_d_NN)[0]
      #else:
        #self.model['fit_parameter']['joint']['distance']['NN'][i,:] = 0
      if self.model[key_counts][:,i,2].sum() > counts_thr:
        dat = joint_hist_norm_dist[:,i,2]/joint_hist_norm_dist[:,i,2].sum()*nbins/self.para['d_thr']
        self.model['fit_parameter']['joint']['distance']['nNN'][i,:] = curve_fit(self.model['fit_function']['distance']['nNN_joint'],d_arr,dat,bounds=bounds_d_nNN)[0]
      #else:
        #self.model['fit_parameter']['joint']['distance']['nNN'][i,:] = 0
        
      #try:
        #p0 = np.hstack([self.model['p_same']['single']['fp_correlation'][i],self.model['fit_parameter']['joint']['distance']['NN'][i,:],self.model['fit_parameter']['joint']['distance']['nNN'][i,:]])
        #p0[np.isnan(p0)] = 0
        #dat = joint_hist_norm_dist[:,i,0]/joint_hist_norm_dist[:,i,0].sum()*nbins/self.para['d_thr']
        #self.model['fit_parameter']['joint']['distance']['all'][i,:] = curve_fit(self.model['fit_function']['distance']['all'],d_arr,dat,bounds=bounds_d,p0=p0)[0]
      #except:
        #1
        
      ### to fp-correlation: NN - reverse lognormal, nNN - reverse lognormal
      if self.model[key_counts][i,:,1].sum() > counts_thr:
        dat = joint_hist_norm_corr[i,:,1]/joint_hist_norm_corr[i,:,1].sum()*nbins
        self.model['fit_parameter']['joint']['fp_correlation']['NN'][i,:] = curve_fit(self.model['fit_function']['fp_correlation']['NN'],fp_arr,dat,bounds=bounds_corr_NN)[0]
      #else:
        #self.model['fit_parameter']['joint']['fp_correlation']['NN'][i,:] = 0
        
      if self.model[key_counts][i,:,2].sum() > counts_thr:
        dat = joint_hist_norm_corr[i,:,2]/joint_hist_norm_corr[i,:,2].sum()*nbins
        self.model['fit_parameter']['joint']['fp_correlation']['nNN'][i,:] = curve_fit(self.model['fit_function']['fp_correlation']['nNN_joint'],fp_arr,dat,bounds=bounds_corr_nNN)[0]
      #else:
        #self.model['fit_parameter']['joint']['fp_correlation']['nNN'][i,:] = 0
        
      #try:
        #p0 = np.hstack([self.model['p_same']['single']['distance'][i],self.model['fit_parameter']['joint']['fp_correlation']['NN'][i,:],self.model['fit_parameter']['joint']['fp_correlation']['nNN'][i,:]])
        ##p0[np.isnan(p0)] = 0
        #dat = joint_hist_norm_corr[i,:,1]/joint_hist_norm_corr[i,:,0].sum()*nbins
        #self.model['fit_parameter']['joint']['fp_correlation']['all'][i,:] = curve_fit(self.model['fit_function']['fp_correlation']['all'],fp_arr,dat,bounds=bounds_corr,p0=p0,method='dogbox')[0]
      #except:
        #1
      
    
    arrays = {'distance':         d_arr,
              'fp_correlation':   fp_arr}
    
    ## smooth parameter functions
    for key in ['distance','fp_correlation']:
      for pop in ['NN','nNN']:#,'all']
        #for ax in range(self.model['fit_parameter']['joint'][key][pop].shape(1)):
        self.model['fit_parameter']['joint'][key][pop] = sp.ndimage.median_filter(self.model['fit_parameter']['joint'][key][pop],[5,1])
        self.model['fit_parameter']['joint'][key][pop] = sp.ndimage.gaussian_filter(self.model['fit_parameter']['joint'][key][pop],[1,0])
        
        for ax in range(self.model['fit_parameter']['joint'][key][pop].shape[1]):
          ## find first/last index, at which parameter has a non-nan value
          nan_idx = np.isnan(self.model['fit_parameter']['joint'][key][pop][:,ax])
          if nan_idx[0]:  ## interpolate beginning
            idx = np.where(~nan_idx)[0][:20]
            y_arr = self.model['fit_parameter']['joint'][key][pop][idx,ax]
            f_interp = np.polyfit(arrays[key][idx],y_arr,1)
            fit_fun = np.poly1d(f_interp)
            
            self.model['fit_parameter']['joint'][key][pop][:idx[0],ax] = fit_fun(arrays[key][:idx[0]])
          if nan_idx[-1]:  ## interpolate beginning
            idx = np.where(~nan_idx)[0][-20:]
            y_arr = self.model['fit_parameter']['joint'][key][pop][idx,ax]
            f_interp = np.polyfit(arrays[key][idx],y_arr,1)
            fit_fun = np.poly1d(f_interp)
            
            self.model['fit_parameter']['joint'][key][pop][idx[-1]+1:,ax] = fit_fun(arrays[key][idx[-1]+1:])
    
    
    #np.nanmean(self.model['fit_parameter']['joint']['fp_correlation']['nNN'])
    #idxes = np.any(np.isnan(self.model['fit_parameter']['joint']['fp_correlation']['nNN']),1)
    #self.model['fit_parameter']['joint']['fp_correlation']['nNN'][idxes,:] = self.model['fit_parameter']['single']['fp_correlation']['nNN']
    #idxes = np.any(np.isnan(self.model['fit_parameter']['joint']['fp_correlation']['NN']),1)
    #self.model['fit_parameter']['joint']['fp_correlation']['NN'][idxes,:] = self.model['fit_parameter']['single']['fp_correlation']['nNN']
    
    ## define probability density functions
    self.model['pdf']['joint'] = np.zeros((2,nbins,nbins))
    for n in range(nbins):
      if not np.any(np.isnan(self.model['fit_parameter']['joint']['fp_correlation']['NN'][n,:])):
        f_NN = fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['joint']['fp_correlation']['NN'][n,:])
        self.model['pdf']['joint'][0,n,:] = f_NN*self.model['p_same']['single']['distance'][n]
    for n in range(nbins):
      if not np.any(np.isnan(self.model['fit_parameter']['joint']['fp_correlation']['nNN'][n,:])):
        self.model['pdf']['joint'][1,n,:] = fun_wrapper(self.model['fit_function']['fp_correlation']['nNN_joint'],fp_arr,self.model['fit_parameter']['joint']['fp_correlation']['nNN'][n,:])*(1-self.model['p_same']['single']['distance'][n])
    
    #self.model['pdf']['joint'] = np.zeros((2,nbins,nbins))
    #for n in range(nbins):
      #if not np.any(np.isnan(self.model['fit_parameter']['joint']['distance']['NN'][n,:])):
        #f_NN = fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['joint']['distance']['NN'][n,:])
        #self.model['pdf']['joint'][0,:,n] = f_NN*self.model['p_same']['single']['fp_correlation'][n]
    #for n in range(nbins):
      #if not np.any(np.isnan(self.model['fit_parameter']['joint']['distance']['nNN'][n,:])):
        #f_nNN = fun_wrapper(self.model['fit_function']['distance']['nNN_joint'],d_arr,self.model['fit_parameter']['joint']['distance']['nNN'][n,:])
        #self.model['pdf']['joint'][1,:,n] = f_nNN*(1-self.model['p_same']['single']['fp_correlation'][n])
    
    ## obtain probability of being same neuron
    self.model['p_same']['joint'] = self.model['pdf']['joint'][0,...]/np.nansum(self.model['pdf']['joint'],0)
    #if count_thr > 0:
      #self.model['p_same']['joint'] *= np.minimum(self.model[key_counts][...,0],count_thr)/count_thr
    #sp.ndimage.filters.gaussian_filter(self.model['p_same']['joint'],2,output=self.model['p_same']['joint'])
    
    
  
  def set_functions(self,model='new'):
    
    functions = {}
    ## define some possible fitting functions
    functions['lognorm'] = lambda x,sigma,mu : 1/(x*sigma*np.sqrt(2*np.pi))*np.exp(-(np.log(x)-mu)**2/(2*sigma**2))
    functions['lognorm_reverse'] = lambda x,sigma,mu : 1/((1-x)*sigma*np.sqrt(2*np.pi))*np.exp(-(np.log(1-x)-mu)**2/(2*sigma**2))
    functions['lognorm_shifted'] = lambda x,sigma,mu,s : 1/((x+s)*sigma*np.sqrt(2*np.pi))*np.exp(-(np.log(x+s)-mu)**2/(2*sigma**2))
    functions['lognorm_reverse_shifted'] = lambda x,sigma,mu,s : 1/((-x+s)*sigma*np.sqrt(2*np.pi))*np.exp(-(np.log(-x+s)-mu)**2/(2*sigma**2))
    functions['gauss'] = lambda x,sigma,mu : 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))
    
    phi = lambda x : 1/np.sqrt(2*np.pi)*np.exp(-1/2*x**2)
    psi = lambda x : 1/2*(1 + sp.special.erf(x/np.sqrt(2)))
    functions['truncated_lognorm'] = lambda x,sigma,mu : 1/sigma * phi((x-mu)/sigma) / (psi((1-mu)/sigma) - psi((0-mu)/sigma))
    functions['truncated_lognorm_reverse'] = lambda x,sigma,mu : 1/sigma * phi((1-x-mu)/sigma) / (psi((1-mu)/sigma) - psi((0-mu)/sigma))
    
    functions['beta'] = lambda x,a,b : x**(a-1)*(1-x)**(b-1) / (math.gamma(a)*math.gamma(b)/math.gamma(a+b))
    functions['linear_sigmoid'] = lambda x,m,sig_slope,sig_center : m*x/(1+np.exp(-sig_slope*(x-sig_center)));
    
    ## set functions for model
    self.model['fit_function']['distance']['NN'] = functions['lognorm']
    self.model['fit_function']['distance']['nNN'] = functions['linear_sigmoid']
    self.model['fit_function']['distance']['nNN_joint'] = functions['gauss']
    self.model['fit_function']['distance']['all'] = lambda x,p,sigma,mu,m,sig_slope,sig_center : p*functions['lognorm'](x,sigma,mu) + (1-p)*functions['linear_sigmoid'](x,m,sig_slope,sig_center)
    
    self.model['fit_function']['fp_correlation']['NN'] = functions['lognorm_reverse']
    if model == 'new':
      self.model['fit_function']['fp_correlation']['nNN'] = functions['gauss']
      self.model['fit_function']['fp_correlation']['all'] = lambda x,p,sigma1,mu1,sigma2,mu2 : p*functions['lognorm_reverse'](x,sigma1,mu1) + (1-p)*functions['gauss'](x,sigma2,mu2)
    else:
      self.model['fit_function']['fp_correlation']['nNN'] = functions['beta']
      self.model['fit_function']['fp_correlation']['all'] = lambda x,p,sigma1,mu1,a,b : p*functions['lognorm_reverse'](x,sigma1,mu1) + (1-p)*functions['beta'](x,a,b)
    self.model['fit_function']['fp_correlation']['nNN_joint'] = functions['gauss']
    
    
  def plot_model(self,animate=False,model='new'):
    
    if ((self.para['model'] == 'old') | (self.para['model']=='both')) & (model=='old'):
      key_counts = 'counts_old'
    else:
      key_counts = 'counts'
    
    nbins = self.model[key_counts].shape[0]
    
    d_arr = np.linspace(0,self.para['d_thr'],nbins+1)[1:]
    d_w = np.diff(d_arr)[0]
    
    fp_arr = np.linspace(0,1,nbins+1)[:-1]
    fp_w = np.diff(fp_arr)[0]
    X, Y = np.meshgrid(fp_arr, d_arr)
    
    mean_corr_NN, var_corr_NN = mean_of_trunc_lognorm(self.model['fit_parameter']['joint']['fp_correlation']['NN'][:,1],self.model['fit_parameter']['joint']['fp_correlation']['NN'][:,0],[0,1])
    #mean_corr_NN, var_corr_NN = 
    
    if model == 'old':
      a = self.model['fit_parameter']['joint']['fp_correlation']['nNN'][:,0]
      b = self.model['fit_parameter']['joint']['fp_correlation']['nNN'][:,1]
      #mean_corr_nNN = a/(a+b)
      #var_corr_nNN = a*b/((a+b)**2*(a+b+1))
      mean_corr_nNN = b
      var_corr_nNN = a
    else:
      #mean_corr_nNN, var_corr_nNN = mean_of_trunc_lognorm(self.model['fit_parameter']['joint']['fp_correlation']['nNN'][:,1],self.model['fit_parameter']['joint']['fp_correlation']['nNN'][:,0],[0,1])
      mean_corr_nNN = self.model['fit_parameter']['joint']['fp_correlation']['nNN'][:,1]
      var_corr_nNN = self.model['fit_parameter']['joint']['fp_correlation']['nNN'][:,0]
    mean_corr_NN = 1-mean_corr_NN
    #mean_corr_nNN = 1-mean_corr_nNN
    
    
    if (model=='old') | (model=='new'):
      
      plt.figure(figsize=(8,4))
      ax_phase = plt.axes([0.3,0.13,0.2,0.4])
      #ax_phase.imshow(self.model[key_counts][:,:,0],extent=[0,1,0,self.para['d_thr']],aspect='auto',clim=[0,0.25*self.model[key_counts][:,:,0].max()],origin='lower')
      NN_ratio = self.model[key_counts][:,:,1]/self.model[key_counts][:,:,0]
      NN_ratio[self.model[key_counts][:,:,0]<50] = np.NaN
      ax_phase.imshow(NN_ratio,extent=[0,1,0,self.para['d_thr']],aspect='auto',clim=[0,0.5],origin='lower',cmap='RdYlGn')
      nlev = 3
      col = (np.ones((nlev,3)).T*np.linspace(0,1,nlev)).T
      p_levels = ax_phase.contour(X,Y,self.model['p_same']['joint'],levels=[0.05,0.5,0.95],colors=col)
      ax_phase.set_xlim([0,1])
      ax_phase.set_ylim([0,self.para['d_thr']])
      ax_phase.tick_params(axis='x',which='both',bottom=True,top=True,labelbottom=False,labeltop=True)
      ax_phase.tick_params(axis='y',which='both',left=True,right=True,labelright=False,labelleft=True)
      ax_phase.yaxis.set_label_position("right")
      #ax_phase.xaxis.tick_top()
      ax_phase.set_xlabel('correlation',fontsize=14)
      ax_phase.set_ylabel('distance',fontsize=14)
      
      
      ax_dist = plt.axes([0.05,0.13,0.175,0.4])
      ax_dist.barh(d_arr,self.model[key_counts][...,0].sum(1).flat,self.para['d_thr']/nbins,facecolor='k',alpha=0.5,orientation='horizontal')
      ax_dist.barh(d_arr,self.model[key_counts][...,2].sum(1),d_w,facecolor='salmon',alpha=0.5)
      ax_dist.barh(d_arr,self.model[key_counts][...,1].sum(1),d_w,facecolor='lightgreen',alpha=0.5)
      ax_dist.invert_xaxis()
      #h_d_move = ax_dist.bar(d_arr,np.zeros(nbins),d_w,facecolor='k')
      
      model_distance_all = (fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['single']['distance']['NN'])*self.model[key_counts][...,1].sum() + fun_wrapper(self.model['fit_function']['distance']['nNN'],d_arr,self.model['fit_parameter']['single']['distance']['nNN'])*self.model[key_counts][...,2].sum())*d_w
      
      ax_dist.plot(fun_wrapper(self.model['fit_function']['distance']['all'],d_arr,self.model['fit_parameter']['single']['distance']['all'])*self.model[key_counts][...,0].sum()*d_w,d_arr,'k')
      ax_dist.plot(model_distance_all,d_arr,'k--')
      
      ax_dist.plot(fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['single']['distance']['all'][1:3])*self.model['fit_parameter']['single']['distance']['all'][0]*self.model[key_counts][...,0].sum()*d_w,d_arr,'g')
      ax_dist.plot(fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['single']['distance']['NN'])*self.model[key_counts][...,1].sum()*d_w,d_arr,'g:')
      
      ax_dist.plot(fun_wrapper(self.model['fit_function']['distance']['nNN'],d_arr,self.model['fit_parameter']['single']['distance']['all'][3:])*(1-self.model['fit_parameter']['single']['distance']['all'][0])*self.model[key_counts][...,0].sum()*d_w,d_arr,'r')
      ax_dist.plot(fun_wrapper(self.model['fit_function']['distance']['nNN'],d_arr,self.model['fit_parameter']['single']['distance']['nNN'])*self.model[key_counts][...,2].sum()*d_w,d_arr,'r',linestyle=':')
      ax_dist.set_ylim([0,self.para['d_thr']])
      ax_dist.set_xlabel('counts',fontsize=14)
      ax_dist.spines['left'].set_visible(False)
      ax_dist.spines['top'].set_visible(False)
      ax_dist.tick_params(axis='y',which='both',left=False,right=True,labelright=False,labelleft=False)
      
      ax_corr = plt.axes([0.3,0.63,0.2,0.35])
      ax_corr.bar(fp_arr,self.model[key_counts][...,0].sum(0).flat,1/nbins,facecolor='k',alpha=0.5)
      ax_corr.bar(fp_arr,self.model[key_counts][...,2].sum(0),1/nbins,facecolor='salmon',alpha=0.5)
      ax_corr.bar(fp_arr,self.model[key_counts][...,1].sum(0),1/nbins,facecolor='lightgreen',alpha=0.5)
      
      f_NN = fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['NN'])
      f_nNN = fun_wrapper(self.model['fit_function']['fp_correlation']['nNN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['nNN'])
      model_fp_correlation_all = (f_NN*self.model[key_counts][...,1].sum() + f_nNN*self.model[key_counts][...,2].sum())*fp_w
      #ax_corr.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['all'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['all'])*self.model[key_counts][...,0].sum()*fp_w,'k')
      ax_corr.plot(fp_arr,model_fp_correlation_all,'k')
      
      #ax_corr.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['all'][1:3])*self.model['fit_parameter']['single']['fp_correlation']['all'][0]*self.model[key_counts][...,0].sum()*fp_w,'g')
      ax_corr.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['NN'])*self.model[key_counts][...,1].sum()*fp_w,'g')
      
      #ax_corr.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['nNN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['all'][3:])*(1-self.model['fit_parameter']['single']['fp_correlation']['all'][0])*self.model[key_counts][...,0].sum()*fp_w,'r')
      ax_corr.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['nNN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['nNN'])*self.model[key_counts][...,2].sum()*fp_w,'r')
      
      ax_corr.set_ylabel('counts',fontsize=14)
      ax_corr.set_xlim([0,1])
      ax_corr.spines['right'].set_visible(False)
      ax_corr.spines['top'].set_visible(False)
      ax_corr.tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=False,labeltop=False)
      
      #ax_parameter = 
      p_steps, rates = self.RoC(100)
      idx = np.where(p_steps == 0.5)[0]
    
      ax_cum = plt.axes([0.63,0.63,0.35,0.35])
      ax_cum.plot([0,1],[0.05,0.05],'grey',linestyle='--')
      ax_cum.plot([0.5,1],[0.95,0.95],'grey',linestyle='--')
      
      ax_cum.plot(rates['cumfrac']['joint'],p_steps[:-1],'grey',label='Joint')
      ax_cum.plot(rates['cumfrac']['distance'],p_steps[:-1],'k',label='Distance')
      ax_cum.plot(rates['cumfrac']['fp_correlation'],p_steps[:-1],'lightgrey',label='Correlation')
      ax_cum.set_ylabel('$p_{same}$',fontsize=14)
      ax_cum.set_xlabel('cumulative fraction',fontsize=14)
      #ax_cum.legend(fontsize=10,frameon=False)
      uncertain = {}
      idx_low = np.where(p_steps>0.05)[0][0]
      idx_high = np.where(p_steps<0.95)[0][-1]
      for key in rates['cumfrac'].keys():
        uncertain[key] = (rates['cumfrac'][key][idx_high] - rates['cumfrac'][key][idx_low])#/(1-rates['cumfrac'][key][idx_high+1])
      
      ax_uncertain = plt.axes([0.7,0.8,0.075,0.1])
      ax_uncertain.bar(2,uncertain['distance'],facecolor='k')
      ax_uncertain.bar(3,uncertain['joint'],facecolor='grey')
      if model=='old':
        ax_uncertain.bar(1,uncertain['fp_correlation'],facecolor='lightgrey')
        ax_uncertain.set_xlim([0.5,3.5])
        ax_uncertain.set_xticks(range(1,4))
        ax_uncertain.set_xticklabels(['Corr.','Dist.','Joint'],rotation=60,fontsize=10)
      else:
        ax_uncertain.set_xlim([1.5,3.5])
        ax_uncertain.set_xticks(range(2,4))
        ax_uncertain.set_xticklabels(['Dist.','Joint'],rotation=60,fontsize=10)
      ax_uncertain.set_ylim([0,0.1])
      ax_uncertain.spines['right'].set_visible(False)
      ax_uncertain.spines['top'].set_visible(False)
      ax_uncertain.set_title('uncertain fraction',fontsize=10)
      
      #ax_rates = plt.axes([0.83,0.6,0.15,0.3])
      #ax_rates.plot(rates['fp']['joint'],p_steps[:-1],'r',label='false positive rate')
      #ax_rates.plot(rates['tp']['joint'],p_steps[:-1],'g',label='true positive rate')
      
      #ax_rates.plot(rates['fp']['distance'],p_steps[:-1],'r--')
      #ax_rates.plot(rates['tp']['distance'],p_steps[:-1],'g--')
      
      #ax_rates.plot(rates['fp']['fp_correlation'],p_steps[:-1],'r:')
      #ax_rates.plot(rates['tp']['fp_correlation'],p_steps[:-1],'g:')
      #ax_rates.legend()
      #ax_rates.set_xlabel('rate')
      #ax_rates.set_ylabel('$p_{same}$')
      
      ax_RoC = plt.axes([0.63,0.13,0.35,0.3])
      ax_RoC.plot(rates['fp']['joint'],rates['tp']['joint'],'grey',label='Joint')
      ax_RoC.plot(rates['fp']['distance'],rates['tp']['distance'],'k',label='Distance')
      ax_RoC.plot(rates['fp']['fp_correlation'],rates['tp']['fp_correlation'],'lightgrey',label='Correlation')
      ax_RoC.plot(rates['fp']['joint'][idx],rates['tp']['joint'][idx],'kx')
      ax_RoC.plot(rates['fp']['distance'][idx],rates['tp']['distance'][idx],'kx')
      ax_RoC.plot(rates['fp']['fp_correlation'][idx],rates['tp']['fp_correlation'][idx],'kx')
      ax_RoC.set_ylabel('true positive',fontsize=14)
      ax_RoC.set_xlabel('false positive',fontsize=14)
      ax_RoC.spines['right'].set_visible(False)
      ax_RoC.set_xlim([0,0.5])
      ax_RoC.set_ylim([0,1])
      
      
      ax_fp = plt.axes([0.88,0.17,0.075,0.1])
      ax_fp.bar(2,rates['fp']['distance'][idx],facecolor='k')
      ax_fp.bar(3,rates['fp']['joint'][idx],facecolor='grey')
      ax_fp.set_xticks([])
      if model=='old':
        ax_fp.bar(1,rates['fp']['fp_correlation'][idx],facecolor='lightgrey')
        ax_fp.set_xlim([0.5,3.5])
      else:
        ax_fp.set_xlim([1.5,3.5])
        
      ax_fp.set_ylim([0,0.05])
      ax_fp.spines['right'].set_visible(False)
      ax_fp.spines['top'].set_visible(False)
      ax_fp.set_title('false positives',fontsize=10)
      
      ax_tp = plt.axes([0.72,0.17,0.075,0.1])
      ax_tp.bar(2,rates['tp']['distance'][idx],facecolor='k')
      ax_tp.bar(3,rates['tp']['joint'][idx],facecolor='grey')
      ax_tp.set_xticks([])
      if model=='old':
        ax_tp.bar(1,rates['tp']['fp_correlation'][idx],facecolor='lightgrey')
        ax_tp.set_xlim([0.5,3.5])
      else:
        ax_tp.set_xlim([1.5,3.5])
      ax_tp.set_ylim([0.8,1])
      ax_tp.spines['right'].set_visible(False)
      ax_tp.spines['top'].set_visible(False)
      ax_tp.set_ylabel('fraction',fontsize=10)
      ax_tp.set_title('true positives',fontsize=10)
      
      #plt.tight_layout()
      plt.show(block=False)
      ext = 'png'
      path = pathcat([self.para['pathMouse'],'Sheintuch_matching_%s.%s'%(model,ext)])
      plt.savefig(path,format=ext,dpi=300)
      
      #ax_cvc = plt.axes([0.65,0.1,0.2,0.4])
      #idx = self.session_data['fp_corr_max']>0
      #ax_cvc.scatter(self.session_data['fp_corr_max'][idx].toarray().flat,self.session_data['fp_corr'][idx].toarray().flat,c='k',marker='.')
      #ax_cvc.plot([0,1],[0,1],'r--')
      #ax_cvc.set_xlim([0,1])
      #ax_cvc.set_ylim([0,1])
      #ax_cvc.set_xlabel('shifted correlation')
      #ax_cvc.set_ylabel('unshifted correlation')
      
      #plt.show(block=False)
    #return
    ## plot fit parameters
    plt.figure()
    plt.subplot(221)
    plt.plot(fp_arr,self.model['fit_parameter']['joint']['distance']['NN'][:,1],'g',label='lognorm $\mu$')
    plt.plot(fp_arr,self.model['fit_parameter']['joint']['distance']['nNN'][:,1],'r',label='gauss $\mu$')
    
    #plt.plot(fp_arr,self.model['fit_parameter']['joint']['distance']['all'][:,2],'g--')
    #plt.plot(fp_arr,self.model['fit_parameter']['joint']['distance']['all'][:,5],'r--')
    plt.legend()
    
    plt.subplot(223)
    plt.plot(fp_arr,self.model['fit_parameter']['joint']['distance']['NN'][:,0],'g',label='lognorm $\sigma$')
    plt.plot(fp_arr,self.model['fit_parameter']['joint']['distance']['nNN'][:,0],'r',label='gauss $\sigma$')
    #plt.plot(fp_arr,self.model['fit_parameter']['joint']['distance']['nNN'][:,1],'r--',label='dist $\gamma$')
    plt.legend()
    
    plt.subplot(222)
    plt.plot(d_arr,mean_corr_NN,'g',label='lognorm $\mu$')#self.model['fit_parameter']['joint']['fp_correlation']['NN'][:,1],'g')#
    plt.plot(d_arr,mean_corr_nNN,'r',label='lognorm $\mu$')
    plt.legend()
    
    plt.subplot(224)
    plt.plot(d_arr,var_corr_NN,'g',label='lognorm $\sigma$')
    plt.plot(d_arr,var_corr_nNN,'r',label='lognorm $\sigma$')
    plt.legend()
    plt.tight_layout()
    plt.show(block=False)
    
    #return
    fig = plt.figure()
    plt.subplot(322)
    plt.plot(fp_arr,self.model['fit_parameter']['joint']['distance']['NN'][:,1],'g')
    plt.plot(fp_arr,self.model['fit_parameter']['joint']['distance']['nNN'][:,1],'r')
    plt.xlim([0,1])
    plt.ylim([0,self.para['d_thr']])
    
    plt.subplot(321)
    plt.plot(d_arr,mean_corr_NN,'g')#self.model['fit_parameter']['joint']['fp_correlation']['NN'][:,1],'g')#
    plt.plot(d_arr,mean_corr_nNN,'r')#self.model['fit_parameter']['joint']['fp_correlation']['nNN'][:,1],'r')#
    plt.xlim([0,self.para['d_thr']])
    plt.ylim([0.5,1])
    
    
    plt.subplot(323)
    plt.bar(d_arr,self.model[key_counts][...,0].sum(1).flat,self.para['d_thr']/nbins,facecolor='k',alpha=0.5)
    plt.bar(d_arr,self.model[key_counts][...,2].sum(1),d_w,facecolor='r',alpha=0.5)
    plt.bar(d_arr,self.model[key_counts][...,1].sum(1),d_w,facecolor='g',alpha=0.5)
    h_d_move = plt.bar(d_arr,np.zeros(nbins),d_w,facecolor='k')
    
    model_distance_all = (fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['single']['distance']['NN'])*self.model[key_counts][...,1].sum() + fun_wrapper(self.model['fit_function']['distance']['nNN'],d_arr,self.model['fit_parameter']['single']['distance']['nNN'])*self.model[key_counts][...,2].sum())*d_w
    
    plt.plot(d_arr,fun_wrapper(self.model['fit_function']['distance']['all'],d_arr,self.model['fit_parameter']['single']['distance']['all'])*self.model[key_counts][...,0].sum()*d_w,'k')
    plt.plot(d_arr,model_distance_all,'k--')
    
    plt.plot(d_arr,fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['single']['distance']['all'][1:3])*self.model['fit_parameter']['single']['distance']['all'][0]*self.model[key_counts][...,0].sum()*d_w,'g')
    plt.plot(d_arr,fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['single']['distance']['NN'])*self.model[key_counts][...,1].sum()*d_w,'g--')
    
    plt.plot(d_arr,fun_wrapper(self.model['fit_function']['distance']['nNN'],d_arr,self.model['fit_parameter']['single']['distance']['all'][3:])*(1-self.model['fit_parameter']['single']['distance']['all'][0])*self.model[key_counts][...,0].sum()*d_w,'r')
    plt.plot(d_arr,fun_wrapper(self.model['fit_function']['distance']['nNN'],d_arr,self.model['fit_parameter']['single']['distance']['nNN'])*self.model[key_counts][...,2].sum()*d_w,'r--')
    plt.xlim([0,self.para['d_thr']])
    plt.xlabel('distance')
    
    plt.subplot(324)
    plt.bar(fp_arr,self.model[key_counts][...,0].sum(0).flat,1/nbins,facecolor='k',alpha=0.5)
    plt.bar(fp_arr,self.model[key_counts][...,2].sum(0),1/nbins,facecolor='r',alpha=0.5)
    plt.bar(fp_arr,self.model[key_counts][...,1].sum(0),1/nbins,facecolor='g',alpha=0.5)
    h_fp_move = plt.bar(fp_arr,np.zeros(nbins),1/nbins,facecolor='k')
    
    model_fp_correlation_all = (fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['NN'])*self.model[key_counts][...,1].sum() + fun_wrapper(self.model['fit_function']['fp_correlation']['nNN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['nNN'])*self.model[key_counts][...,2].sum())*fp_w
    plt.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['all'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['all'])*self.model[key_counts][...,0].sum()*fp_w,'k')
    plt.plot(fp_arr,model_fp_correlation_all,'k--')
    
    plt.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['all'][1:3])*self.model['fit_parameter']['single']['fp_correlation']['all'][0]*self.model[key_counts][...,0].sum()*fp_w,'g')
    plt.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['NN'])*self.model[key_counts][...,1].sum()*fp_w,'g--')
    
    plt.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['nNN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['all'][3:])*(1-self.model['fit_parameter']['single']['fp_correlation']['all'][0])*self.model[key_counts][...,0].sum()*fp_w,'r')
    plt.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['nNN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['nNN'])*self.model[key_counts][...,2].sum()*fp_w,'r--')
    
    plt.xlabel('correlation')
    plt.xlim([0,1])
    
    ax_d = plt.subplot(326)
    d_bar1 = ax_d.bar(d_arr,self.model[key_counts][:,0,0],self.para['d_thr']/nbins,facecolor='k',alpha=0.5)
    d_bar2 = ax_d.bar(d_arr,self.model[key_counts][:,0,2],self.para['d_thr']/nbins,facecolor='r',alpha=0.5)
    d_bar3 = ax_d.bar(d_arr,self.model[key_counts][:,0,1],self.para['d_thr']/nbins,facecolor='g',alpha=0.5)
    
    ### now, plot model stuff
    d_model_nNN, = ax_d.plot(d_arr,fun_wrapper(self.model['fit_function']['distance']['nNN_joint'],d_arr,self.model['fit_parameter']['joint']['distance']['nNN'][0,:]),'r')
    d_model_NN, = ax_d.plot(d_arr,fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['joint']['distance']['NN'][0,:]),'g')
    
    h_d = [d_bar1,d_bar3,d_bar2,h_d_move,d_model_NN,d_model_nNN]
    #h_d = d_bar1
    ax_d.set_xlabel('distance')
    ax_d.set_xlim([0,self.para['d_thr']])
    ax_d.set_ylim([0,self.model[key_counts][...,0].max()*1.1])
    
    
    ax_fp = plt.subplot(325)
    fp_bar1 = ax_fp.bar(fp_arr,self.model[key_counts][0,:,0],1/nbins,facecolor='k',alpha=0.5)
    fp_bar2 = ax_fp.bar(fp_arr,self.model[key_counts][0,:,2],1/nbins,facecolor='r',alpha=0.5)
    fp_bar3 = ax_fp.bar(fp_arr,self.model[key_counts][0,:,1],1/nbins,facecolor='g',alpha=0.5)
    
    ### now, plot model stuff
    fp_model_nNN, = ax_fp.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['nNN_joint'],fp_arr,self.model['fit_parameter']['joint']['fp_correlation']['nNN'][0,:]),'r')
    fp_model_NN, = ax_fp.plot(fp_arr,fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['joint']['fp_correlation']['NN'][0,:]),'g')
    
    
    h_fp = [fp_bar1,fp_bar3,fp_bar2,h_fp_move,fp_model_NN,fp_model_nNN]
    ax_fp.set_xlabel('corr')
    ax_fp.set_xlim([0,1])
    ax_fp.set_ylim([0,self.model[key_counts][...,0].max()*1.1])
    
  
    def update_distr(i,h_d,h_fp):
      
      n = i%self.model[key_counts].shape[0]
      for k in range(3):
        [h.set_height(dat) for h,dat in zip(h_d[k],self.model[key_counts][:,n,k])]
        [h.set_height(dat) for h,dat in zip(h_fp[k],self.model[key_counts][n,:,k])]
      
      d_move = np.zeros(self.model[key_counts].shape[0])
      fp_move = np.zeros(self.model[key_counts].shape[1])
      
      d_move[n] = self.model[key_counts][n,:,0].sum()
      fp_move[n] = self.model[key_counts][:,n,0].sum()
      
      [h.set_height(dat) for h,dat in zip(h_d[3],d_move)]
      [h.set_height(dat) for h,dat in zip(h_fp[3],fp_move)]
      
      self.model['p_same']['single']['distance'][n]*self.model[key_counts][n,:,0].sum()
      (1-self.model['p_same']['single']['distance'][n])*self.model[key_counts][n,:,0].sum()
      
      
      h_fp[4].set_ydata(fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['joint']['fp_correlation']['NN'][n,:])*self.model['p_same']['single']['distance'][n]*self.model[key_counts][n,:,0].sum()*fp_w)
      h_fp[5].set_ydata(fun_wrapper(self.model['fit_function']['fp_correlation']['nNN_joint'],fp_arr,self.model['fit_parameter']['joint']['fp_correlation']['nNN'][n,:])*(1-self.model['p_same']['single']['distance'][n])*self.model[key_counts][n,:,0].sum()*fp_w)
      
      self.model['p_same']['single']['fp_correlation'][n]*self.model[key_counts][:,n,0].sum()
      (1-self.model['p_same']['single']['fp_correlation'][n])*self.model[key_counts][:,n,0].sum()
      
      h_d[4].set_ydata(fun_wrapper(self.model['fit_function']['distance']['NN'],d_arr,self.model['fit_parameter']['joint']['distance']['NN'][n,:])*self.model['p_same']['single']['fp_correlation'][n]*self.model[key_counts][:,n,0].sum()*d_w)
      
      h_d[5].set_ydata(fun_wrapper(self.model['fit_function']['distance']['nNN_joint'],d_arr,self.model['fit_parameter']['joint']['distance']['nNN'][n,:])*(1-self.model['p_same']['single']['fp_correlation'][n])*self.model[key_counts][:,n,0].sum()*d_w)
      #print(tuple(h_d[0]))
      return tuple(h_d[0]) + tuple(h_d[1]) + tuple(h_d[2]) + tuple(h_d[3]) + (h_d[4],) + (h_d[5],) + tuple(h_fp[0]) + tuple(h_fp[1]) + tuple(h_fp[2]) + tuple(h_fp[3]) + (h_fp[4],) + (h_fp[5],)
    
    plt.tight_layout()
    if animate or False:
      
      #Writer = animation.writers['ffmpeg']
      #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=900)
      
      anim = animation.FuncAnimation(fig, update_distr, fargs=(h_d,h_fp),frames=nbins,interval=100, blit=True)
      svPath = pathcat([self.para['pathMouse'],'animation_single_models_%s.gif'%model])
      anim.save(svPath, writer='imagemagick',fps=15)#writer)
      print('animation saved at %s'%svPath)
      #anim
      plt.show()
    else:
      update_distr(20,h_d,h_fp)
      plt.show(block=False)
    #return
    
    
    plt.figure()
    plt.subplot(221)
    plt.imshow(self.model[key_counts][:,:,0],extent=[0,1,0,self.para['d_thr']],aspect='auto',clim=[0,self.model[key_counts][:,:,0].max()],origin='lower')
    nlev = 3
    col = (np.ones((nlev,3)).T*np.linspace(0,1,nlev)).T
    p_levels = plt.contour(X,Y,self.model['p_same']['joint'],levels=[0.05,0.5,0.95],colors=col)
    #plt.colorbar(p_levels)
    #plt.imshow(self.model[key_counts][...,0],extent=[0,1,0,self.para['d_thr']],aspect='auto',origin='lower')
    plt.subplot(222)
    plt.imshow(self.model[key_counts][...,0],extent=[0,1,0,self.para['d_thr']],aspect='auto',origin='lower')
    plt.subplot(223)
    plt.imshow(self.model[key_counts][...,1],extent=[0,1,0,self.para['d_thr']],aspect='auto',origin='lower')
    plt.subplot(224)
    plt.imshow(self.model[key_counts][...,2],extent=[0,1,0,self.para['d_thr']],aspect='auto',origin='lower')
    plt.show(block=False)
    
    
    
    
    
    plt.figure()
    W_NN = self.model[key_counts][...,1].sum() / self.model[key_counts][...,0].sum()
    print(W_NN)
    #W_NN = 0.5
    W_nNN = 1-W_NN
    plt.subplot(211)
    pdf_NN = fun_wrapper(self.model['fit_function']['fp_correlation']['NN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['NN'])*fp_w
    pdf_nNN = fun_wrapper(self.model['fit_function']['fp_correlation']['nNN'],fp_arr,self.model['fit_parameter']['single']['fp_correlation']['nNN'])*fp_w
    pdf_all = pdf_NN*W_NN+pdf_nNN*W_nNN
    
    plt.plot(fp_arr,pdf_NN*W_NN,'g')
    plt.plot(fp_arr,pdf_nNN*W_nNN,'r')
    plt.plot(fp_arr,pdf_all)
    
    plt.subplot(212)
    plt.plot(fp_arr,pdf_NN*W_NN/pdf_all,'k')
    plt.ylim([0,1])
    
    plt.show(block=False)
    
    
    
    X, Y = np.meshgrid(fp_arr, d_arr)
    
    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,8),subplot_kw={'projection':'3d'})
    #ax = plt.subplot(221,projection='3d')
    prob = ax1.plot_surface(X,Y,self.model['pdf']['joint'][0,...],cmap='jet')
    prob.set_clim(0,6)
    ax1.set_xlabel('corr')
    ax1.set_ylabel('d')
    ax1.set_zlabel('model')
    #ax1 = plt.subplot(222,projection='3d')
    prob = ax2.plot_surface(X,Y,self.model['pdf']['joint'][1,...],cmap='jet')
    prob.set_clim(0,6)
    ax2.set_xlabel('corr')
    ax2.set_ylabel('d')
    ax2.set_zlabel('model')
    
    prob = ax3.plot_surface(X,Y,self.model['p_same']['joint'],cmap='jet')
    prob.set_clim(0,1)
    ax3.set_zlim([0,1])
    ax3.set_xlabel('corr')
    ax3.set_ylabel('d')
    ax3.set_zlabel('model')
    
    #ax = plt.subplot(224,projection='3d')
    prob = ax4.plot_surface(X,Y,self.model[key_counts][...,0],cmap='jet')
    #prob = ax.bar3d(X.flatten(),Y.flatten(),np.zeros((nbins,nbins)).flatten(),np.ones((nbins,nbins)).flatten()*fp_w,np.ones((nbins,nbins)).flatten()*d_w,self.model[key_counts][...,0].flatten(),cmap='jet')
    #prob.set_clim(0,1)
    ax4.set_xlabel('corr')
    ax4.set_ylabel('d')
    ax4.set_zlabel('occurence')
    plt.tight_layout()
    axes = [ax1,ax2,ax3,ax4]
    
    def rotate_view(i,axes,fixed_angle=30):
      for ax in axes:
        ax.view_init(fixed_angle,(i*2)%360)
      #return ax
    
    if animate:
      anim = animation.FuncAnimation(fig, rotate_view, fargs=(axes,30), frames=180, interval=100, blit=False)
      svPath = pathcat([self.para['pathMouse'],'animation_p_same.gif'])
      anim.save(svPath, writer='imagemagick',fps=15)#writer)
      print('animation saved at %s'%svPath)
      #anim
      #plt.show()
    else:
      rotate_view(100,axes,fixed_angle=30)
      plt.show(block=False)
    #anim.save('animation.mp4', writer=writer)
    #print('animation saved')
    
    print('proper weighting of bin counts')
    print('smoothing by gaussian')
  
  def RoC(self,steps,model='new'):
    
    if ((self.para['model'] == 'old') | (self.para['model']=='both')) & (model=='old'):
      key_counts = 'counts_old'
    else:
      key_counts = 'counts'
    
    p_steps = np.linspace(0,1,steps+1)
    
    rates = {'tp':      {},
             'tn':      {},
             'fp':      {},
             'fn':      {},
             'cumfrac': {}}
    
    for key in rates.keys():
      rates[key] = {'joint':np.zeros(steps),
                    'distance':np.zeros(steps),
                    'fp_correlation':np.zeros(steps)}
    
    nTotal = self.model[key_counts][...,0].sum()
    for i in range(steps):
      p = p_steps[i]
      
      for key in ['joint','distance','fp_correlation']:
        
        if key == 'joint':
          idxes_negative = self.model['p_same']['joint'] < p
          idxes_positive = self.model['p_same']['joint'] >= p
          
          tp = self.model[key_counts][idxes_positive,1].sum()
          tn = self.model[key_counts][idxes_negative,2].sum()
          fp = self.model[key_counts][idxes_positive,2].sum()
          fn = self.model[key_counts][idxes_negative,1].sum()
          
          rates['cumfrac']['joint'][i] = self.model[key_counts][idxes_negative,0].sum()/nTotal
        elif key == 'distance':
          idxes_negative = self.model['p_same']['single']['distance'] < p
          idxes_positive = self.model['p_same']['single']['distance'] >= p
          
          tp = self.model[key_counts][idxes_positive,:,1].sum()
          tn = self.model[key_counts][idxes_negative,:,2].sum()
          fp = self.model[key_counts][idxes_positive,:,2].sum()
          fn = self.model[key_counts][idxes_negative,:,1].sum()
          
          rates['cumfrac']['distance'][i] = self.model[key_counts][idxes_negative,:,0].sum()/nTotal
        else:
          idxes_negative = self.model['p_same']['single']['fp_correlation'] < p
          idxes_positive = self.model['p_same']['single']['fp_correlation'] >= p
          
          tp = self.model[key_counts][:,idxes_positive,1].sum()
          tn = self.model[key_counts][:,idxes_negative,2].sum()
          fp = self.model[key_counts][:,idxes_positive,2].sum()
          fn = self.model[key_counts][:,idxes_negative,1].sum()
          
          rates['cumfrac']['fp_correlation'][i] = self.model[key_counts][:,idxes_negative,0].sum()/nTotal
          
        rates['tp'][key][i] = tp/(fn+tp)
        rates['tn'][key][i] = tn/(fp+tn)
        rates['fp'][key][i] = fp/(fp+tn)
        rates['fn'][key][i] = fn/(fn+tp)
    
    return p_steps, rates
    
  def save_registration(self,suffix=''):
    
    results = {'assignment':self.matches,
               'p_matched':self.p_matched,
               'p_same':self.data['p_same'],
               'nA':self.data['nA']}
    
    results['cm'] = np.zeros(results['assignment'].shape + (2,))
    
    for s in range(self.nS):
      idx_c = np.where(~np.isnan(results['assignment'][:,s]))
      idx_n = results['assignment'][idx_c,s].astype('int')
      if s>0:
        results['p_same'][s] = sp.sparse.csr_matrix(results['p_same'][s])
      #if s==self.nS-1:
      results['cm'][idx_c,s,:] = self.data['cm'][s][idx_n,:]
      #else:
        #results['cm'][idx_c,s,:] = self.data['cm'][s][idx_c,:]
      
    
    pathSv = pathcat([self.para['pathMouse'],'matching/Sheintuch_registration_%s%s.pkl'%(os.path.splitext(self.para['fp_file'])[0],suffix)])
    pickleData(results,pathSv,'save')
  
  def load_registration(self,suffix=''):
    
    pathLd = pathcat([self.para['pathMouse'],'matching/Sheintuch_registration_%s%s.pkl'%(os.path.splitext(self.para['fp_file'])[0],suffix)])
    return pickleData([],pathLd,'load')
    
  def save_model(self,suffix=''):
    
    pathSv = pathcat([self.para['pathMouse'],'matching/Sheintuch_model_%s%s.pkl'%(os.path.splitext(self.para['fp_file'])[0],suffix)])
    
    results = {}
    for key in ['p_same','fit_parameter','pdf','counts','counts_old']:
      results[key] = self.model[key]
    pickleData(results,pathSv,'save')
    
  
  def load_model(self,suffix=''):
    pathLd = pathcat([self.para['pathMouse'],'matching/Sheintuch_model_%s%s.pkl'%(os.path.splitext(self.para['fp_file'])[0],suffix)])
    results = pickleData([],pathLd,'load')
    for key in results.keys():
      self.model[key] = results[key]
    self.para['nbins'] = self.model['p_same']['joint'].shape[0]
    

def mean_of_trunc_lognorm(mu,sigma,trunc_loc):
  
  alpha = (trunc_loc[0]-mu)/sigma
  beta = (trunc_loc[1]-mu)/sigma
  
  phi = lambda x : 1/np.sqrt(2*np.pi)*np.exp(-1/2*x**2)
  psi = lambda x : 1/2*(1 + sp.special.erf(x/np.sqrt(2)))
  
  trunc_mean = mu + sigma * (phi(alpha) - phi(beta))/(psi(beta) - psi(alpha))
  trunc_var = np.sqrt(sigma**2 * (1 + (alpha*phi(alpha) - beta*phi(beta))/(psi(beta) - psi(alpha)) - ((phi(alpha) - phi(beta))/(psi(beta) - psi(alpha)))**2))
  
  return trunc_mean,trunc_var

def norm_nrg(a_):

  a = a_.copy()
  dims = a.shape
  a = a.reshape(-1, order='F')
  indx = np.argsort(a, axis=None)[::-1]
  cumEn = np.cumsum(a.flatten()[indx]**2)
  cumEn /= cumEn[-1]
  a = np.zeros(np.prod(dims))
  a[indx] = cumEn
  return a.reshape(dims, order='F')
