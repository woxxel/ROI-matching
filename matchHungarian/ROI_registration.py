''' this code is largely based on the "roi.py" method of the CaImAn package - many of the things therefore are NOT my work and are copy-pasted from the original code. I have changed the code, however, to an extent that makes it easier to open a new file for it.
      I have changed:
        - alignment of sessions and footprints now per default via projection image of footprints, not from extra provided template
        - maintaining sparse structures when possible, for better memory performance
        - calculation of scores
    
    last changed: Alexander Schmidt - 12.04.2020
'''

import cv2, sys, progressbar, os, scipy, time
from tqdm import *
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Programs/PC_analysis')
from utils import calculate_img_correlation, get_shift_and_flow, com, pathcat
from scipy.io import loadmat

class match_ROIs:

  def __init__(self,basePath,mouse):
  
    self.para['pathMouse'] = pathcat([basePath,mouse])
    self.para['pxtomu'] = 530.68/512
    
    
  def run_registration(s1,s2):
    
    #pathData1 = pathcat([self.para['pathMouse'],'Session%02d'%s1,'results_redetect.mat'])
    #ld = loadmat(pathData1,variable_names=['A'])
    #self.A_ref = ld['A']
    
    #pathData2 = pathcat([self.para['pathMouse'],'Session%02d'%s2,'results_redetect.mat'])
    #ld = loadmat(pathData2,variable_names=['A'])
    #self.A2 = ld['A']
    
    #self.shift_correction()
    
    p_same = calculate_p(self.A1, self.A2, d_thr, dims=dims, pxtomu=pxtomu, enclosed_thr=enclosed_thr,model='joint')
    
    self.p_same = hungarian_ROI_registration(A1,A2,(512,512))
  
  #def shift_correction(self):
      #t_start = time.time()
      
      #if 'csc_matrix' not in str(type(self.A_ref)):
          #self.A_ref = scipy.sparse.csc_matrix(self.A_ref)
      #if 'csc_matrix' not in str(type(self.A2)):
          #self.A2 = scipy.sparse.csc_matrix(self.A2)
      
      #if align_flag:  # first align ROIs from session 2 to the template from session 1
        
        #(x_shift,y_shift),flow,(x_grid,y_grid) = get_shift_and_flow(self.A_ref,self.A2,dims,projection=1,plot_bool=True)
        
        #if use_opt_flow:    ## for each pixel, find according position in other map
          #x_remap = (x_grid - x_shift + flow[:,:,0])
          #y_remap = (y_grid - y_shift + flow[:,:,1])
        #else:
          #x_remap = (x_grid - x_shift)
          #y_remap = (y_grid - y_shift)
          
        #self.A2 = scipy.sparse.hstack([scipy.sparse.csc_matrix(cv2.remap(img.reshape(dims), x_remap,y_remap, cv2.INTER_CUBIC).reshape(-1,1)) for img in self.A2.toarray().T])
        
      #self.A_ref = scipy.sparse.vstack([a.multiply(a>max_thr*a.max())/a.sum() if a.sum()>0 else a for a in self.A_ref.T]).T
      #self.A2 = scipy.sparse.vstack([a.multiply(a>max_thr*a.max())/a.sum() if a.sum()>0 else a for a in self.A2.T]).T 
      #print("Correcting shifts done , t = %5.3fs"%(time.time()-t_start))
      
    
    
def hungarian_ROI_registration(A1, A2, dims, template1=None, template2=None, 
                               align_flag=True, use_opt_flow = True, 
                               max_thr = 0, thresh_cost=.7, d_thr=10, 
                               enclosed_thr=None, print_assignment=False, 
                               plot_results=False, Cn=None, cmap='viridis'):
    """
    Register ROIs across different sessions using an intersection over union 
    metric and the Hungarian algorithm for optimal matching

    Args:
        A1: ndarray or csc_matrix  # pixels x # of components
            ROIs from session 1

        A2: ndarray or csc_matrix  # pixels x # of components
            ROIs from session 2

        dims: list or tuple
            dimensionality of the FOV

        template1: ndarray dims
            template from session 1 - if 'None', template is calculated from projection of A1

        template2: ndarray dims
            template from session 2 - if 'None', template is calculated from projection of A2

        align_flag: bool
            align the templates before matching
        
        use_opt_flow: bool
            use dense optical flow to align templates
        
        max_thr: scalar
            max threshold parameter before binarization    
        
        thresh_cost: scalar
            maximum distance considered

        d_thr: scalar
            max distance between centroids

        enclosed_thr: float
            if not None set distance to at most the specified value when ground 
            truth is a subset of inferred

        print_assignment: bool
            print pairs of matched ROIs

        plot_results: bool
            create a plot of matches and mismatches

        Cn: ndarray
            background image for plotting purposes

        cmap: string
            colormap for background image

    Returns:
        matched_ROIs1: list
            indeces of matched ROIs from session 1

        matched_ROIs2: list
            indeces of matched ROIs from session 2

        non_matched1: list
            indeces of non-matched ROIs from session 1

        non_matched2: list
            indeces of non-matched ROIs from session 2

        performance:  list
            (precision, recall, accuracy, f_1 score) with A1 taken as ground truth

        A2: csc_matrix  # pixels x # of components
            ROIs from session 2 aligned to session 1

    """
    
    
    
    matches, costs = find_matches(1-p_same, print_assignment=print_assignment)
    matches = matches[0]
    costs = costs[0]
    ##print("Assignment done , t = %5.3fs"%(time.time()-t_start))
    
    ##%% store indeces
    #idx_tp = np.where(np.array(costs) < thresh_cost)[0]
    #if len(idx_tp) > 0:
        #matched_ROIs1 = matches[0][idx_tp]    # ground truth
        #matched_ROIs2 = matches[1][idx_tp]   # algorithm - comp
        #non_matched1 = np.setdiff1d(
            #list(range(D[0].shape[0])), matches[0][idx_tp])
        #non_matched2 = np.setdiff1d(
            #list(range(D[0].shape[1])), matches[1][idx_tp])
        #TP = np.sum(np.array(costs) < thresh_cost) * 1.
        
        #shifts_matched = [shifts[id1,id2,:] for id1, id2 in zip(matched_ROIs1,matched_ROIs2)]
    #else:
        #TP = 0.
        #plot_results = False
        #matched_ROIs1 = []
        #matched_ROIs2 = []
        #non_matched1 = list(range(D[0].shape[0]))
        #non_matched2 = list(range(D[0].shape[1]))
        
        #shifts_matched = []

    ##%% compute precision and recall

    #FN = D[0].shape[0] - TP
    #FP = D[0].shape[1] - TP
    #TN = 0

    #performance = dict()
    #performance['recall'] = old_div(TP, (TP + FN))
    #performance['precision'] = old_div(TP, (TP + FP))
    #performance['accuracy'] = old_div((TP + TN), (TP + FP + FN + TN))
    #performance['f1_score'] = 2 * TP / (2 * TP + FP + FN)
    #logging.info(performance)
    
    #if plot_results:
        #if Cn is None:
            #if template1 is not None:
                #Cn = template1
            #elif template2 is not None:
                #Cn = template2
            #else:
                #Cn = np.reshape(A1.sum(1) + A2.sum(1), dims, order='F')

        #masks_1 = np.reshape(A1.toarray(), dims + (-1,),
                             #order='F').transpose(2, 0, 1)
        #masks_2 = np.reshape(A2.toarray(), dims + (-1,),
                             #order='F').transpose(2, 0, 1)
##        try : #Plotting function
        #level = 0.98
        #plt.figure(figsize=(15,12))
        #plt.rcParams['pdf.fonttype'] = 42
        #font = {'family': 'Myriad Pro',
                #'weight': 'regular',
                #'size': 10}
        #plt.rc('font', **font)
        #lp, hp = np.nanpercentile(Cn, [5, 95])
        #plt.subplot(1, 2, 1)
        #plt.imshow(Cn, vmin=lp, vmax=hp, cmap=cmap)
        #[plt.contour(norm_nrg(mm), levels=[level], colors='w', linewidths=1)
         #for mm in masks_1[matched_ROIs1]]
        #[plt.contour(norm_nrg(mm), levels=[level], colors='r', linewidths=1)
         #for mm in masks_2[matched_ROIs2]]
        #plt.title('Matches')
        #plt.axis('off')
        #plt.subplot(1, 2, 2)
        #plt.imshow(Cn, vmin=lp, vmax=hp, cmap=cmap)
        #[plt.contour(norm_nrg(mm), levels=[level], colors='w', linewidths=1)
         #for mm in masks_1[non_matched1]]
        #[plt.contour(norm_nrg(mm), levels=[level], colors='r', linewidths=1)
         #for mm in masks_2[non_matched2]]
        #plt.title('Mismatches')
        #plt.axis('off')
        #plt.draw()
        #plt.pause(1)
##        except Exception as e:
##            logging.warning("not able to plot precision recall usually because we are on travis")
##            logging.warning(e)

    #return matched_ROIs1, matched_ROIs2, non_matched1, non_matched2, performance, A2, scores, shifts_matched
  
  
def calculate_p(A1, A2, d_thr, dims=None, pxtomu=1, enclosed_thr=None,model='joint'):
    """
    Compute distance matrix based on an intersection over union metric. Matrix are compared in order,
    with matrix i compared with matrix i+1

    Args:
        A1,A2: 2D-array (px x N)
            The reference (A1) and new (A2) set of spatial footprints to compare
        
        d_thr: float
            maximum distance among centroids allowed between components. This corresponds to a distance
            at which two components are surely disjoined
        
        dims: tuple as (y,x)-value
            dimension of the imaging window
        
        enclosed_thr: float
            if not None set distance to at most the specified value when ground truth is a subset of inferred
    
    Returns:
        D: score matrix
    
    
    Raises:
        Exception: 'Nan value produced. Error in inputs'

    """
    
    # numbers of neurons in A1, A2
    n1 = A1.shape[1]
    n2 = A2.shape[1]
    
    # compute centroid positions
    cm1 = com(A1, dims[0], dims[1])*pxtomu
    cm2 = com(A2, dims[0], dims[1])*pxtomu
    
    D = scipy.spatial.distance.cdist(cm1,cm2)
    
    c = scipy.sparse.lil_matrix((n1,n2))
    
    #if enclosed_thr is not None:
      #gt_val = A1.T.dot(A1).diagonal()
    
    for i in tqdm(range(n1),desc='calculating footprint correlations'):
        for j in np.where(D[i,:]<d_thr)[0]:
            
            if A2[:,j].count_nonzero() > 10:
                c[i,j], shift_tmp = calculate_img_correlation(A1[:,i], A2[:,j], dims, crop=True, binary='half')
                if np.isnan(D[i, j]) | np.isnan(c[i, j]):
                    #D[i,j] = 1;
                    raise Exception('Nan value produced. Error in inputs')
        #else:
          #c[i, j] = np.NaN
    # D is (1 - spatial correlation of footprints) and normalized to values 0-1 to allow for Hungarian algorithm
    
    model = pickleData([],pathModel,'load')
    
    nbins = model['p_same']['joint'].shape[0]
    d_arr = np.linspace(0,d_thr,nbins+1)[:-1]
    fp_arr = np.linspace(0,1,nbins+1)[:-1]
    
    d_w = np.diff(d_arr)[0]
    fp_w = np.diff(fp_arr)[0]
    
    close_neighbours = D < d_thr
    
    D_idx = (D[close_neighbours] / d_w).astype('int')
    c_idx = (c[close_neighbours] / fp_w).astype('int')
    
    p_same = np.zeros((n1,n2))
    p_same[close_neighbours] = model['p_same']['joint'][D_idx,c_idx]
    
    return p_same