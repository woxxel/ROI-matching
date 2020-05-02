import numpy as np
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
import scipy as sp
import sys, os
from scipy.io import loadmat
from tqdm import *
from past.utils import old_div
from utils import pathcat, com, calculate_img_correlation, pickleData

sys.path.append('/home/wollex/Data/Science/PhD/Programs/CaImAn')
import caiman as cm

class CNMF_plots:
  
  def __init__(self):
    pathCNMF = '/home/wollex/Data/Science/PhD/Data/tmp/cnmf_results_762_10.pkl'
    ld = pickleData([],pathCNMF,'load')
    self.estimates = ld.estimates
    print(self.estimates)
  
    nCells,T = self.estimates.C.shape
  
    ### motion correction plot ?
    
    ### load memmap
    
    
    #pathMemmap = '/home/wollex/Data/Science/PhD/Data/tmp/memmap762_10__d1_512_d2_512_d3_1_order_C_frames_8989_.mmap'
    #Yr, dims, T = cm.load_memmap(pathMemmap)
    #images = np.reshape(Yr.T, [T] + list(dims), order='F')
    
    #Yr = np.transpose(np.reshape(images, (T, -1), order='F'))
    #print(Yr.shape)
    #print(Y.shape)
    
  def get_memmap(self,fname,sv_dir):
    
    #self.fname_memmap = pathcat([sv_dir,'memmap_%s_%d_'%('raw',1)])
    #if not os.path.exists(self.fname_memmap):
      #self.fname_memmap = cm.save_memmap([fname], base_name=self.fname_memmap, save_dir=sv_dir, n_chunks=200, order='C', dview=None)
      #print('writing done')
    
    Yr, dims, T = cm.load_memmap(self.fname_memmap)
    self.Y = np.reshape(Yr.T, (T,) + (512,512), order='F')
    
    #self.shift = np.zeros((T,2))
    #for t in tqdm(range(T-1)):
      #_, self.shift[t,:] = calculate_img_correlation(Y[t,100:400,100:400],Y[t+1,100:400,100:400],dims=(300,300),plot_bool=False)
    
    
    
  def plot(self,nCells=100):
    
    plt.rcParams['font.size'] = 12
    
    fr = 15
    plt.figure(figsize=(7,4))
    t_arr = np.linspace(0,8989/15,8989)
    ax_motion = plt.subplot(231)
    ax_motion.plot(t_arr,np.cumsum(self.shift[:,0],0))
    
    ax_motion2 = plt.subplot(234)
    ax_motion2.plot(t_arr,np.cumsum(self.shift[:,1],0))
    
    
    ax_sn = plt.subplot(132)
    ax_g = plt.subplot(133)
    np.random.seed()
    n_arr = [1840,1564,434,2477,2432,1441,2697,2426,800,1877]
    
    ax_sn.fill_between([0.25*fr,0.5*fr],[0,0],[0.00015],alpha=0.2)
    for n in n_arr:#np.random.randint(0,self.estimates.C.shape[0],10):
      fluor = self.estimates.C[n,:] + self.estimates.YrA[n,:]
      sn = GetSn(fluor, ax_sn, range_ff=[0.25,0.5], fr=fr, method='logmexp')
      g = estimate_time_constant(fluor, ax_g, p=1, sn=sn, lags=30, fudge_factor=1.)
    
    g = np.zeros(nCells)
    sn = np.zeros(nCells)
    for n in tqdm(range(nCells)):
      fluor = self.estimates.C[n,:] + self.estimates.YrA[n,:]
      fluor -= fluor.min()
      sn[n] = GetSn(fluor, None, range_ff=[0.25,0.5], fr=fr, method='logmexp')
      g[n] = estimate_time_constant(fluor, None, p=1, sn=sn[n], lags=5, fudge_factor=1.)
      sn[n] /= np.median(fluor)
    
    ax_sn.set_xlabel('f [Hz]')
    ax_sn.set_ylabel('PSD')
    ax_sn.set_yticks([])
    ax_sn.set_ylim([-0.00001,0.00015])
    
    ax_g.set_xlabel('$\Delta t$')
    ax_g.set_yscale('log')
    ax_g.set_ylabel('Acov')
    ax_g.set_yticks([])
    ax_g.set_ylim([10**(-6),2*10**(-4)])
    
    #ax = plt.subplot(234)
    #ax.scatter(sn,-1/(fr*np.log(g)),s=1,c='k')
    #ax.set_xlabel('noise')
    
    plt.show(block=False)

    #g, sn = estimate_parameters(fluor, p=1, sn=None, g=None, range_ff=[.25, .5],
                                      #method='logmexp', lags=5, fudge_factor=1.)
    ## noise estimation plot
    
    ## AR, timeconstant estimation plot
    
    ## ROI contour plot
    
    ## ROI statistics (sizes, overlaps)
    
  
  
  
class post_CNMF:
  def __init__(self,basePath,mouse,s,dataSet='OnACID'):
    pxtomu = 530.68/512
    
    pathMouse = pathcat([basePath,mouse])
    pathSession = pathcat([pathMouse,'Session%02d'%s])
    pathData = pathcat([pathSession,'results_OnACID.mat'])
    ld = loadmat(pathData,variable_names=['Cn'])
    
    pathData = pathcat([pathSession,'results_%s.mat'%dataSet])
    self.data = loadmat(pathData,variable_names=['A','C','SNR','r_values'])
    self.data['Cn'] = ld['Cn'].T
    
    self.cm = com(self.data['A'],512,512) * pxtomu
    self.D_ROIs = sp.spatial.distance.pdist(self.cm)
    self.D_ROIs_mat = sp.spatial.distance.squareform(self.D_ROIs)
    
    nCells = self.data['C'].shape[0]
    self.corrcoef = np.zeros((nCells,nCells))
    for n in tqdm(range(nCells)):
      for nn in range(n+1,nCells):
        self.corrcoef[n,nn] = np.corrcoef(self.data['C'][n,:],self.data['C'][nn,:])[0,1]
    
    self.Acorr = np.zeros((nCells,nCells))
    for n in tqdm(range(nCells)):
      idx_close=np.where(self.D_ROIs_mat[n,:]<20)[0]
      for nn in idx_close:
        self.Acorr[n,nn],_ = calculate_img_correlation(self.data['A'][:,n],self.data['A'][:,nn],shift=False)
    
    self.Asz = (self.data['A']>0).sum(0)
      
  def plot(self,nCells=None):
    
    nCells = self.data['C'].shape[0] if nCells is None else nCells
    plt.figure(figsize=(7,5))
    ### activity distribution
    ax_ROI = plt.axes([0.1,0.45,0.4,0.5])
    ax_ROI.imshow(self.data['Cn'],origin='lower')
    [ax_ROI.contour((a/a.max()).reshape(512,512).toarray(), levels=[0.2], colors='k', linewidths=0.2) for a in self.data['A'][:,:nCells].T]
    sbar = ScaleBar(530.68/512 *10**(-6),location='lower right')
    ax_ROI.add_artist(sbar)
    ax_ROI.set_xticks([])
    ax_ROI.set_yticks([])
    
    ### Calcium trace examples
    ax_Ca = plt.axes([0.55,0.65,0.4,0.3])
    t_arr = np.linspace(0,8989/15,8989)
    
    N = 6
    color = iter(plt.cm.rainbow(np.linspace(0,1,N)))
    for i,n in enumerate(np.random.randint(0,self.data['C'].shape[0],N)):
      col = next(color)
      c = self.data['C'][n,:]
      ax_Ca.plot(t_arr,c/c.max()+i,c=col,linewidth=0.5)
      a = self.data['A'][:,n].reshape(512,512).toarray()
      ax_ROI.contour(a/a.max(),levels=[0.2],colors=[col],linewidth=0.5)
    ax_Ca.set_xlabel('t [s]')
    ax_Ca.set_ylabel('Ca$^{2+}$')
    ax_Ca.set_yticks([])
    ax_Ca.spines['right'].set_visible(False)
    ax_Ca.spines['top'].set_visible(False)
    
    ax_Asz = plt.axes([0.1,0.1,0.15,0.15])
    ax_Asz.hist(self.Asz.flat,np.linspace(100,400,51),facecolor='k')
    ax_Asz.set_xlabel('size (px)')
    ax_Asz.set_ylabel('count')
    ax_Asz.set_yticks([])
    ax_Asz.spines['right'].set_visible(False)
    ax_Asz.spines['top'].set_visible(False)
    
    #### SNR, r value, CNN stuff
    idxes = np.triu_indices(self.D_ROIs_mat.shape[0])
    
    #### distance vs correlation plot
    
    ax_D = plt.axes([0.55,0.4,0.25,0.1])
    ax_D.hist(self.D_ROIs,np.linspace(0,700,101),facecolor='k')
    ax_D.set_ylabel('count')
    ax_D.set_yticks([])
    ax_D.set_xlabel('distance [$\mu$m]')
    ax_D.set_xlim([0,700])
    ax_D.spines['right'].set_visible(False)
    ax_D.spines['top'].set_visible(False)
    
    #plt.rcParams['scatter.marker'] = '.'
    ax_D_corr = plt.axes([0.55,0.1,0.25,0.25])
    ax_D_corr.plot(self.D_ROIs_mat[idxes],self.corrcoef[idxes],'.',color='k',markersize=1,markeredgewidth=0,markerfacecolor='k')
    ax_D_corr.set_xlabel('distance [$\mu$m]')
    ax_D_corr.set_ylabel('corr')
    ax_D_corr.set_xlim([0,700])
    ax_D_corr.set_ylim([-0.3,1])
    
    ax_corr = plt.axes([0.4,0.1,0.095,0.25])
    ax_corr.hist(self.corrcoef[idxes],np.linspace(-0.3,1,101),facecolor='k',orientation='horizontal')
    ax_corr.invert_xaxis()
    ax_corr.set_xticks([])
    ax_corr.set_yticks([])
    ax_corr.set_xlabel('count')
    ax_corr.set_ylim([-0.3,1])
    ax_corr.spines['left'].set_visible(False)
    ax_corr.spines['top'].set_visible(False)
    
    ax_Ac = plt.axes([0.85,0.4,0.075,0.1])
    ax_Ac.plot(self.Acorr[idxes],self.D_ROIs_mat[idxes],'.',color='r',markersize=1,markeredgewidth=0,markerfacecolor='r')
    ax_Ac.set_ylim([0,20])
    ax_Ac.yaxis.tick_right()
    ax_Ac.yaxis.set_label_position("right")
    ax_Ac.set_ylabel('d [$\mu$m]')
    
    ax_Acc = plt.axes([0.85,0.1,0.075,0.25])
    ax_Acc.plot(self.Acorr[idxes],self.corrcoef[idxes],'.',color='r',markersize=1,markeredgewidth=0,markerfacecolor='r')
    ax_Acc.set_ylim([-0.3,1])
    ax_Acc.yaxis.set_label_position("right")
    ax_Acc.set_ylabel('Ca$^{2+}$ corr')
    ax_Acc.set_xlabel('$c_{fp}$')
    
    plt.tight_layout()
    plt.show(block=False)
    
    ext = 'png'
    path = '/home/wollex/Data/Science/PhD/Thesis/pics/Methods/postCNMF.%s'%ext
    plt.savefig(path,format=ext,dpi=300)
  
  ### 
  
  
def estimate_time_constant(fluor, ax, p=2, sn=None, lags=5, fudge_factor=1.):
    """
    Estimate AR model parameters through the autocovariance function

    Args:
        fluor        : nparray
            One dimensional array containing the fluorescence intensities with
            one entry per time-bin.
    
        p            : positive integer
            order of AR system
    
        sn           : float
            noise standard deviation, estimated if not provided.
    
        lags         : positive integer
            number of additional lags where he autocovariance is computed
    
        fudge_factor : float (0< fudge_factor <= 1)
            shrinkage factor to reduce bias

    Returns:
        g       : estimated coefficients of the AR process
    """

    #if sn is None:
        #sn = GetSn(fluor)
    
    lags += p
    xc = axcov(fluor, lags)
    if not (ax is None):
      ax.plot(np.linspace(-2,2,2*lags+1),xc)
    xc = xc[:, np.newaxis]
    
    A = sp.linalg.toeplitz(xc[lags + np.arange(lags)],
                              xc[lags + np.arange(p)]) - sn**2 * np.eye(lags, p)
    g = np.linalg.lstsq(A, xc[lags + 1:], rcond=None)[0]
    gr = np.roots(np.concatenate([np.array([1]), -g.flatten()]))
    gr = old_div((gr + gr.conjugate()), 2.)
    np.random.seed(45) # We want some variability below, but it doesn't have to be random at
                       # runtime. A static seed captures our intent, while still not disrupting
                       # the desired identical results from runs.
    gr[gr > 1] = 0.95 + np.random.normal(0, 0.01, np.sum(gr > 1))
    gr[gr < 0] = 0.15 + np.random.normal(0, 0.01, np.sum(gr < 0))
    g = np.poly(fudge_factor * gr)
    g = -g[1:]

    return g.flatten()
  
  
def GetSn(fluor, ax, range_ff=[0.25, 0.5], fr=15, method='logmexp'):
    """    
    Estimate noise power through the power spectral density over the range of large frequencies    

    Args:
        fluor    : nparray
            One dimensional array containing the fluorescence intensities with
            one entry per time-bin.
    
        range_ff : (1,2) array, nonnegative, max value <= 0.5
            range of frequency (x Nyquist rate) over which the spectrum is averaged  
    
        method   : string
            method of averaging: Mean, median, exponentiated mean of logvalues (default)

    Returns:
        sn       : noise standard deviation
    """

    ff, Pxx = sp.signal.welch(fluor,fs=fr)
    
    ind1 = ff > range_ff[0]
    ind2 = ff < range_ff[1]
    ind = np.logical_and(ind1, ind2)
    Pxx_ind = Pxx[ind]
    sn = {
        'mean': lambda Pxx_ind: np.sqrt(np.mean(old_div(Pxx_ind, 2))),
        'median': lambda Pxx_ind: np.sqrt(np.median(old_div(Pxx_ind, 2))),
        'logmexp': lambda Pxx_ind: np.sqrt(np.exp(np.mean(np.log(old_div(Pxx_ind, 2)))))
    }[method](Pxx_ind)
    
    if not (ax is None):
      ax.plot(ff[:],Pxx[:])
    
    return sn
  
def axcov(data, maxlag=5):
    """
    Compute the autocovariance of data at lag = -maxlag:0:maxlag

    Args:
        data : array
            Array containing fluorescence data
    
        maxlag : int
            Number of lags to use in autocovariance calculation

    Returns:
        axcov : array
            Autocovariances computed from -maxlag:0:maxlag
    """

    data = data - np.mean(data)
    T = len(data)
    bins = np.size(data)
    xcov = np.fft.fft(data, np.power(2, nextpow2(2 * bins - 1)))
    xcov = np.fft.ifft(np.square(np.abs(xcov)))
    xcov = np.concatenate([xcov[np.arange(xcov.size - maxlag, xcov.size)],
                           xcov[np.arange(0, maxlag + 1)]])
    return np.real(old_div(xcov, T))
  

def nextpow2(value):
    """
    Find exponent such that 2^exponent is equal to or greater than abs(value).

    Args:
        value : int

    Returns:
        exponent : int
    """

    exponent = 0
    avalue = np.abs(value)
    while avalue > np.power(2, exponent):
        exponent += 1
    return exponent
