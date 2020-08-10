import numpy as np
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib import rc
import scipy as sp
from scipy.stats import norm

import sys, os, time
from scipy.io import loadmat
from tqdm import *
from past.utils import old_div
from utils import pathcat, com, calculate_img_correlation, pickleData

sys.path.append('/home/wollex/Data/Science/PhD/Programs/CaImAn')
import caiman as cm

class CNMF_plots:

  def __init__(self,pathCNMF='/home/wollex/Data/Science/PhD/Data/tmp/cnmf_results_762_10.pkl'):

    if not (pathCNMF==''):
      ld = pickleData([],pathCNMF,'load')
      self.estimates = ld.estimates
      print(self.estimates)

      nCells,T = self.estimates.C.shape

    ### motion correction plot ?

    ### load memmap


    pathMemmap = '/home/wollex/Data/Science/PhD/Data/tmp/memmap762_10__d1_512_d2_512_d3_1_order_C_frames_8989_.mmap'
    #Yr, dims, T = cm.load_memmap(pathMemmap)
    #images = np.reshape(Yr.T, [T] + list(dims), order='F')

    #Yr = np.transpose(np.reshape(images, (T, -1), order='F'))
    #print(Yr.shape)
    #print(Y.shape)

  def get_memmap(self,fname='/home/wollex/Data/Science/PhD/Data/tmp/memmap762_10__d1_512_d2_512_d3_1_order_C_frames_8989_.mmap',sv_dir=''):

    self.fname_memmap = fname
    #self.fname_memmap = pathcat([sv_dir,'memmap_%s_%d_'%('raw',1)])
    #if not os.path.exists(self.fname_memmap):
      #self.fname_memmap = cm.save_memmap([fname], base_name=self.fname_memmap, save_dir=sv_dir, n_chunks=200, order='C', dview=None)
      #print('writing done')

    self.Yr, dims, T = cm.load_memmap(self.fname_memmap)
    self.Y = np.reshape(self.Yr.T, (T,) + (512,512), order='F')

    #t_start = time.time()
    #r = self.xy_align(range(T))
    #self.shift = r['shift']
    #t_end = time.time()
    #print('alignment done - time taken: %5.3g'%(t_end-t_start))

  def run_movie(self,T=10,f=15):

    plt.figure(figsize=(20,10))
    ax_raw = plt.subplot(121)
    ax_aligned = plt.subplot(122)

    N = (512,512)
    im_raw = ax_raw.imshow(np.zeros(N),origin='lower')
    im_aligned = ax_aligned.imshow(np.zeros(N),origin='lower')

    [ax_aligned.contour((a/a.max()).reshape(512,512).toarray(), levels=[0.2], colors='w', linewidths=1) for a in self.Ain[:,:500].T]
    [ax_aligned.contour((a/a.max()).reshape(512,512).toarray(), levels=[0.2], colors='r', linewidths=1) for a in self.Aout[:,:500].T]

    plt.show(block=False)
    self.shift = self.shift.astype('int')
    im_raw.set_clim(0,1)
    im_aligned.set_clim(0,1)
    #def anim_movie(t,im_raw,im_aligned):

    for t in range(T):
      tmp = self.Y[t,...]/self.Y[t,...].max()
      tmp_aligned = np.zeros(N)*np.NaN
      tmp_aligned[max(0,self.shift[t,0]):N[0]+min(0,self.shift[t,0]),
            max(0,self.shift[t,1]):N[1]+min(0,self.shift[t,1])] = \
            tmp[-min(0,self.shift[t,0]):N[0]-max(0,self.shift[t,0]),
                         -min(0,self.shift[t,1]):N[1]-max(0,self.shift[t,1])]

      im_raw.set_data(tmp)
      im_aligned.set_data(tmp_aligned)


      plt.title(self.shift[t,:])
      plt.draw()
      plt.waitforbuttonpress(1/f)


  def plot_r_val(self,Y,A,C,b,f, Athresh=0.1, Npeaks=5, tB=-3, tA=10, thres=0.3):
    K, _ = np.shape(C)
    A = csc_matrix(A)
    AA = (A.T * A).toarray()
    nA = np.sqrt(np.array(A.power(2).sum(0)))
    AA = old_div(AA, np.outer(nA, nA.T))
    AA -= np.eye(K)

    LOC = find_activity_intervals(C, Npeaks=Npeaks, tB=tB, tA=tA, thres=thres)
    rval = np.zeros(K)

    significant_samples: List[Any] = []
    for i in range(K):
        if (i + 1) % 200 == 0:         # Show status periodically
            logging.info('Components evaluated:' + str(i))
        if LOC[i] is not None:
            atemp = A[:, i].toarray().flatten()
            atemp[np.isnan(atemp)] = np.nanmean(atemp)
            ovlp_cmp = np.where(AA[:, i] > Athresh)[0]
            indexes = set(LOC[i])
            for _, j in enumerate(ovlp_cmp):
                if LOC[j] is not None:
                    indexes = indexes - set(LOC[j])

            if len(indexes) == 0:
                indexes = set(LOC[i])
                logging.warning('Component {0} is only active '.format(i) +
                                'jointly with neighboring components. Space ' +
                                'correlation calculation might be unreliable.')

            indexes = np.array(list(indexes)).astype(np.int)
            px = np.where(atemp > 0)[0]
            if px.size < 3:
                logging.warning('Component {0} is almost empty. '.format(i) + 'Space correlation is set to 0.')
                rval[i] = 0
                significant_samples.append({0})
            else:
                ysqr = np.array(Y[px, :])
                ysqr[np.isnan(ysqr)] = np.nanmean(ysqr)
                mY = np.mean(ysqr[:, indexes], axis=-1)
                significant_samples.append(indexes)
                rval[i] = scipy.stats.pearsonr(mY, atemp[px])[0]

        else:
            rval[i] = 0
            significant_samples.append(0)

    return rval, significant_samples


  def xy_align(self,stacks_to_do):
    '''{
    stacks_to_do: range of frames

    Calculates rigid x-y shift for each frame by maximizing correlation
    between pairs of images recursively.
    r = running data structure with:
        .I mean image
        .n number of frames in mean
        .T x-y shifts for the frames
    #max_shift = search radius to maximize correlation
    Based on:
    https://scanbox.wordpress.com/2014/03/20/recursive-image-alignment-and-statistics/
    }'''
    nStacks = len(stacks_to_do)
    if nStacks==1:
      # if only one stack, assign initial values to parameters
      r = {'I':       np.squeeze(self.Y[stacks_to_do,...]),
           'n':       1,
           'shift':   np.zeros(2)}
    else:
      # split into two groups and run again recursively
      #print(stacks_to_do[:nStacks//2])
      #print(stacks_to_do[nStacks//2:])
      r_input = self.xy_align(stacks_to_do[:nStacks//2])
      r_goal = self.xy_align(stacks_to_do[nStacks//2:])

      r = {'I':[],
           'n':[],
           'shift':[]}
      # use the average of the selected channel to get alignment
      N = r_input['I'].shape
      y_bounds = np.where(~np.isnan(r_input['I'][:,N[1]//2] + r_goal['I'][:,N[1]//2]))[0][[0,-1]] + (0,1)
      x_bounds = np.where(~np.isnan(r_input['I'][N[0]//2,:] + r_goal['I'][N[0]//2,:]))[0][[0,-1]] + (0,1)

      x_range = np.abs(x_bounds - N[1]//2).min()
      x_idxes = [N[1]//2 - x_range,N[1]//2 + x_range]

      y_range = np.abs(y_bounds - N[0]//2).min()
      y_idxes = [N[0]//2 - y_range,N[0]//2 + y_range]

      _,shift = calculate_img_correlation(r_input['I'][y_idxes[0]:y_idxes[1],x_idxes[0]:x_idxes[1]],r_goal['I'][y_idxes[0]:y_idxes[1],x_idxes[0]:x_idxes[1]],dims=(y_idxes[1]-y_idxes[0],x_idxes[1]-x_idxes[0]),plot_bool=False)   ## shift in y,x order
      shift = -shift
      shift = shift.astype('int')
      # shift input image (fill empty spaces with NaN)
      if (np.abs(shift).sum())>0:
        tmp = np.zeros(N)*np.NaN
        tmp[max(0,shift[0]):N[0]+min(0,shift[0]),
            max(0,shift[1]):N[1]+min(0,shift[1])] = \
            r_input['I'][-min(0,shift[0]):N[0]-max(0,shift[0]),
                         -min(0,shift[1]):N[1]-max(0,shift[1])]
        r_input['I'] = tmp

      # combine images, number of stacks and transforms
      r['n'] = r_input['n']+r_goal['n']
      r['I'] = (r_input['n']*r_input['I'] + r_goal['n']*r_goal['I'])/r['n']
      r['shift'] = np.vstack([np.ones((r_input['n'],1))*shift + r_input['shift'] , r_goal['shift']])

    return r

  def plot_SNR(self,n):
    t_arr = np.linspace(0,8989/15,8989)

    C = self.estimates.C[n,:]+self.estimates.YrA[n,:]
    baseline=np.median(C)
    trace = C.copy()
    trace -= baseline
    trace *= -1*(trace <= 0)
    N_s = (trace>0).sum()
    noise = np.sqrt((trace**2).sum()/(N_s*(1-2/np.pi)))

    plt.figure(figsize=(4,2.5))
    ax = plt.subplot(111)
    tlen=20
    tstart=220
    tend=tstart+tlen
    #print(t_range)

    peak = C[tstart*15:tend*15].max()
    peak_pos = t_arr[np.where(C==peak)[0]]

    y_arr = np.linspace(-4*noise,baseline+4*noise,200)
    x1 = norm.pdf(y_arr,loc=baseline,scale=noise)
    x1 = x1/x1.max()*3+peak_pos
    x2 = peak_pos*np.ones(200)

    #plt.plot(x_offset,A_0,'ko')
    ax.fill_betweenx(y_arr,x1,x2,facecolor='r',alpha=1,edgecolor=None)

    ax.plot(t_arr,C,'k')
    ax.plot(peak_pos,peak,'ro',markersize='5')

    ax.annotate('SNR = $p(C|\mathcal{N}(c_0,\sigma))$',(peak_pos,peak),xytext=(peak_pos+2,peak*0.8),arrowprops={'arrowstyle':'->'},fontsize=14)

    ax.set_xlim([tstart,tend])
    ax.set_xticks(np.linspace(tstart,tend,tlen//10+1))
    ax.set_yticks([])
    ax.set_xlabel('time [s]')
    ax.set_ylabel('Ca$^{2+}$')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.show(block=False)

    ext = 'png'
    path = '/home/wollex/Data/Science/PhD/Thesis/pics/Methods/SNR_explained.%s'%ext
    plt.savefig(path,format=ext,dpi=300)

  def plot(self,nCells=None):

    nCells = self.estimates.C.shape[0] if nCells is None else nCells
    plt.rcParams['font.size'] = 12

    fr = 15
    plt.figure(figsize=(7,2.5))
    t_arr = np.linspace(0,8989/15,8989)
    ax_motion = plt.axes([0.1,0.575,0.25,0.35])
    ax_motion.plot(t_arr,self.shift[:,0],'k')
    ax_motion.set_ylabel('$\Delta x$ [px]')
    ax_motion.set_xticks([])
    ax_motion.set_ylim([-50,10])
    ax_motion.spines['right'].set_visible(False)
    ax_motion.spines['top'].set_visible(False)
    ax_motion.spines['bottom'].set_visible(False)

    ax_motion2 = plt.axes([0.1,0.2,0.25,0.35])
    #ax_motion.plot(t_arr,self.shift[:,1],'r')
    ax_motion2.plot(t_arr,self.shift[:,1],'r')
    ax_motion2.set_ylabel('$\Delta y$ [px]')
    ax_motion2.set_xlabel('time [s]')
    ax_motion2.set_ylim([-50,10])
    ax_motion2.spines['right'].set_visible(False)
    ax_motion2.spines['top'].set_visible(False)
    ax_motion2.spines['bottom'].set_visible(False)

    ax_sn = plt.axes([0.475,0.65,0.25,0.3])
    ax_g = plt.axes([0.475,0.2,0.25,0.275])
    np.random.seed()
    #n_arr = [1840,1564,434,2477,2432,1441,2697,2426,800,1877]
    #n_arr = np.random.randint(0,self.estimates.C.shape[0],5)
    #print(n_arr)
    n_arr = [2370,1704,821,560,681]
    ax_sn.fill_between([0.25*fr,0.5*fr],[0,0],[0.00015],alpha=0.2)
    for n in n_arr:#np.random.randint(0,self.estimates.C.shape[0],10):
      fluor = self.estimates.C[n,:] + self.estimates.YrA[n,:]
      sn = GetSn(fluor, ax_sn, range_ff=[0.25,0.5], fr=fr, method='logmexp')
      g = estimate_time_constant(fluor, ax_g, p=1, sn=sn, lags=30, fudge_factor=1.)
      #print(-1/(fr*np.log(g)))

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
    ax_sn.set_ylim([-0.00001,0.000175])
    ax_sn.spines['right'].set_visible(False)
    ax_sn.spines['top'].set_visible(False)

    ax_g.set_xlabel('$\Delta t$ [s]')
    ax_g.set_yscale('log')
    ax_g.set_ylabel('Cov$_{y}$')
    ax_g.set_yticks([])
    ax_g.set_ylim([10**(-6),10**(-3)])
    ax_g.spines['right'].set_visible(False)
    ax_g.spines['top'].set_visible(False)

    ax = plt.axes([0.85,0.2,0.125,0.75])
    ax2 = ax.twiny()
    ax3 = ax.twinx()
    ax2.hist(-1/(fr*np.log(g)),np.linspace(0,2,51),facecolor='b',alpha=0.5,density=True,orientation='horizontal')
    ax2.invert_xaxis()
    ax2.set_xticks([])
    ax3.hist(sn,np.linspace(0,1,51),facecolor='r',alpha=0.5,density=True)
    ax3.set_yticks([])
    ax3.invert_yaxis()
    ax.plot(sn,-1/(fr*np.log(g)),'.',markersize=1,markeredgewidth=0,markerfacecolor='k')
    ax.set_xlabel('noise/median')
    ax.set_ylabel('$\\tau_{\gamma}$ [s]')
    ax.set_ylim([0,2])
    ax.set_xlim([0,1])

    plt.tight_layout()
    plt.show(block=False)

    ext = 'png'
    path = '/home/wollex/Data/Science/PhD/Thesis/pics/Methods/preCNMF.%s'%ext
    plt.savefig(path,format=ext,dpi=300)




class post_CNMF:
  def __init__(self,basePath,mouse,s,dataSet='OnACID',SNR_thr=2,rval_thr=0):
    pxtomu = 530.68/512

    pathMouse = pathcat([basePath,mouse])
    pathSession = pathcat([pathMouse,'Session%02d'%s])
    pathData = pathcat([pathSession,'results_OnACID.mat'])
    ld = loadmat(pathData,variable_names=['Cn'])

    pathData = pathcat([pathSession,'results_%s.mat'%dataSet])
    self.data = loadmat(pathData,variable_names=['A','C','SNR','r_values'])
    self.data['Cn'] = ld['Cn'].T

    self.idx_good = (self.data['SNR']>2).T & (self.data['r_values']>0).T
    self.cm = com(self.data['A'],512,512) * pxtomu
    self.cm[np.squeeze(~self.idx_good),:] = np.NaN
    self.D_ROIs = sp.spatial.distance.pdist(self.cm)
    self.D_ROIs_mat = sp.spatial.distance.squareform(self.D_ROIs)


    nCells = self.data['C'].shape[0]
    self.corrcoef = np.zeros((nCells,nCells))
    for n in tqdm(np.where(self.idx_good)[0]):
      for nn in range(n+1,nCells):
          if self.idx_good[nn]:
              self.corrcoef[n,nn] = np.corrcoef(self.data['C'][n,:],self.data['C'][nn,:])[0,1]

    self.Acorr = np.zeros((nCells,nCells))
    for n in tqdm(np.where(self.idx_good)[0]):#range(nCells)):
      idx_close=np.where(self.D_ROIs_mat[n,:]<20)[0]
      for nn in idx_close:
        self.Acorr[n,nn],_ = calculate_img_correlation(self.data['A'][:,n],self.data['A'][:,nn],shift=False)

    self.Asz = (self.data['A']>0).sum(0)

  def plot(self,nCells=None):

    rc('font',size=10)
    rc('axes',labelsize=12)
    rc('xtick',labelsize=8)
    rc('ytick',labelsize=8)

    nCells = self.data['C'].shape[0] if nCells is None else nCells
    plt.figure(figsize=(7,5),dpi=300)

    idx_border = np.any(self.cm < 5,1) | np.any(self.cm > 507,1)
    idx_good = np.squeeze((~np.isnan(self.cm[:,0])) & (self.data['SNR']>3) & (self.data['r_values']>0)) & (~idx_border)
    idx_decent = np.squeeze((~np.isnan(self.cm[:,0])) & (self.data['SNR']>2) & (self.data['SNR']<3) & (~idx_border) & (self.data['r_values']>0))

    ### activity distribution
    ax_ROI = plt.axes([0.05,0.425,0.4,0.55])
    ax_ROI.imshow(self.data['Cn'],origin='lower',clim=[self.data['Cn'].min(),self.data['Cn'].max()])#,cmap='viridis')
    [ax_ROI.contour((a/a.max()).reshape(512,512).toarray(), levels=[0.3], colors='k', linewidths=[0.3]) for a in self.data['A'][:,idx_good].T]
    [ax_ROI.contour((a/a.max()).reshape(512,512).toarray(), levels=[0.3], colors='k', linewidths=[0.3], linestyles=['dotted']) for a in self.data['A'][:,idx_decent].T]
    sbar = ScaleBar(530.68/512 *10**(-6),location='lower right')
    ax_ROI.add_artist(sbar)
    ax_ROI.set_xticks([])
    ax_ROI.set_yticks([])

    ### Calcium trace examples
    ax_Ca = plt.axes([0.55,0.675,0.4,0.3])
    t_arr = np.linspace(0,8989/15,8989)

    N = 6
    color = iter(plt.cm.rainbow(np.linspace(0,1,N)))
    # for i,n in enumerate(np.random.randint(0,self.data['C'].shape[0],N)):
    for i,n in enumerate(np.random.choice(np.where(idx_good)[0],N,replace=False)):
      col = next(color)
      c = self.data['C'][n,:]
      p = ax_Ca.plot(t_arr,c/c.max()+i,linewidth=0.5)#,c=col,
      col = p[0].get_color()
      a = self.data['A'][:,n].reshape(512,512).toarray()
      ax_ROI.contour(a/a.max(),colors=[col],levels=[0.3],linewidths=[1.3])#
    ax_Ca.set_xlabel('t [s]')
    ax_Ca.set_ylabel('Ca$^{2+}$ activity')
    ax_Ca.set_yticks([])
    ax_Ca.spines['right'].set_visible(False)
    ax_Ca.spines['top'].set_visible(False)

    ax_Asz = plt.axes([0.05,0.1,0.25,0.25])
    ax_Asz.hist(self.Asz.flat,np.linspace(75,450,51),facecolor='k')
    ax_Asz.set_xlabel('size (px)')
    ax_Asz.set_ylabel('count')
    ax_Asz.set_yticks([])
    ax_Asz.spines['right'].set_visible(False)
    ax_Asz.spines['top'].set_visible(False)

    #### SNR, r value, CNN stuff
    # idxes = np.triu_indices(self.D_ROIs_mat.shape[0])
    idxes = np.ones(self.D_ROIs_mat.shape,'bool')
    idxes[~idx_good,:] = False
    idxes[:,~idx_good] = False
    idxes[np.tril_indices(self.D_ROIs_mat.shape[0])] = False


    #### distance vs correlation plot
    ax_D = plt.axes([0.55,0.41,0.25,0.1])
    ax_D.hist(self.D_ROIs_mat[idxes],np.linspace(0,700,701),facecolor='k')
    ax_D.set_ylabel('count')
    ax_D.set_yticks([])
    # ax_D.set_xlabel('distance [$\mu$m]')
    ax_D.set_xlim([0,700])
    ax_D.spines['right'].set_visible(False)
    ax_D.spines['top'].set_visible(False)

    #plt.rcParams['scatter.marker'] = '.'
    ax_D_corr = plt.axes([0.55,0.1,0.25,0.25])
    ax_D_corr.plot(self.D_ROIs_mat[idxes],self.corrcoef[idxes],'.',color='k',markersize=1,markeredgewidth=0,markerfacecolor='k')

    nsteps = 701
    D_arr = np.linspace(0,700,nsteps)
    corr_mean = np.zeros(nsteps)*np.NaN
    for i in range(nsteps-1):
        idx_D = (self.D_ROIs_mat > D_arr[i]) & (self.D_ROIs_mat < D_arr[i+1])
        corr_mean[i] = self.corrcoef[idxes & idx_D].mean()
    ax_D_corr.plot(D_arr,corr_mean,'r-',label='$\left\langle c_{Ca} \\right\\rangle (d)$')
    ax_D_corr.tick_params(axis='y',which='both',left=True,right=True)
    ax_D_corr.tick_params(axis='x',which='both',bottom=True,top=True,labelbottom=True,labeltop=False)
    ax_D_corr.set_xlabel('distance [$\mu$m]')
    ax_D_corr.legend(fontsize=10,loc='upper right')
    # ax_D_corr.set_ylabel('corr')
    ax_D_corr.set_xlim([0,700])
    ax_D_corr.set_ylim([-0.3,1])



    ax_corr = plt.axes([0.42,0.1,0.075,0.25])
    ax_corr.hist(self.corrcoef[idxes],np.linspace(-0.3,1,101),facecolor='k',orientation='horizontal')
    ax_corr.invert_xaxis()
    ax_corr.set_xticks([])
    ax_corr.tick_params(axis='y',which='both',left=False,right=True,labelright=False,labelleft=False)
    ax_corr.set_xlabel('count')
    ax_corr.set_ylim([-0.3,1])
    ax_corr.spines['left'].set_visible(False)
    ax_corr.spines['top'].set_visible(False)

    ax_Ac = plt.axes([0.85,0.41,0.075,0.1])
    ax_Ac.plot(self.Acorr[idxes],self.D_ROIs_mat[idxes],'.',color='k',markersize=1,markeredgewidth=0,markerfacecolor='k')
    ax_Ac.set_ylim([0,20])
    ax_Ac.yaxis.tick_right()
    ax_Ac.yaxis.set_label_position("right")
    ax_Ac.set_ylabel('d [$\mu$m]')
    ax_Ac.spines['top'].set_visible(False)

    ax_Acc = plt.axes([0.85,0.1,0.075,0.25])
    ax_Acc.plot(self.Acorr[idxes],self.corrcoef[idxes],'.',color='k',markersize=1,markeredgewidth=0,markerfacecolor='k')
    ax_Acc.set_ylim([-0.3,1])
    ax_Acc.yaxis.set_label_position("right")
    ax_Acc.tick_params(axis='y',which='both',left=True,right=True,labelright=False,labelleft=False)
    ax_Acc.tick_params(axis='x',which='both',bottom=True,top=True,labelbottom=True,labeltop=False)
    ax_Acc.set_ylabel('$c_{Ca}$')
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
    #gr[gr > 1] = 0.95 + np.random.normal(0, 0.01, np.sum(gr > 1))
    #gr[gr < 0] = 0.15 + np.random.normal(0, 0.01, np.sum(gr < 0))
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


def find_activity_intervals(C, Npeaks: int = 5, tB=-3, tA=10, thres: float = 0.3):# -> List:
    # todo todocument
    import peakutils
    K, T = np.shape(C)
    L: List = []
    for i in range(K):
        if np.sum(np.abs(np.diff(C[i, :]))) == 0:
            L.append([])
            logging.debug('empty component at:' + str(i))
            continue
        indexes = peakutils.indexes(C[i, :], thres=thres)
        srt_ind = indexes[np.argsort(C[i, indexes])][::-1]
        srt_ind = srt_ind[:Npeaks]
        L.append(srt_ind)

    LOC = []
    for i in range(K):
        if len(L[i]) > 0:
            interval = np.kron(L[i], np.ones(int(np.round(tA - tB)), dtype=int)) + \
                np.kron(np.ones(len(L[i]), dtype=int), np.arange(tB, tA))
            interval[interval < 0] = 0
            interval[interval > T - 1] = T - 1
            LOC.append(np.array(list(set(interval))))
        else:
            LOC.append(None)

    return LOC
