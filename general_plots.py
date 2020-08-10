import cv2, time
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sstats
from scipy import ndimage
from scipy.stats.distributions import t
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

from scipy.io import loadmat
from utils import fit_plane, z_from_point_normal_plane, rotation_matrix, bootstrap_data, calculate_img_correlation, get_shift_and_flow, pathcat

pathPlots = "/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Results/pics/general_plots/"


def plot_optical_flow(pathMouse,sessions=None):

  ## and display
  dims=(512,512)
  pathSession1 = '%sSession%02d/results_OnACID.mat' % (pathMouse,sessions[0])
  ROIs1_ld = loadmat(pathSession1)
  pathSession2 = '%sSession%02d/results_OnACID.mat' % (pathMouse,sessions[1])
  ROIs2_ld = loadmat(pathSession2)

  shift,flow,(x_grid,y_grid) = get_shift_and_flow(ROIs1_ld['A'],ROIs2_ld['A'],dims,projection=1,plot_bool=True)
  print('shift: ',shift)

  angles = find_and_rotate_from_flow(flow)

  #reverse_flow = cv2.calcOpticalFlowFarneback(Cn2_norm, Cn_norm, None,0.5,5,128,3,7,1.5,0)
  #print(reverse_flow.shape)

def find_and_rotate_from_flow(flow):
  dims = flow.shape
  x = np.hstack([np.ones((dims[0],1)),np.arange(dims[0]).reshape(dims[0],1)])
  x_grid, y_grid, z_grid = np.meshgrid(np.arange(0., dims[0]).astype(np.float32),
                                       np.arange(0., dims[1]).astype(np.float32),
                                       np.arange(0,1).astype(np.float32))

  W = sstats.norm.pdf(range(dims[0]),dims[0]/2,dims[0]/(2*1.96))
  W /= W.sum()
  W = np.sqrt(np.diag(W))
  x_w = np.dot(W,x)
  flow_w = np.dot(flow[:,:,1],W)
  x0,res,rank,s = np.linalg.lstsq(x_w,flow_w)

  d = -x0[0,:]/x0[1,:]
  #W = np.sqrt(np.diag(1/res))
  r = sstats.linregress(range(dims[0]),d)

  tilt_ax = r.intercept+r.slope*range(512)

  dist_mat = np.abs((r.slope*x_grid[:,:,0]-y_grid[:,:,0]+r.intercept)/np.sqrt(r.slope**2+1**2))
  slope_normal = np.array([-r.slope,1])
  slope_normal /= np.linalg.norm(slope_normal)
  f_perp = np.dot(flow[:,:,:2],slope_normal)

  ## need two cases to capture both, flows away from and towards axis (rotation to a more even vs rotation to a more skewed plane) - not entirely precise, only holds for relatively small angles theta
  h_dat = np.sign(f_perp)*np.sin(np.arccos((dist_mat - np.abs(f_perp))/dist_mat))*dist_mat

  data = np.stack([x_grid[...,0],y_grid[...,0],h_dat],2)
  data = data.reshape(dims[0]*dims[1],3)
  mask = ~np.isnan(data[:,2])

  #(p,n),tmp = fit_plane(data[mask,:],anchor=[dims[0]/2-1,tilt_ax[int(dims[0]/2-1)],0])
  (p,n),tmp = fit_plane(data[mask,:])
  h_plane = z_from_point_normal_plane(x_grid[...,0],y_grid[...,0],p,n)

  #par,par_std = bootstrap_data(fit_plane,data[mask,:],100)

  ## get axis (vanishing z-component) from this:
  x_ax = np.linspace(0,dims[0]-1,dims[0])
  y_ax = n[0]/n[1]*(p[0]-x_ax) + p[1] + n[2]/n[1]*p[2]

  intercept = n.dot(p)
  slope = -n[0]/n[1]

  dist_mat = np.abs((slope*x_grid[:,:,0]-y_grid[:,:,0]+intercept)/np.sqrt(slope**2+1**2))
  slope_normal = np.array([-slope,1])
  slope_normal /= np.linalg.norm(slope_normal)
  f_perp = np.dot(flow[:,:,:2],slope_normal)

  T = np.eye(4)
  T[:3,3] = -p


  angles = np.arccos(n)/(2*np.pi)*360
  print('angles:',angles)

  #print(par[1])
  #angles_bs = np.arccos(par[1])/(2*np.pi)*360
  #angles_std = par_std[1]/np.sqrt(1-par[1]**2)

  #print('angles_bs: ',angles_bs)
  #print('angles_std: ',angles_std)


  plt.figure()
  ax = plt.subplot(221)
  im = ax.imshow(h_dat,origin='lower',cmap='jet',clim=[-60,60])
  plt.colorbar(im)
  ax.plot(d,'b:')
  ax.plot(tilt_ax,'b--')
  ax.plot(x_ax,y_ax,'b-')
  ax.set_xlim([0,dims[0]])
  ax.set_ylim([0,dims[1]])

  ax2 = plt.subplot(222)
  im2 = ax2.imshow(h_plane,origin='lower',cmap='jet',clim=[-60,60])
  ax2.plot(d,'b:')
  ax2.plot(tilt_ax,'b--')
  ax2.plot(x_ax,y_ax,'b-')
  ax2.plot(p[0],p[1],'rx')
  ax2.set_xlim([0,dims[0]])
  ax2.set_ylim([0,dims[1]])
  plt.colorbar(im2)

  plt.subplot(223)
  plt.imshow(h_dat-h_plane,origin='lower',cmap='jet',clim=[-10,10])
  plt.colorbar()

  idxes = 1
  ax = plt.subplot(224,projection='3d')
  ax.plot_surface(x_grid[::50,::50,0],y_grid[::50,::50,0],h_plane[::50,::50],color='k')
  ax.plot_surface(x_grid[::idxes,::idxes,0],y_grid[::idxes,::idxes,0],h_dat[::idxes,::idxes],color='r',alpha=0.5)
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')

  plt.show(block=False)

  pl_bool=False
  if pl_bool:

    idxes=30

    T = np.eye(4)
    tx = dims[0]/2 - 1
    ty = tilt_ax[int(tx)]
    T[:3,3] = [-tx,-ty,0]

    nG = int(dims[0]/idxes +1)
    grid = np.zeros((nG,nG,4))
    grid[:,:,0] = x_grid[::idxes,::idxes,0]
    grid[:,:,1] = y_grid[::idxes,::idxes,0]
    grid[:,:,2] = z_grid[::idxes,::idxes,0]
    grid[:,:,3] = 1
    grid = np.einsum('jk,...k',T,grid)



    phi = np.arctan(slope)
    phi_deg = phi/(2*np.pi)*360
    phi_std = std_err/(2*np.pi)*360
    rot_mat = cv2.getRotationMatrix2D((dims[0]/2-1,tilt_ax[int(dims[0]/2-1)]),phi_deg,1.0)
    flow_x_rot = cv2.warpAffine(flow[:,:,0], rot_mat, dims[:2], flags=cv2.INTER_LINEAR,borderValue=np.NaN)
    flow_y_rot = cv2.warpAffine(flow[:,:,1], rot_mat, dims[:2], flags=cv2.INTER_LINEAR,borderValue=np.NaN)

    flow = np.pad(flow,[[0,0],[0,0],[0,1]])

    flow_x_rot2 = np.cos(phi)*flow_x_rot + np.sin(phi)*flow_y_rot
    flow_y_rot2 = -np.sin(phi)*flow_x_rot + np.cos(phi)*flow_y_rot

    fig = plt.figure(figsize=(12,8))
    ax_Q = plt.subplot(1,2,1,projection='3d')
    hQ = ax_Q.quiver(grid[:,:,0], grid[:,:,1], grid[:,:,2], flow[::idxes,::idxes,0], flow[::idxes,::idxes,1],flow[::idxes,::idxes,2],color='r',pivot='tail')
    ax_Q.plot(np.linspace(-tx,tx,dims[0]),d,np.zeros(512),'b:')
    ax_Q.plot(np.linspace(-tx,tx,dims[0]),tilt_ax,np.zeros(512),'b-')
    ax_Q.set_xlim([-300,300])
    ax_Q.set_ylim([-300,300])
    ax_Q.set_zlim([-300,300])

    ax_Q2 = plt.subplot(222)
    hQ2 = ax_Q2.quiver(grid[:,:,0], grid[:,:,1], flow[::idxes,::idxes,0], flow[::idxes,::idxes,1],color='r',pivot='tail')
    ax_Q2.set_xlim([-300,300])
    ax_Q2.set_ylim([-300,300])
    #print(hQ2.get_offsets())

    ax_D = plt.subplot(224)
    hD = ax_D.plot(0,0,'rx')
    ax_D.set_xlim([-180,180])
    #ax_D.set_ylim([0,10])

    anim = animation.FuncAnimation(fig, update_quiver, fargs=(hQ,hQ2,slope, grid[:,:,:3], flow, idxes),interval=10, blit=False)

    fig.tight_layout()
    plt.show(block=False)


    plt.figure(figsize=(10,12))
    plt.subplot(321)
    plt.quiver(x_grid[::idxes,::idxes], y_grid[::idxes,::idxes], flow[::idxes,::idxes,0], flow[::idxes,::idxes,1], angles='xy', scale_units='xy', scale=1, headwidth=4,headlength=4, width=0.002, units='width')
    plt.plot(d,'r--')
    plt.plot(tilt_ax,'r:')
    plt.ylim([0,dims[0]])

    plt.subplot(322)
    plt.quiver(x_grid[::idxes,::idxes], y_grid[::idxes,::idxes], flow_x_rot2[::idxes,::idxes], flow_y_rot2[::idxes,::idxes], angles='xy', scale_units='xy', scale=1, headwidth=4,headlength=4, width=0.002, units='width')
    plt.ylim([0,dims[0]])

    #plt.subplot(322)
    #plt.quiver(x_grid[::idxes,::idxes], y_grid[::idxes,::idxes], flow_rot[::idxes,::idxes,0], flow_rot[::idxes,::idxes,1], angles='xy', scale_units='xy', scale=1, headwidth=4,headlength=4, width=0.002, units='width')
    #plt.ylim([0,dims[0]])
    y_grid_dist = y_grid[:,:,0]-tilt_ax[int(dims[0]/2-1)]
    h = np.sqrt(y_grid_dist**2 - (y_grid_dist - flow_y_rot2)**2)
    h_1 = np.sqrt((y_grid_dist - flow_y_rot2)**2 - y_grid_dist**2)

    h2 = np.sin(angles[2]/360*2*np.pi)*dist_mat
    plt.subplot(323)
    plt.imshow(h,origin='lower')
    plt.colorbar()
    plt.subplot(324)
    plt.imshow(h_1,origin='lower')
    plt.colorbar()

    #h_infer = np.sqrt(y_grid_dist**2 - (y_grid_dist*np.cos(theta))**2)
    #plt.subplot(324)
    #plt.imshow(h2,origin='lower')

    #plt.imshow(h_infer,origin='lower')
    plt.colorbar()

    plt.subplot(325)
    plt.plot(alpha[1,:],'kx')
    plt.plot([0,512],np.ones(2)*np.mean(alpha[1,:]),'r-')

    #plt.subplot(326)
    #plt.imshow(h-h_infer,origin='lower')
    #plt.colorbar()

    plt.show(block=False)
  return angles

  #return (phi_deg,phi_std),(theta_deg,theta_std)


def update_quiver(num,hQ,hQ2,slope,grid,flow,idxes):

  #print('-------- phi: %5.3f ---------'%phi_arr[i])
  phi = (num+180)%360-180
  R = rotation_matrix([1,slope,0],phi,degree=True)
  flow_rot = np.einsum('jk,...k',R,flow)
  flow_sum = (flow_rot**2).sum(1).sum(0)


  grid_rot = np.einsum('jk,...k',R,grid)

  dims = grid.shape
  hQ2.set_offsets(grid_rot[:,:,:2].reshape(dims[0]*dims[1],2))
  hQ2.set_UVC(flow_rot[::idxes,::idxes,0],flow_rot[::idxes,::idxes,1])

  dim=grid.shape

  grid_plot = grid_rot.reshape(dim[0]*dim[1],dim[2])
  flow_plot = flow_rot[::idxes,::idxes,:].reshape(dim[0]*dim[1],dim[2])

  segs = [[[x[0],x[1],x[2]],[x[0]+f[0],x[1]+f[1],x[2]+f[2]]] for (x,f) in zip(grid_plot,flow_plot)]

  hQ.set_segments(segs)
  plt.title('phi = %4.2f'%phi)

  return hQ

def plot_correlation_image(pathMouse,session=None):

  plt.close('all')

  ## and display
  pathSession1 = '%sSession%02d/results_OnACID.mat' % (pathMouse,session)
  ROIs1_ld = loadmat(pathSession1)

  Cn = ROIs1_ld['Cn']
  Cn -= Cn.min()
  Cn /= Cn.max()

  plt.figure(figsize=(5,4))
  plt.imshow(Cn)
  plt.xlabel('x [px]')
  plt.ylabel('y [px]')
  plt.show(block=False)

  plt.savefig('%scorrelation_img.png'%pathPlots)

def plot_session_shifts(basePath,mouse,sessions=None,dataSet='OnACID',pltSave=False):

  pathMouse = pathcat([basePath,mouse])
  pathSession1 = pathcat([pathMouse,'Session%02d/results_%s.mat' % (sessions[0],dataSet)])
  ROIs1_ld = loadmat(pathSession1,variable_names=['A','Cn'])

  pathSession2 = pathcat([pathMouse,'Session%02d/results_%s.mat' % (sessions[1],dataSet)])
  ROIs2_ld = loadmat(pathSession2,variable_names=['A','Cn'])

  A1 = ROIs1_ld['A']
  A2 = ROIs2_ld['A']
  #dims = Cn.shape

  get_shift_and_flow(A1,A2,plot_bool=True)
  #Cn2 = cv2.remap(ROIs2_ld['Cn'].astype(np.float32), x_remap, y_remap, cv2.INTER_NEAREST)

  plt.figure()
  plt.subplot(1,2,1)
  plt.imshow(ROIs1_ld['Cn'],origin='lower')
  plt.subplot(1,2,2)
  plt.imshow(ROIs2_ld['Cn'],origin='lower')
  plt.show(block=False)

  if pltSave:
    plt.savefig('%ssession_shift.png'%pathPlots)
