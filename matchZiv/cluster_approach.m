%%% clustering of ROIs from several sessions
%%% 
%%% program works according to several steps:
%%%   1. loading all ROI data and calculate centroid distances and footprint correlation
%%%   2. from all data: generate model of distance/footprint correlation distribution (Ziv, 2017)
%%%   ?3. on patches (due to memory): calculate network of ROIs: weights = prob of being the same neuron
%%%   ?4. on patches (due to memory): partitioning of networks into clusters (Spohrns 2012/ Newman 2006)
%%%         don't think that'll work... not made for such super sparse networks
%%%         rather: prepartition to join all >0.5 prob into one cluster, independent of double assignment, etc
%%%         than: manual refinement /check of clusters with low cluster score & resolving double assignments
%%%         -> need proper GUI for that!
%%%         -> also in there: merging / splitting of ROIs (need to keep original .h5 file??)
%%%   ?5. merging results
%%%   ?6. searching for merges / split of ROIs and applying
%%%   ?7. display results with enabling manual refinement of matching
%%%


function [data,xdata,ROI_cluster] = cluster_approach(mouse,nSes,data,xdata)
  
  pathMouse = pathcat('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/',sprintf('%d',mouse));
  
  %%% 1. loading data
  if nargin < 3
    if nargin < 2
      nSes = [];
    end
    [data, footprints] = match_loadSessions(pathMouse,nSes);
  end
  nSes = length(data);
  
  
  %%% 2. creating model of distributions
  if nargin < 4
    [xdata, histo, para] = match_analyzeData(data,nSes,12);      %% calculate distances and footprint correlation
    [model,histo] = match_buildModel(xdata,histo,para,nSes,pathMouse);
  %    [ROC] = estimate_model_accuracy(histo,model,para,pathMouse);
    
    %% and assigning probabilities to each (close) pair
    [xdata] = match_assign_prob(xdata,data,model,para);
  end
  
  
  %%% 3. initial registration to clusters
  % add data.assign
  %   .num
  %   .cluster
  %   .something??
%    match_initial_registration(xdata,data)
  disp('ROI_clustering...')
  for s = 1:nSes
    data(s).registered = false(data(s).nROI,1);
  end
  nCluster = 0;
  ROI_cluster = struct;
  tic
  for s = 1:nSes
    data(s).status = zeros(data(s).nROI,1);
    
    for sm = 1:nSes
      
      for n = 1:data(s).nROI
        
        if sm == s
          continue
        end
        
        if ~data(s).registered(n)   %% add new ROI_cluster if not already belonging to one
          nCluster = nCluster + 1;
          data(s).cluster(n).ID = nCluster;
          data(s).registered(n) = true;
          ROI_cluster(data(s).cluster(n).ID).list = zeros(nSes,1);
          ROI_cluster(data(s).cluster(n).ID).list(s) = n;
          ROI_cluster(data(s).cluster(n).ID).ct = 1;
        end
        
        cluster_candidates = find(xdata(s,sm).p_same_joint(n,:)>0.5);    %% all ROIs in sm that are candidates to be same as ROI (s,n)
        
        for m = cluster_candidates
          if ~data(sm).registered(m)
            data(sm).cluster(m).ID = data(s).cluster(n).ID;
            data(sm).registered(m) = true;
          end
          
          for c = data(s).cluster(n).ID         %% in case, there are several from one session that fall into one cluster
            if ~ismember(m,ROI_cluster(c).list(sm,:))
              idx = nnz(ROI_cluster(c).list(sm,:))+1;
              ROI_cluster(c).list(sm,idx) = m;
              if idx == 1
                ROI_cluster(c).ct = ROI_cluster(c).ct + 1;
              end
            end
          end
        end
      end
    end
  end
  
  nMatches = [ROI_cluster.ct];
  disp(sprintf('number of ROI_clusters: %d',nCluster))
  disp(sprintf('number of real ROI_clusters: %d',sum(nMatches > 1)))
  toc
  disp('pre-clustering done')
  
  
  
  extents_prototype = ones(2,2);
  extents_prototype(:,1) = data(1).imSize;
  
  for c = 1:nCluster
    %% calculate mean centroid
%      c
    [s_lst l_lst] = find(ROI_cluster(c).list);
    cluster_centroid = [];
    
    ROI_cluster(c).extents = extents_prototype;
    ROI_cluster(c).A_cluster = sparse(data(1).imSize(1),data(1).imSize(2));
    
    for i = 1:length(s_lst)
      
      n = ROI_cluster(c).list(s_lst(i),l_lst(i));
      s = s_lst(i);
      
      A_tmp = reshape(data(s).A(:,n),data(1).imSize(1),data(1).imSize(2));
      
      if ~nnz(A_tmp)
        ROI_cluster(c).list(s_lst(i),l_lst(i)) = 0;
        data(s).status(n) = NaN;
%          disp(sprintf('out: %d',c))
      else
        
        cluster_centroid = [cluster_centroid; data(s).centroid(n,:)];
        [y_idx x_idx] = find(A_tmp);
        
        %% only needed if CNMF is run on patches
        ROI_cluster(c).extents(1,1) = max(1,min(ROI_cluster(c).extents(1,1),min(y_idx)));
        ROI_cluster(c).extents(1,2) = min(data(s).imSize(1),max(ROI_cluster(c).extents(1,2),max(y_idx)));
        ROI_cluster(c).extents(2,1) = max(1,min(ROI_cluster(c).extents(2,1),min(x_idx)));
        ROI_cluster(c).extents(2,2) = min(data(s).imSize(2),max(ROI_cluster(c).extents(2,2),max(x_idx)));
        ROI_cluster(c).A_cluster = ROI_cluster(c).A_cluster + A_tmp;
      end
        
    end
    
    ROI_cluster(c).distances = squareform(pdist(cluster_centroid));
    ROI_cluster(c).distances(logical(eye(length(cluster_centroid)))) = NaN;
    ROI_cluster(c).mean_centroid = mean(cluster_centroid,1);
    
%      ROI_cluster(c).centroid = [sum((1:data(1).imSize(1))*ROI_cluster(c).A_cluster),sum(ROI_cluster(c).A_cluster*(1:data(1).imSize(2))')];
    
  end
  
  cluster.centroids = cat(1,ROI_cluster(c).mean_centroid);
  size(cluster.centroids)
  
  disp('now plotting')
  figure
  
  load('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/884/Session01/reduced_MF1_LK1.mat','max_im')
  imagesc(max_im)
  
  plt_handles = struct;
  hold on
  for c = 1:nCluster
    if ~isempty(ROI_cluster(c).mean_centroid)
      plt_handles(c).centroid = plot(ROI_cluster(c).mean_centroid(2),ROI_cluster(c).mean_centroid(1),'ko');
      
    end
  end
  hold off
  
    
%    pathData = sprintf('%s%d/cluster_registered.mat',basePath,mouse);
%    save(pathData,'ROI_cluster','data','xdata','model','para','-v7.3')
%    disp(sprintf('data saved under %s',pathData))
  
%    figure
%    histogram(nMatches)
  
  
  %%% 3.1. and display
  
  
  
  
  
  
  
  
  
  
  
%    %%% 3. on patches (due to memory): calculate network of ROIs: weights = prob of being the same neuron
%    %% -> with sparse matrices might be possible on whole -> however, shouldn't give any advantage, just slows it down (eigenvalue computation...)
%    
%    
%    
%    %% create patches
%    patch_overlap = 8;
%    npatches = 4;
%    np = sqrt(npatches);
%    [X_patches,Y_patches] = meshgrid(linspace(0,data(1).imSize(2),np+1),linspace(0,data(1).imSize(2),np+1));
%    
%    for x = 1:np
%      
%      for y = 1:np
%        
%        X_patch = [X_patches(y,x)-patch_overlap, X_patches(y,x+1)+patch_overlap];
%        Y_patch = [Y_patches(y,x)-patch_overlap, Y_patches(y+1,x)+patch_overlap];
%        
%        centroid_patch = [];
%        
%        for s = 1:nSes
%          
%          inpatch_x = all([data(s).centroid(:,2) < X_patch(2),data(s).centroid(:,2) > X_patch(1)],2);
%          inpatch_y = all([data(s).centroid(:,1) < Y_patch(2),data(s).centroid(:,1) > Y_patch(1)],2);
%          
%          inpatch = all([inpatch_x,inpatch_y],2);
%          
%          centroid_patch = [centroid_patch;data(s).centroid(inpatch,:)];
%        end
%        
%  %        size(centroid_patch)
%        subplot(np,np,x+(np-y)*np)
%        scatter(centroid_patch(:,2),centroid_patch(:,1),5,'ko')
%        xlim([X_patch(1),X_patch(2)])
%        ylim([Y_patch(1),Y_patch(2)])
%        
%        dist_arr = squareform(pdist(centroid_patch));   %% how to properly access by indices?
%        
%  %        figure
%  %        imagesc(dist_arr)
%  %        colorbar
%        
%          
%        
%      end
%    end
%    %% iterate through 
%  %    for s = 1:nSes
%    
  
  
  
  
  
  
end