%% check validity of ROI detection algorithm

function check_cluster(cluster,data,xdata,model,para,p_thr,idx)
  
  close all
  %% pick a random cluster,if none given
  if nargin < 5 || isempty(idx)
    idx = randi(length(cluster))
  end
  disp(sprintf('displaying cluster #%d',idx))
  
  cluster(idx).list
%    cluster(idx).list = cluster(idx).list(:,1);
%    cluster(idx).list
  
  %% if more than one entry:
  if length(find(cluster(idx).list)) > 1
    
    %% initialize:
    
    nSes = size(cluster(idx).list,1);
    width = size(cluster(idx).list,2);
    
    % get positions of entries in cluster list
    [y_idx x_idx] = find(cluster(idx).list);
    entries = length(y_idx);
    
    
    res = struct;
    
    res.corr = zeros(nSes,nSes);
    res.dist = zeros(nSes,nSes);
    res.prob = zeros(nSes,nSes);
    res.fp_corr_oneway = zeros(nSes,nSes,2);
    
    res.corr(:) = nan;
    res.dist(:) = nan;
    res.prob(:) = nan;
    res.fp_corr_oneway(:) = NaN;
    
    for i = 1:entries;
      s = y_idx(i);
      l = x_idx(i);
      n = cluster(idx).list(s,l);
      for w = 1:width
        for sm=1:nSes
          m = cluster(idx).list(sm,w);
          if m && sm~=s
            res.corr(s,sm) = xdata(s,sm).fp_corr(n,m);
            res.dist(s,sm) = xdata(s,sm).dist(n,m);
            res.prob(s,sm) = xdata(s,sm).p_same_joint(n,m);
            
            res.fp_corr_oneway(s,sm,1) = get_1w_corr(data,[s,n],[sm,m]);
            res.fp_corr_oneway(s,sm,2) = get_1w_corr(data,[sm,m],[s,n]);
          end
        end
      end
    end
    
    cleanup = true;
    
%      plot_check_cluster(cluster(idx),res,data,nSes,17);
    
    ROI_score_old = get_ROI_score(res,0)
    plot_ROIs(cluster(idx),data,res,0)
    
    %% first, check cluster removal, then check single neurons
%      s_rm = find(nanmean(res.prob) < 0.5)
%      ROI_score_new = get_ROI_score(res,s_rm)
%      disp(sprintf('removing ROIs from cluster to enhance cluster score from %4.2g to %4.2g',ROI_score_old,ROI_score_new))
%      s_rm
%      plot_ROIs(cluster(idx),data,res,s_rm)
%      plot_check_cluster(cluster(idx),res,data,nSes,0);
%      for s=s_rm
%        res.prob(s,:) = NaN;
%        res.prob(:,s) = NaN;
%        res.fp_corr_oneway(s,:,:) = NaN;
%        res.fp_corr_oneway(:,s,:) = NaN;
%        cluster(idx).list(s,:) = 0;
%      end
    
    while cleanup
      ROI_score_old = get_ROI_score(res,0);
      disp(sprintf('old ROI score: %6.4g',ROI_score_old))
      
      ROI_score_test = zeros(nSes,1);
      for s = 1:nSes
        if nnz(cluster(idx).list(s,:)) > 1  %% something's wrong with merge status
          disp('test for merging...')
          %% should test effect of merging vs removing extra ROIs (or all)
          ROI_score_new = get_ROI_score_merge(cluster(idx),data,model,para,res,s)   %% merge data.A from both ROIs, recalculate prob and fp_1way and calculate score (by calling get_ROI_score with temporary res struct)
          disp(sprintf('new ROI score: %6.4g vs %6.4g',ROI_score_new,ROI_score_old))
          if ROI_score_new > ROI_score_old
            disp(sprintf('merging ROI %d & %d in session #%d from cluster to enhance cluster score from %6.4g to %6.4g',cluster(idx).list(s,1),cluster(idx).list(s,2),s,ROI_score_old,ROI_score_new))
            plot_ROIs(cluster(idx),data,res,s)
          end
          
          
%            disp('pause')
%            pause(5)
          
        else
          ROI_score_test(s) = get_ROI_score(res,s);
        end
      end
      
      [ROI_score_new,s_rm] = max(ROI_score_test);
      disp(sprintf('new ROI score: %6.4g (when removing ROI from session %d)',ROI_score_new,s_rm))
      
      
      if ROI_score_new > ROI_score_old
        
        n_rm = cluster(idx).list(s_rm,1);
%          [s_rm n_rm]
%          plot_ROIs(cluster(idx),data,res,s_rm)
%          plot_check_cluster(cluster(idx),res,data,nSes,s_rm);
        res.prob(s_rm,:) = NaN;
        res.prob(:,s_rm) = NaN;
        res.fp_corr_oneway(s_rm,:,:) = NaN;
        res.fp_corr_oneway(:,s_rm,:) = NaN;
        cluster(idx).list(s_rm,:) = 0;
        
        disp(sprintf('removing ROI in session #%d from cluster to enhance cluster score from %4.2g to %4.2g',s_rm,ROI_score_old,ROI_score_new))
      else
        cleanup=false;
      end
%        cleanup=false;
      
    end
    plot_ROIs(cluster(idx),data,res,0)
    
    
  end
  
end



function [corr_1w] = get_1w_corr(data,N,M)
  s = N(1);
  n = N(2);
  sm = M(1);
  m = M(2);
  A_idx = find(data(s).A(:,n));
  corr_1w = full(dot(data(s).A(A_idx,n),data(sm).A(A_idx,m))/(data(s).norm(n)*norm(data(sm).A(A_idx,m))));
end



function [ROI_score] = get_ROI_score(res,s_rm)
  
  nSes = size(res.prob,1);
  N = sum(~isnan(nanmean(res.prob)));
  
  if s_rm>0
    for s=s_rm
      res.prob(s,:) = NaN;
      res.prob(:,s) = NaN;
      res.fp_corr_oneway(s,:,:) = NaN;
      res.fp_corr_oneway(:,s,:) = NaN;
    end
    N=N-length(s_rm);
  end
  
  p_mean = nanmean(res.prob(:));
  p_var = nanvar(res.prob(:));
  p_thr = min(0.9,p_mean-2*p_var);
  
  fp_mean = squeeze(nanmean(res.fp_corr_oneway,1));
  fp_max = nanmax(fp_mean,[],2);
  
  frac_1w_mean = nanmean(fp_max);
  frac_1w_var = nanvar(fp_max);
  
%    ROI_score = 1/3*(p_mean^(1+p_var) + frac_1w_mean^(1+frac_1w_var) + nanmin(res.prob(:)))*N/nSes;
%    ROI_score = (p_mean^(1+p_var) * frac_1w_mean^(1+frac_1w_var) * nanmin(res.prob(:)))^(1/3)*N/nSes;
  ROI_score = 1/5*(p_mean^(1+p_var) + frac_1w_mean^(1+frac_1w_var) + nanmin(res.prob(:)) + 2*N/nSes);
  
%    disp(sprintf('s: %d ROI score: %6.4g',s_rm(1), ROI_score))
end



function [ROI_score] = get_ROI_score_merge(cluster,data,model,para,res,s_merge)
  
  microns_per_pixel = 530.684/512;
  
  imSize = [512,512];
  nSes = size(res.prob,1);
  %% merging data.A together (all combinations, if >2 exist?
  A_merge = sparse(size(data(s_merge).A(:,1),1),1);
  
  for n = cluster.list(s_merge,:)
    if n > 0
      A_merge = A_merge + data(s_merge).A(:,n);
    end
    %% obtain reference values of prob etc?
  end
  A_test = bwconncomp(full(reshape(A_merge,imSize(1),imSize(2)))>0,8);
  if A_test.NumObjects > 1
    disp('dont merge, ROIs are not connected!')
%      A_test.NumObjects
    merge_status = false;
%      return
  end
  
  A_area = nnz(A_merge);
  
%    disp(sprintf('area of the merged ROI: %d',nnz(A_merge)))
  
  A_norm = norm(A_merge);
  A_tmp = reshape(A_merge/sum(A_merge(:)),imSize(1),imSize(2));
  
  A_centroid = zeros(nSes,2);
  for s = 1:nSes
    %% compare to the ROI with highest matching probability (for now) - is it always #1?
    n = cluster.list(s,1);
    
    A_centroid(s,:) = [sum((1:imSize(1))*A_tmp),sum(A_tmp*(1:imSize(2))')];
    
    %% recalculating res. structure values (prob & fp_corr_oneway)
    dist(s) = microns_per_pixel*sqrt((A_centroid(s,1) - data(s).centroid(n,1)).^2 + (A_centroid(s,2) - data(s).centroid(n,2)).^2);
    fp_corr(s) = full(dot(A_merge,data(s).A(:,n))/(A_norm*data(s).norm(n)));
    
    idx_dist = max(1,ceil(para.nbins*dist(s)/para.dist_max));
    idx_corr = max(1,ceil(para.nbins*fp_corr(s)/para.corr_max));
    
    res.prob(s,s_merge) = model.p_same_joint(idx_dist,idx_corr);
    res.prob(s_merge,s) = model.p_same_joint(idx_dist,idx_corr);
    
    %% now add 1w correlation
  end
  
%    [fp_corr',res.corr(s_merge,:)']
  ROI_score = get_ROI_score(res,s_merge);

end



function plot_ROIs(cluster,data,res,s_rm)
  
%    
  
  [y_idx x_idx] = find(cluster.list);
  entries = length(y_idx);
  nSes = size(res.prob,1);
  
  figure('position',[100 100 700 500])
%    ax1 = subplot(5,3,[1,4]);
  hold on
  
  check_idx = [];
  for i = 1:entries
    s = y_idx(i);
    l = x_idx(i);
    c = ones(3,1)*4*s/(5.*nSes);
    n = cluster.list(s,l);
    A_tmp = reshape(data(s).A(:,n),512,512);
    
%        labelstr = sprintf('s=%02d,   corr=%4.2f/%4.2f,   dist=%4.2f/%4.2f,   p_{same}=%4.2f/%4.2f',s,cluster(idx).fp_corr(s,l),min_corr(s),cluster(idx).dist(s,l),max_dist(s),cluster(idx).p_same(s,l),min_prob(s));
    labelstr = sprintf('s=%02d',s);
%      
%        disp(sprintf('skipping neuron %d from session %d',n,s))
%      if all(isnan(res.prob(s,:)))
%        disp(sprintf('not plotting session %d',s))
%      elseif nanmean(res.prob(s,:))<p_thr%ismember(s,s_red)

    if ismember(s,s_rm)
      contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color','r','DisplayName',labelstr)
    else
      contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color',c,'DisplayName',labelstr)
    end
  end
  hold off
  
  r = 12;
  xlim([1,512])
  ylim([1,512])
  
  if s_rm(1) > 0
    s=s_rm(1);
  else
    s = y_idx(1);
  end
  
  n = cluster.list(s,1);
  
  xlim([data(s).centroid(n,2)-r data(s).centroid(n,2)+r])
  ylim([data(s).centroid(n,1)-r data(s).centroid(n,1)+r])
  
  legend('Location','eastoutside')
  
end




function plot_check_cluster(cluster,res,data,nSes,s_red)
  
  cmrg = [linspace(1,0,50)',linspace(0,1,50)',linspace(0,0,50)'];
  
  [y_idx x_idx] = find(cluster.list);
  entries = length(y_idx);
  
  res.frac_1w = res.fp_corr_oneway(:,:,1)./res.fp_corr_oneway(:,:,2);
  
  
  figure('position',[100 100 1200 1000])
  ax1 = subplot(5,3,[1,4]);
  hold on
  
  check_idx = [];
  for i = 1:entries
    s = y_idx(i);
    l = x_idx(i);
    c = ones(3,1)*4*s/(5.*nSes);
    n = cluster.list(s,l);
    A_tmp = reshape(data(s).A(:,n),512,512);
    
%        labelstr = sprintf('s=%02d,   corr=%4.2f/%4.2f,   dist=%4.2f/%4.2f,   p_{same}=%4.2f/%4.2f',s,cluster(idx).fp_corr(s,l),min_corr(s),cluster(idx).dist(s,l),max_dist(s),cluster(idx).p_same(s,l),min_prob(s));
    labelstr = sprintf('s=%02d',s);
%      
%        disp(sprintf('skipping neuron %d from session %d',n,s))
    if all(isnan(res.prob(s,:)))
      disp(sprintf('not plotting session %d',s))
    elseif ismember(s,s_red)
      contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color','r','DisplayName',labelstr)
%        check_idx = [check_idx s];
    else
      contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color',c,'DisplayName',labelstr)
    end
  end
  hold off
  
  s = y_idx(1);
  n = cluster.list(s,x_idx(1));
  r = 12;
  xlim([1,512])
  ylim([1,512])
  xlim([data(s).centroid(n,2)-r data(s).centroid(n,2)+r])
  ylim([data(s).centroid(n,1)-r data(s).centroid(n,1)+r])
  legend('Location','eastoutside')
  
  clims = [0.5,1];
  ax2 = subplot(5,3,[2,5]);
  imagesc(res.dist)
  colormap(cmrg)
  set(gca,'clim',[0,5])
  colorbar('peer',gca)
  title('centroid distance')
  
  ax3 = subplot(5,3,[3,6]);
  imagesc(res.corr)
  colormap(cmrg)
  set(gca,'clim',clims)
  colorbar('peer',gca)
  title('footprint correlation')
  
  ax4 = subplot(5,3,[7,10]);
  imagesc(res.prob)
  colormap(cmrg)
  set(gca,'clim',clims)
  colorbar('peer',gca)
  title('p_{same}')
  
  subplot(5,3,[8,11])
  imagesc(res.fp_corr_oneway(:,:,1))
  colormap(cmrg)
  set(gca,'clim',[0.5,1])
%    colorbar('peer',gca)
  title('footprint correlation (1-way)')
  
  subplot(5,3,[9,12])
  imagesc(res.fp_corr_oneway(:,:,2))
  colormap(cmrg)
  set(gca,'clim',[0.5,1])
%    colorbar('peer',gca)
  title('footprint correlation (1-way)')
  
  subplot(5,3,13)
  hold on
  plot(nanmean(res.frac_1w),'kx')
%    plot(nanmean(res.prob),'kx','DisplayName','p_{same}')
%    plot(nanmean(res.fp_corr_oneway(:,:,1)),'rD','DisplayName','1-way fp')
%    plot(nanmean(res.fp_corr_oneway(:,:,2)),'rs','DisplayName','1-way fp')
%    plot(nanmin(res.prob),'gx','DisplayName','min')
  hold off
%    legend('Location','southeast')
  xlim([1,nSes])
  ylim([0.6 1.3])
  xlabel('Session #')
  ylabel('fraction oneway fp correlation')
  
  
  subplot(5,3,14)
  hold on
  plot(nanmean(res.prob),'kx','DisplayName','p_{same}')
  plot(nanmean(res.fp_corr_oneway(:,:,1)),'rD','DisplayName','1-way fp')
  plot(nanmean(res.fp_corr_oneway(:,:,2)),'rs','DisplayName','1-way fp')
  plot(nanmin(res.prob),'gx','DisplayName','min')
  hold off
  legend('Location','southeast')
  xlim([1,nSes])
  ylim([0.5 1])
  xlabel('Session #')
  ylabel('average fp correlation')
  
  subplot(5,3,15)
  hold on
  plot(nanvar(res.prob),'kx','DisplayName','p_{same}')
  plot(nanvar(res.fp_corr_oneway(:,:,1)),'rD','DisplayName','1-way fp')
  plot(nanvar(res.fp_corr_oneway(:,:,2)),'rs','DisplayName','1-way fp')
  hold off
  legend()
  xlim([1,nSes])
  ylim([0 0.1])
  xlabel('Session #')
  ylabel('variance fp correlation')
  
  
  
end