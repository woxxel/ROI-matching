%% check validity of ROI detection algorithm

function check_cluster(cluster,data,xdata,p_thr,idx,idx_red)
  
  close all
  if nargin < 5 || isempty(idx)
    idx = randi(length(cluster))
  end
  
  disp(sprintf('displaying cluster #%d',idx))
  cluster(idx).list
  cluster(idx).list = cluster(idx).list(:,1);
  cluster(idx).list
  
  if length(find(cluster(idx).list)) > 1
    
    [y_idx x_idx] = find(cluster(idx).list);
    
    if nargin==6
      s_red = idx_red;
    else
      s_red = 0;
    end
    
    entries = length(y_idx);
    
    nSes = size(cluster(idx).list,1);
    width = size(cluster(idx).list,2);
    
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
            
            A_idx = find(data(s).A(:,n));
            res.fp_corr_oneway(s,sm,1) = full(dot(data(s).A(A_idx,n),data(sm).A(A_idx,m))/(data(s).norm(n)*norm(data(sm).A(A_idx,m))));
            
            A_idx = find(data(sm).A(:,m));
            res.fp_corr_oneway(s,sm,2) = full(dot(data(sm).A(A_idx,m),data(s).A(A_idx,n))/(data(sm).norm(m)*norm(data(s).A(A_idx,n))));
          end
        end
      end
    end
    
    cleanup = true;
    [check_idx,res] = plot_check_cluster(cluster(idx),res,data,nSes,0);
%      size(res.prob)
    
    get_ROI_score(res,0)
    %% first, check cluster removal, then check single neurons
    s_rm = find(nanmean(res.prob) < 0.5)
    get_ROI_score(res,s_rm)
    
    while cleanup
      ROI_score_old = get_ROI_score(res,0);
      disp(sprintf('old ROI score: %6.4g',ROI_score_old))
      
      ROI_score_test = zeros(nSes,1);
      for s = 1:nSes
        ROI_score_test(s) = get_ROI_score(res,s);
      end
      
      [ROI_score_new,s_rm] = max(ROI_score_test);
      disp(sprintf('new ROI score: %6.4g (when removing ROI from session %d)',ROI_score_new,s_rm))
      
      
      
      if ROI_score_new > ROI_score_old*1.01
        
        n_rm = cluster(idx).list(s_rm,1);
        [s_rm n_rm]
        plot_ROIs(cluster(idx),data,res,s_rm,n_rm)
        
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
    [check_idx,res] = plot_check_cluster(cluster(idx),res,data,nSes,0);
%      check_idx_old = check_idx;
%      while ~isempty(check_idx)
      
      %% here, check impact of the removal of each single ROI on overall cluster
      %% accept only removal, if overall ROI score is enhanced! (and only the "best" impact at each time)
%        for c_idx = check_idx_old
        
        
%        end
%        check_idx_old = check_idx; 
%      end
%      [check_idx,res] = plot_check_cluster(cluster(idx),res,data,nSes,check_idx)
%      [check_idx,res] = plot_check_cluster(cluster(idx),res,data,nSes,check_idx)
    
%      sprintf('mean correlation %4.2f',nanmean(res.prob(:)))
%      sprintf('variance %6.4f',nanvar(res.prob(:)))
    
    
  end
  
end



function [ROI_score] = get_ROI_score(res,s_rm)
  
  nSes = size(res.prob,1);
  mask = true(nSes);
  N = sum(~isnan(nanmean(res.prob)));
  
  if s_rm>0
%      mask(s_rm,:) = false;
%      mask(:,s_rm) = false;
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
%    disp(sprintf('mean correlation %4.2f, variance %6.4f, threshold: %4.2f',p_mean,p_var, p_thr))
  
  fp_corr_1w_1 = res.fp_corr_oneway(:,:,1);
  fp_corr_1w_2 = res.fp_corr_oneway(:,:,2);
  
  res.frac_1w = fp_corr_1w_1./fp_corr_1w_2;
  res.frac_1w(res.frac_1w>1) = 1./res.frac_1w(res.frac_1w>1);
  frac_1w_mean = nanmean(res.frac_1w(:));
  frac_1w_var = nanvar(res.frac_1w(:));
  
%    ROI_score_raw = 1/3*(p_mean^(1+p_var) + frac_1w_mean^(1+frac_1w_var) + nanmin(res.prob(:)));
  [p_mean^(1+p_var) frac_1w_mean^(1+frac_1w_var) nanmin(res.prob(:))]
%    ROI_score = 1/3*(p_mean^(1+p_var) + frac_1w_mean^(1+frac_1w_var) + nanmin(res.prob(:)))*N/nSes;
  
  ROI_score = (p_mean^(1+p_var) * frac_1w_mean^(1+frac_1w_var) * nanmin(res.prob(:)))^(1/3)*N/nSes;
  
  disp(sprintf('ROI score: %4.2g',ROI_score))
%    res.prob
end


function plot_ROIs(cluster,data,res,s_rm,n_rm)
  
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

    if (s_rm==s && n_rm==n)
      contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color','r','DisplayName',labelstr)
    else
      contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color',c,'DisplayName',labelstr)
    end
  end
  hold off
  
  r = 12;
  xlim([1,512])
  ylim([1,512])
  xlim([data(s_rm).centroid(n_rm,2)-r data(s_rm).centroid(n_rm,2)+r])
  ylim([data(s_rm).centroid(n_rm,1)-r data(s_rm).centroid(n_rm,1)+r])
  legend('Location','eastoutside')
  
end




function [check_idx,res] = plot_check_cluster(cluster,res,data,nSes,check_idx_old)
  
  cmrg = [linspace(1,0,50)',linspace(0,1,50)',linspace(0,0,50)'];
  
  [y_idx x_idx] = find(cluster.list);
  entries = length(y_idx);
  
  for s = 1:nSes
    if ismember(s,check_idx_old)
      res.prob(s,:) = nan;
      res.prob(:,s) = nan;
      res.dist(s,:) = nan;
      res.dist(:,s) = nan;
      res.corr(s,:) = nan;
      res.corr(:,s) = nan;
      res.fp_corr_oneway(s,:,:) = nan;
      res.fp_corr_oneway(:,s,:) = nan;
    end
  end
  p_mean = nanmean(res.prob(:));
  p_var = nanvar(res.prob(:));
  p_thr = p_mean-4*p_var;
  p_thr = 0.95;
  
  disp(sprintf('mean correlation %4.2f, variance %6.4f, threshold: %4.2f',p_mean,p_var, p_thr))
  
%    fp_1w_mean = nanmean(reshape(res.fp_corr_oneway(:,:,1),[],1));
%    fp_1w_var = nanvar(reshape(res.fp_corr_oneway(:,:,1),[],1));
  res.frac_1w = res.fp_corr_oneway(:,:,1)./res.fp_corr_oneway(:,:,2);
  
%    1/3*(p_mean^(1+p_var) + fp_1w_mean^(1+fp_1w_var) + nanmin(res.prob(:)))
%    ROI_score = 1/3*(p_mean^(1+p_var) + res.frac_1w + nanmin(res.prob(:)))*sum(~isnan(nanmean(res.prob)))/nSes;
  
  
  
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
%      elseif nanmean(res.prob(s,:))<p_thr%ismember(s,s_red)
%        contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color','r','DisplayName',labelstr)
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