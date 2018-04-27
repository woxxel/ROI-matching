

function plot_ROI_cluster(ROI_cluster,score,data,s_red)
  
  nSes = size(ROI_cluster.list,1);
  
  cmrg = [linspace(1,0,50)',linspace(0,1,50)',linspace(0,0,50)'];
  
  [y_idx x_idx] = find(ROI_cluster.list);
  entries = length(y_idx);
  
%    score.frac_1w = score.fp_corr_oneway(:,:,1)./score.fp_corr_oneway(:,:,2);
  
  
  figure('position',[100 100 800 600])
  ax1 = subplot(1,2,1);
  hold on
  
  check_idx = [];
  for i = 1:entries
    s = y_idx(i);
    l = x_idx(i);
    c = ones(3,1)*4*s/(5.*nSes);
    n = ROI_cluster.list(s,l);
    A_tmp = reshape(data(s).A(:,n),512,512);
    
%        labelstr = sprintf('s=%02d,   corr=%4.2f/%4.2f,   dist=%4.2f/%4.2f,   p_{same}=%4.2f/%4.2f',s,ROI_cluster.fp_corr(s,l),min_corr(s),ROI_cluster.dist(s,l),max_dist(s),ROI_cluster.p_same(s,l),min_prob(s));
    labelstr = sprintf('s=%02d',s);
%      
%        disp(sprintf('skipping neuron %d from session %d',n,s))
    if all(isnan(score.prob(s,:)))
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
  n = ROI_cluster.list(s,x_idx(1));
  r = 12;
  xlim([1,512])
  ylim([1,512])
  xlim([data(s).centroid(n,2)-r data(s).centroid(n,2)+r])
  ylim([data(s).centroid(n,1)-r data(s).centroid(n,1)+r])
  legend('Location','eastoutside')
  
%    clims = [0.5,1];
%    ax2 = subplot(5,3,[2,5]);
%    imagesc(score.dist)
%    colormap(cmrg)
%    set(gca,'clim',[0,5])
%    colorbar('peer',gca)
%    title('centroid distance')
  
%    ax3 = subplot(5,3,[3,6]);
%    imagesc(score.corr)
%    colormap(cmrg)
%    set(gca,'clim',clims)
%    colorbar('peer',gca)
%    title('footprint correlation')
  
%    ax4 = subplot(5,3,[7,10]);
%    imagesc(score.prob)
%    colormap(cmrg)
%    set(gca,'clim',clims)
%    colorbar('peer',gca)
%    title('p_{same}')
%    
%    subplot(5,3,[8,11])
%    imagesc(score.fp_corr_oneway(:,:,1))
%    colormap(cmrg)
%    set(gca,'clim',[0.5,1])
%  %    colorbar('peer',gca)
%    title('footprint correlation (1-way)')
%    
%    subplot(5,3,[9,12])
%    imagesc(score.fp_corr_oneway(:,:,2))
%    colormap(cmrg)
%    set(gca,'clim',[0.5,1])
%  %    colorbar('peer',gca)
%    title('footprint correlation (1-way)')
  
  
  fp_mean = squeeze(nanmean(score.fp_corr_oneway,1));
  fp_max = nanmax(fp_mean,[],2);
  
  frac_1w_mean = nanmean(fp_max);
  frac_1w_var = nanvar(fp_max);
  
  
  subplot(3,2,2)
  hold on
  plot(fp_max,'kx')
  plot(fp_mean(:,1),'rD')
  plot(fp_mean(:,2),'rs')
%    plot(nanmean(score.prob),'kx','DisplayName','p_{same}')
%    plot(nanmean(score.fp_corr_oneway(:,:,1)),'rD','DisplayName','1-way fp')
%    plot(nanmean(score.fp_corr_oneway(:,:,2)),'rs','DisplayName','1-way fp')
%    plot(nanmin(score.prob),'gx','DisplayName','min')
  hold off
%    legend('Location','southeast')
  xlim([1,nSes])
  ylim([0.5 1])
  xlabel('Session #')
  ylabel('fraction oneway fp correlation')
  
  
  subplot(3,2,4)
  hold on
  plot(nanmean(score.prob),'kx','DisplayName','p_{same}')
%    plot(nanmean(score.fp_corr_oneway(:,:,1)),'rD','DisplayName','1-way fp')
%    plot(nanmean(score.fp_corr_oneway(:,:,2)),'rs','DisplayName','1-way fp')
  plot(nanmin(score.prob),'gx','DisplayName','min')
  hold off
  legend('Location','southeast')
  xlim([1,nSes])
  ylim([0.5 1])
  xlabel('Session #')
  ylabel('average fp correlation')
  
  subplot(3,2,6)
  hold on
  plot(nanvar(score.prob),'kx','DisplayName','p_{same}')
  plot(nanvar(score.fp_corr_oneway(:,:,1)),'rD','DisplayName','1-way fp')
  plot(nanvar(score.fp_corr_oneway(:,:,2)),'rs','DisplayName','1-way fp')
  hold off
  legend()
  xlim([1,nSes])
  ylim([0 0.1])
  xlabel('Session #')
  ylabel('variance fp correlation')
  
  
  
end