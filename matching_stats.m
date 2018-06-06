

function matching_stats(ROI_cluster,data)
  
  nSes = length(ROI_cluster(2).list);
  nCluster = length(ROI_cluster);
  
  nMatches = [ROI_cluster.ct];
  hist_matches = hist(nMatches);
  
  fig = figure('position',[100 100 1400 600]);
  
  subplot(1,3,2)
  
  hold on
  histogram(nMatches(nMatches<3),linspace(1,nSes+1,nSes+1),'FaceColor','y')
  histogram(nMatches(nMatches>=3),linspace(1,nSes+1,nSes+1),'FaceColor','b')
%    hist(nMatches(nMatches<3),linspace(1,nSes,nSes),'y')
  hold off
  xlim([-1,nSes+1])
  ylim([0 300])
  xlabel('# sessions')
  ylabel('# ROIs')
  title(sprintf('# stable ROIs detected: %d',sum(hist_matches(3:end))))
  
  
  
%    hist(nMatches,linspace(1,nSes,nSes))
%    xlim([0,nSes])
%    xlabel('detected in # of sessions')
%    ylabel('# of pairs')
%    text(0.4,0.8,sprintf('# (s>=3): %d',sum(hist_matches(3:end))),'units','normalized','FontSize',10)
%    
  nROI = zeros(nSes,2);
  nROI_match = zeros(nSes,1);
  nROI_range_match = zeros(nSes-1,1);
  nROI_match(1) = nan;
  testROI(nSes) = struct;
  for s=1:nSes
    testROI(s).neurons = zeros(data(s).nROI,1);
  end
  
  
  ROI_cluster_stats = zeros(nCluster,2);    % first entry is number of ROIs in there, 2nd is ROI score
  time_prep = 0;
  time_score = 0;
  for c = 1:nCluster
    ROI_cluster_stats(c,1) = nnz(sum(ROI_cluster(c).list,2));
    ROI_cluster_stats(c,2) = ROI_cluster(c).score;
    
    if ~isempty(ROI_cluster(c).list)
    %% get number of ROIs detected in this session
    
      %% get number of ROIs detected same as last session
      for s = 1:nSes
        if nnz(ROI_cluster(c).list(s,:)) > 0
          for n = ROI_cluster(c).list(s,:)
            testROI(s).neurons(n) = testROI(s).neurons(n) + 1;
          end
          nROI(s,1) = nROI(s,1)+nnz(ROI_cluster(c).list(s,:));
          if nnz(ROI_cluster(c).list) == 1
            nROI(s,2) = nROI(s,2) + 1;
          end
          
          if s > 1
            if ROI_cluster(c).list(s-1)
              nROI_match(s) = nROI_match(s)+1;
            end
            
            for a = 1:s-1
              if ROI_cluster(c).list(s-a)
                nROI_range_match(a) = nROI_range_match(a) + 1;
              end
            end
          end
        end
      end
    end
  end
  
  for s=1:nSes
    disp(sprintf('session: %d:  mean: %6.4g, max: %6.4g',s,mean(testROI(s).neurons),max(testROI(s).neurons)))
  end
  
  
  disp(sprintf('total number of ROIs: %d', sum(nROI(:,1))))
  
  subplot(1,3,1)
  hold on
  plot(nROI(:,1),'k-')
  plot(nROI(:,1)-nROI(:,2),'r-')
  plot(nROI_match,'k--')
  hold off
  xlim([1,nSes])
  limsy=get(gca,'YLim');
  set(gca,'Ylim',[0 limsy(2)]);
  ylabel('total number')
  xlabel('session')
  
  ax3 = axes('position',[0.17,0.2,0.15,0.15]);
  plot(ax3,nROI_match./nROI(:,1))
  xlabel(ax3,'session')
  ylabel(ax3,'% matched')
  xlim(ax3,[1,nSes])
  ylim(ax3,[0.3,1])
  
%    subplot(2,2,4)
%    plot(nROI_range_match)
%    xlabel('\Delta s')
%    ylabel('# matched over \Delta s')
  
%    basePath = '/media/wollex/AS2/';
%    mouse = 884;
%    path = sprintf('%s%d/match_stats_%d.jpg',basePath,mouse,nSes);

%    saveas(fig,path,'jpg')
    
%    disp(sprintf('saved under %s',path))
  
  
  counts = [ROI_cluster.ct];
  scores = [ROI_cluster.score];
  avg_score = zeros(nSes,1);
  
  for ct = 1:nSes
    avg_score(ct) = nanmean(scores(counts==ct));
  end
  
  subplot(2,3,6)
  hold on
  scatter(counts,scores,'kx')
  plot(1:nSes,avg_score,'r-','LineWidth',2)
  hold off
  xlabel('ROI count')
  ylabel('cluster score')
  
  subplot(2,3,3)
  histogram(scores,'FaceColor','r')
  xlabel('cluster score')
  xlim([0.5,1])
  
  
  
  idx_good = find(counts>=floor(nSes/3) & scores >=0.9);
  gaps = zeros(nSes,length(idx_good));
  
  sort_mat = [counts(idx_good)',scores(idx_good)',idx_good'];
  idx_sorted = sortrows(sort_mat);
  for i = 1:length(idx_good)
    gaps(:,i) = ROI_cluster(idx_sorted(i,3)).list==0;
  end
  
  figure
  imagesc(gaps)
  colormap(gca,'gray')
  xlabel('cluster')
  ylabel('session #')
  title('missing ROIs')
  
  
%    fig = figure;
%    hold on
%    histogram(nMatches(nMatches<3),linspace(1,nSes+1,nSes+1),'FaceColor','y')
%    histogram(nMatches(nMatches>=3),linspace(1,nSes+1,nSes+1),'FaceColor','b')
%  %    hist(nMatches(nMatches<3),linspace(1,nSes,nSes),'y')
%    hold off
%    xlim([-1.5,nSes+0.5])
%    ylim([0 550])
%    xlabel('# sessions')
%    ylabel('# ROIs')
%    title(sprintf('# stable ROIs detected: %d',sum(hist_matches(3:end))))
%  %    text(0.4,0.8,sprintf('# (s>=3): %d',sum(hist_matches(3:end))),'units','normalized','FontSize',10)
%    basePath = '/media/wollex/AS2/';
%    mouse = 884;
%    path = sprintf('%s%d/match_stats2_%d.jpg',basePath,mouse,nSes);
%    
%  %    path = sprintf('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/match_stats_%d.jpg',nSes);
%    saveas(fig,path,'jpg')
%      
%    disp(sprintf('saved under %s',path))
end