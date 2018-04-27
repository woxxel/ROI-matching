

function matching_stats(cluster)
  
  nSes = length(cluster(2).list);
  nCluster = length(cluster);
  
  nMatches = [cluster.ct];
  hist_matches = hist(nMatches);
  
  fig = figure('position',[100 100 800 500])
  
  subplot(1,2,2)
  
  hold on
  histogram(nMatches(nMatches<3),linspace(1,nSes+1,nSes+1),'FaceColor','y')
  histogram(nMatches(nMatches>=3),linspace(1,nSes+1,nSes+1),'FaceColor','b')
%    hist(nMatches(nMatches<3),linspace(1,nSes,nSes),'y')
  hold off
  xlim([-1,nSes+1])
  ylim([0 550])
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
  
  nROI_test = 0;
  
  for c = 1:nCluster
    
    if ~isempty(cluster(c).list)
    %% get number of ROIs detected in this session
    
      %% get number of ROIs detected same as last session
      for s = 1:nSes
        if nnz(cluster(c).list(s,:)) > 0
          nROI(s,1) = nROI(s,1)+nnz(cluster(c).list(s,:));
          if nnz(cluster(c).list) == 1
            nROI(s,2) = nROI(s,2) + 1;
          end
          
          if s > 1
            if cluster(c).list(s-1)
              nROI_match(s) = nROI_match(s)+1;
            end
            
            for a = 1:s-1
              if cluster(c).list(s-a)
                nROI_range_match(a) = nROI_range_match(a) + 1;
              end
            end
          end
        end
      end
    end
  end
  
  disp(sprintf('total number of ROIs: %d, %d', nROI_test, sum(nROI(:,1))))
  
  subplot(1,2,1)
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
  
  ax3 = axes('position',[0.25,0.2,0.2,0.2])
  plot(ax3,nROI_match./nROI(:,1))
  xlabel(ax3,'session')
  ylabel(ax3,'% matched')
  xlim(ax3,[1,nSes])
  ylim(ax3,[0,1])
  limsy=get(ax3,'YLim');
  set(ax3,'Ylim',[0 limsy(2)]);
  
%    subplot(2,2,4)
%    plot(nROI_range_match)
%    xlabel('\Delta s')
%    ylabel('# matched over \Delta s')
  
  basePath = '/media/wollex/AS2/';
  mouse = 884;
  path = sprintf('%s%d/match_stats_%d.jpg',basePath,mouse,nSes);

  saveas(fig,path,'jpg')
    
  disp(sprintf('saved under %s',path))
  
  
  
  
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