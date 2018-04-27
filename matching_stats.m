

function matching_stats(cluster)
  
  nSes = length(cluster(2).list);
  nCluster = length(cluster);
  
  nMatches = [cluster.ct];
  hist_matches = hist(nMatches);
  
  fig = figure('position',[100 100 900 900])
  
  subplot(2,2,2)
  hist(nMatches,linspace(1,nSes,nSes))
  xlim([0,nSes])
  xlabel('detected in # of sessions')
  ylabel('# of pairs')
  text(0.4,0.8,sprintf('# (s>=3): %d',sum(hist_matches(3:end))),'units','normalized','FontSize',10)
  
  nROI = zeros(nSes,1);
  nROI_match = zeros(nSes,1);
  nROI_range_match = zeros(nSes-1,1);
  nROI_match(1) = nan;
  
  for c = 1:nCluster
    
    if ~isempty(cluster(c).list)
    %% get number of ROIs detected in this session
    
      %% get number of ROIs detected same as last session
      for s = 1:nSes
        if cluster(c).list(s)
          nROI(s) = nROI(s)+1;
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
  
  subplot(2,2,1)
  hold on
  plot(nROI,'k-')
  plot(nROI_match,'k--')
  hold off
  xlim([1,nSes])
  limsy=get(gca,'YLim');
  set(gca,'Ylim',[0 limsy(2)]);
  ylabel('total number')
  xlabel('session')
  
  subplot(2,2,3)
  plot(nROI_match./nROI)
  xlabel('session')
  ylabel('% matched with previous')
  xlim([1,nSes])
  ylim([0,1])
  limsy=get(gca,'YLim');
  set(gca,'Ylim',[0 limsy(2)]);
  
  subplot(2,2,4)
  plot(nROI_range_match)
  xlabel('\Delta s')
  ylabel('# matched over \Delta s')
  
  basePath = '/media/mizuta/AS2/';
  mouse = 884;
  path = sprintf('%s%d/match_stats_%d.jpg',basePath,mouse,nSes);

  saveas(fig,path,'jpg')
    
  disp(sprintf('saved under %s',path))
  
  
  
  
  fig = figure;
  hold on
  hist(nMatches(nMatches<3),linspace(1,nSes,nSes),'FaceColor','y')
  hist(nMatches(nMatches>=3),linspace(1,nSes,nSes),'FaceColor','b')
%    hist(nMatches(nMatches<3),linspace(1,nSes,nSes),'y')
  hold off
  xlim([-1.5,nSes+0.5])
  ylim([0 550])
  xlabel('# sessions')
  ylabel('# ROIs')
  title(sprintf('# stable ROIs detected: %d',sum(hist_matches(3:end))))
%    text(0.4,0.8,sprintf('# (s>=3): %d',sum(hist_matches(3:end))),'units','normalized','FontSize',10)
  basePath = '/media/mizuta/AS2/';
  mouse = 884;
  path = sprintf('%s%d/match_stats2_%d.jpg',basePath,mouse,nSes);
  
%    path = sprintf('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/match_stats_%d.jpg',nSes);
  saveas(fig,path,'jpg')
    
  disp(sprintf('saved under %s',path))
end