

function [nROI, nROI_match] = matching_stats(cluster)
  
  nSes = length(cluster(2).list);
  nCluster = length(cluster);
  
  nSes
  nCluster
  
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
  
  figure
  subplot(2,1,1)
  hold on
  plot(nROI,'k-')
  plot(nROI_match,'k--')
  hold off
  xlim([1,nSes])
  limsy=get(gca,'YLim');
  set(gca,'Ylim',[0 limsy(2)]);
  ylabel('total number')
  xlabel('session')
  
  subplot(2,1,2)
  plot(nROI_match./nROI)
  xlabel('session')
  ylabel('% matched')
  xlim([1,nSes])
  limsy=get(gca,'YLim');
  set(gca,'Ylim',[0 limsy(2)]);
  figure
  plot(nROI_range_match)
end