

function plot_clusters(clusters,c)
  
  figure('position',[200 200 1800 1200]);
  subplot(1,2,1)
  hold on;
  for s=1:15
    if ~isnan(clusters(c,s).neuron)
      contour(reshape(clusters(c,s).A,512,512),[0.1 0.1],'k')
    end
  end
  hold off
  xlim([clusters(c,1).centroid(2)-10,clusters(c,1).centroid(2)+10])
  ylim([clusters(c,1).centroid(1)-10,clusters(c,1).centroid(1)+10])
  
  for s = 1:15
    if ~isnan(clusters(c,s).neuron)
      subplot(15,2,s*2)
      plot(clusters(c,s).CaTrace,'k')
    end
  end
end