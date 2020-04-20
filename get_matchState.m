

function matchState = get_matchState(data,clusters)
  
  
  matchState = struct('session',struct('neuron',cell(data.nSes,1)),'clusters',struct('list',zeros(data.nSes,1)));
  
  if length(clusters) > 1
    nC = length(clusters);
  else
    nC = clusters
  end
  
  for c = 1:nC
    matchState.clusters(c).list = zeros(data.nSes,1);
  end
  
  for s = 1:data.nSes
    matchState.session(s).neuron = struct('cluster_ID',cell(data.session(s).nROI,1));
    
    for n = 1:data.session(s).nROI
      
      matchState.session(s).neuron(n).cluster_ID = data.session(s).ROI(n).cluster_ID;
      c_arr = matchState.session(s).neuron(n).cluster_ID;
      for c = c_arr
        matchState.clusters(c).list(s,1) = n;
      end
    end
  end
  
end