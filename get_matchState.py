

def get_matchState(data,clusters)
  
  
  matchState = struct('session',struct('neuron',cell(data.nSes,1)),'clusters',struct('list',zeros(data.nSes,1)));
  
  if len(clusters) > 1:
    nC = len(clusters)
  else:
    nC = clusters
  
  for c in range(nC):
    matchState.clusters(c).list = zeros(data.nSes,1)
  
  for s in range(data.nSes):
    matchState.session(s).neuron = struct('cluster_ID',cell(data.session(s).nROI,1));
    
    for n in range(data.session(s).nROI):
      
      matchState.session(s).neuron(n).cluster_ID = data.session(s).ROI(n).cluster_ID;
      c_arr = matchState.session(s).neuron(n).cluster_ID;
      for c in c_arr:
        matchState.clusters(c).list(s,1) = n
  
  return matchState