

function [score] = prepare_ROI_score(ROI_cluster,data,xdata)
  
  nSes = size(ROI_cluster.list,1);
  width= size(ROI_cluster.list,2);
  N = nnz(sum(ROI_cluster.list,2));
  
  [y_idx x_idx] = find(ROI_cluster.list);
  entries = length(y_idx);
  
  score = struct;
  
  score.prob = zeros(nSes,nSes);
  score.fp_corr_oneway = zeros(nSes,nSes,2);
  
  score.prob(:) = nan;
  score.fp_corr_oneway(:) = NaN;
  
  for i = 1:entries;
    s = y_idx(i);
    l = x_idx(i);
    n = ROI_cluster.list(s,l);
    for w = 1:width
      for sm=1:nSes
        m = ROI_cluster.list(sm,w);
        if m && sm~=s
          score.prob(s,sm) = xdata(s,sm).p_same_joint(n,m);
          
          score.fp_corr_oneway(s,sm,1) = get_1w_corr(data,[s,n],[sm,m]);
          score.fp_corr_oneway(s,sm,2) = get_1w_corr(data,[sm,m],[s,n]);
        end
      end
    end
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