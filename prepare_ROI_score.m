

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
  
  mask_1w(entries) = struct('neurons',struct);
  for i = 1:entries;
    s = y_idx(i);
    l = x_idx(i);
    n = ROI_cluster.list(s,l);
    mask_1w(s).neurons(n).idx = find(data(s).A(:,n));
  end
    
  for i = 1:entries;
    s = y_idx(i);
    l = x_idx(i);
    n = ROI_cluster.list(s,l);
    for j = i+1:entries
      sm = y_idx(j);
      l = x_idx(j);
      m = ROI_cluster.list(sm,l);
      if sm~=s
        val_tmp = xdata(s,sm).p_same_joint(n,m);
        score.prob(s,sm) = val_tmp;
        score.prob(sm,s) = val_tmp;
        
        val_tmp = get_1w_corr(data,mask_1w(s).neurons(n).idx,[s,n],[sm,m]);
        score.fp_corr_oneway(s,sm,1) = val_tmp;
        score.fp_corr_oneway(sm,s,2) = val_tmp;
        
        val_tmp = get_1w_corr(data,mask_1w(sm).neurons(m).idx,[sm,m],[s,n]);
        score.fp_corr_oneway(s,sm,2) = val_tmp;
        score.fp_corr_oneway(sm,s,1) = val_tmp;
      end
    end
  end
end





function [corr_1w] = get_1w_corr(data,A_idx,N,M)
  corr_1w = full(dot(data(N(1)).A(A_idx,N(2)),data(M(1)).A(A_idx,M(2)))/(data(N(1)).norm(N(2))*norm(data(M(1)).A(A_idx,M(2)))));
end