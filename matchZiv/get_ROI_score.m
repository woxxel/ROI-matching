

function [ROI_score] = get_ROI_score(ROI_cluster,score,s_rm)
  
  %%% could also take into account:
  %%%   - correlation with direct neighbours
  %%%   - movement to consistent direction
  %%%   - minimum correlation
  %%%   - what else?
  
  nSes = size(ROI_cluster.list,1);
  N = nnz(sum(ROI_cluster.list,2));
  
  if s_rm>0
    for s=s_rm
      score.prob(s,:) = NaN;
      score.prob(:,s) = NaN;
%        score.fp_corr_oneway(s,:,:) = NaN;
%        score.fp_corr_oneway(:,s,:) = NaN;
    end
    N=N-length(s_rm);
  end
  
  p_mean = nanmean(score.prob(:));
  p_var = nanvar(score.prob(:));
  
%    fp_mean = squeeze(nanmean(score.fp_corr_oneway,1));
%    fp_max = nanmax(fp_mean,[],2);
  
%    frac_1w_mean = nanmean(fp_max);
%    frac_1w_var = nanvar(fp_max);
    
  %    [p_mean frac_1w_mean]
  %    [p_mean^(1+p_var) frac_1w_mean^(1+frac_1w_var)]  %nanmean(nanmin(score.prob)) 
    
  %    ROI_score = 1/3*(p_mean^(1+p_var) + frac_1w_mean^(1+frac_1w_var) + nanmin(score.prob(:)));
  %    ROI_score = (p_mean^(1+p_var) * frac_1w_mean^(1+frac_1w_var) * nanmin(score.prob(:)))^(1/3);
  
%    ROI_score = 1/3*(p_mean^(1+p_var) + 2*frac_1w_mean^(1+frac_1w_var));
  ROI_score = p_mean^(1+p_var);

end
