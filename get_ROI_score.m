

function [ROI_score] = get_ROI_score(ROI_cluster,score,mode,s_rm)
  
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
      score.fp_corr_oneway(s,:,:) = NaN;
      score.fp_corr_oneway(:,s,:) = NaN;
    end
    N=N-length(s_rm);
  end
  
  p_mean = nanmean(score.prob(:));
  p_var = nanvar(score.prob(:));
  
  fp_mean = squeeze(nanmean(score.fp_corr_oneway,1));
  fp_max = nanmax(fp_mean,[],2);
  
  frac_1w_mean = nanmean(fp_max);
  frac_1w_var = nanvar(fp_max);
  
  if strcmp(mode,'compare')
    p = 0.15; %% probability-equal of adding or subtracting one ROI (should be ~5% per summand)
    baseline = 1/p;
    ROI_delete = ((baseline-nSes)/baseline + N/baseline);
    
  %    [p_mean frac_1w_mean ROI_delete]
  %    [p_mean^(1+p_var) frac_1w_mean^(1+frac_1w_var) ROI_delete]  %nanmean(nanmin(score.prob)) 
    
  %    ROI_score = 1/3*(p_mean^(1+p_var) + frac_1w_mean^(1+frac_1w_var) + nanmin(score.prob(:)))*N/nSes;
  %    ROI_score = (p_mean^(1+p_var) * frac_1w_mean^(1+frac_1w_var) * nanmin(score.prob(:)))^(1/3)*N/nSes;
    
    ROI_score = 1/4*(p_mean^(1+p_var) + 2*frac_1w_mean^(1+frac_1w_var) + ROI_delete);
  
%    disp(sprintf('s: %d ROI score: %6.4g',s_rm(1), ROI_score))
  elseif strcmp(mode,'final')
    ROI_score = 1/3*(p_mean^(1+p_var) + 2*frac_1w_mean^(1+frac_1w_var));
  end
    

end
