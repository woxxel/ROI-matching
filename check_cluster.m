%% check validity of ROI detection algorithm

function check_ROI_cluster(ROI_cluster,data,xdata,para,p_thr,idx)
  
  close all
  %% pick a random ROI_cluster,if none given
  if nargin < 6 || isempty(idx)
    idx = randi(length(ROI_cluster))
  end
  disp(sprintf('displaying ROI_cluster #%d',idx))
  
  ROI_cluster = ROI_cluster(idx);
  ROI_cluster.list
%    ROI_cluster(idx).list = ROI_cluster(idx).list(:,1);
%    ROI_cluster(idx).list
  
  %% if more than one entry:
  if length(find(ROI_cluster.list)) > 1
    
    %% initialize:    
    nSes = size(ROI_cluster.list,1);
    
    score = prepare_ROI_score(ROI_cluster,data,xdata);
    ROI_cluster.score
%      cleanup = true;
    
    plot_ROI_cluster(ROI_cluster,score,data,0);
    
%      while cleanup
%        ROI_score_old = get_ROI_score(ROI_cluster,score,0);
%        disp(sprintf('old ROI score: %6.4g',ROI_score_old))
%        
%        % get average probability and check all ROIs that are below
%        disp('find ROIs that can be removed')
%        s_test = find(nanmean(score.prob) < nanmean(score.prob(:))-nanvar(score.prob(:)));
%        
%        ROI_score_test = zeros(length(s_test),1);
%        for i = 1:length(s_test)
%          s = s_test(i);
%  %          if nnz(ROI_cluster.list(s,:)) > 1  %% something's wrong with merge status
%  %            disp('test for merging...')
%  %            %% should test effect of merging vs removing extra ROIs (or all)
%  %            ROI_score_new = get_ROI_score_merge(ROI_cluster,data,model,para,res,s)   %% merge data.A from both ROIs, recalculate prob and fp_1way and calculate score (by calling get_ROI_score with temporary res struct)
%  %            disp(sprintf('new ROI score: %6.4g vs %6.4g',ROI_score_new,ROI_score_old))
%  %            if ROI_score_new > ROI_score_old
%  %              disp(sprintf('merging ROI %d & %d in session #%d from ROI_cluster to enhance ROI_cluster score from %6.4g to %6.4g',ROI_cluster.list(s,1),ROI_cluster.list(s,2),s,ROI_score_old,ROI_score_new))
%  %              plot_ROIs(ROI_cluster,data,res,s)
%  %            end
%  %            
%  %          else
%            ROI_score_test(i) = get_ROI_score(ROI_cluster,score,s);
%  %          end
%        end
%        [ROI_score_new,i_rm] = max(ROI_score_test);
%        s_rm = s_test(i_rm);
%        disp(sprintf('new ROI score: %6.4g (when removing ROI from session %d)',ROI_score_new,s_rm))
%        
%        
%        if ROI_score_new > ROI_score_old
%          
%          n_rm = ROI_cluster.list(s_rm,1);
%  %          plot_check_ROI_cluster(ROI_cluster,res,data,nSes,s_rm);
%          score.prob(s_rm,:) = NaN;
%          score.prob(:,s_rm) = NaN;
%          score.fp_corr_oneway(s_rm,:,:) = NaN;
%          score.fp_corr_oneway(:,s_rm,:) = NaN;
%          ROI_cluster.list(s_rm,:) = 0;
%          
%          disp(sprintf('removing ROI in session #%d from ROI_cluster to enhance ROI_cluster score from %4.2g to %4.2g',s_rm,ROI_score_old,ROI_score_new))
%        else
%          cleanup=false;
%        end
%        
%      end
%      plot_ROI_cluster(ROI_cluster,score,data,0);
%    end
  
end



%  function [ROI_score] = get_ROI_score_merge(ROI_cluster,data,model,para,res,s_merge)
%    
%    microns_per_pixel = 530.684/512;
%    
%    imSize = [512,512];
%    nSes = size(res.prob,1);
%    %% merging data.A together (all combinations, if >2 exist?
%    A_merge = sparse(size(data(s_merge).A(:,1),1),1);
%    
%    for n = ROI_cluster.list(s_merge,:)
%      if n > 0
%        A_merge = A_merge + data(s_merge).A(:,n);
%      end
%      %% obtain reference values of prob etc?
%    end
%    A_test = bwconncomp(full(reshape(A_merge,imSize(1),imSize(2)))>0,8);
%    if A_test.NumObjects > 1
%      disp('dont merge, ROIs are not connected!')
%  %      A_test.NumObjects
%      merge_status = false;
%  %      return
%    end
%    
%    A_area = nnz(A_merge);
%    
%  %    disp(sprintf('area of the merged ROI: %d',nnz(A_merge)))
%    
%    A_norm = norm(A_merge);
%    A_tmp = reshape(A_merge/sum(A_merge(:)),imSize(1),imSize(2));
%    
%    A_centroid = zeros(nSes,2);
%    for s = 1:nSes
%      %% compare to the ROI with highest matching probability (for now) - is it always #1?
%      n = ROI_cluster.list(s,1);
%      
%      A_centroid(s,:) = [sum((1:imSize(1))*A_tmp),sum(A_tmp*(1:imSize(2))')];
%      
%      %% recalculating res. structure values (prob & fp_corr_oneway)
%      dist(s) = microns_per_pixel*sqrt((A_centroid(s,1) - data(s).centroid(n,1)).^2 + (A_centroid(s,2) - data(s).centroid(n,2)).^2);
%      fp_corr(s) = full(dot(A_merge,data(s).A(:,n))/(A_norm*data(s).norm(n)));
%      
%      idx_dist = max(1,ceil(para.nbins*dist(s)/para.dist_max));
%      idx_corr = max(1,ceil(para.nbins*fp_corr(s)/para.corr_max));
%      
%      res.prob(s,s_merge) = model.p_same_joint(idx_dist,idx_corr);
%      res.prob(s_merge,s) = model.p_same_joint(idx_dist,idx_corr);
%      
%      %% now add 1w correlation
%    end
%    
%  %    [fp_corr',res.corr(s_merge,:)']
%    ROI_score = get_ROI_score(res,s_merge);
%  
%  end
