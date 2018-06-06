


function [ROI_cluster_final,ROI_data] = cell_matching(basePath,mouse,p_thr)
    
    disp('loading data from file...')
    tic
    pathData = sprintf('%s%d/cluster_registered.mat',basePath,mouse);
    data = load(pathData);
    ROI_cluster = data.ROI_cluster;
    ROI_data = data.data;
    ROI_xdata = data.xdata;
    disp('...done')
    toc
    
    nSes = size(ROI_cluster(1).list,1);
    
    mode = 'threshold';
%      mode = 'other';
    
    
    %% now, go through all clusters and assign surely matching ROIs to each other (p_same>0.95)
    %%% here, implementing footprints in the matching process should help/improve the results quite a bit
    
    %% afterwards, check chance of others belonging to the same cluster or whether chance is larger of them to form an own cluster
    %% for ROIs in same session, check whether merging improves matching probability
    %% also, remove surely matched ROIs in one cluster from others (or rather, track, which ones are matched already
    tic
    merge_ct = 0;
    merge2_ct = 0;
    merge_ct_real = 0;
    c = 1;
    c_final = 0;
    change_ct = 0;
    switch_ct = 0;
    A_thr = 400;
    
    disp('registering')
    while c < length(ROI_cluster)
      
      if mod(c,100)==0
        disp(c)
      end
      %% first, remove matched ROIs from ROI_cluster
      [s w] = find(ROI_cluster(c).list);
      for i = 1:length(s)
        n = ROI_cluster(c).list(s(i),w(i));
        if ROI_data(s(i)).matched(n)
          ROI_cluster(c).list(s(i),w(i)) = 0;
        end
      end
      
      if nnz(ROI_cluster(c).list) < 1   %% only look at ROI_clusters, that actually have some matching possibilities
        ROI_cluster(c) = [];
      else
        c_final = c_final + 1;
        n_ref = 0;
        s_ref = 0;
        
        %% merge status: for every neuron in the final_list, have 3 entries: previous, current and following session match status
        %% match status does not refer to matching to a certain neuron, but rather assigning to this ROI_cluster!
        ROI_cluster_final(c_final) = struct('list',zeros(nSes,1),'merge_status',false(nSes,3),'ct',0,'score',0);
        
        for s = 1:nSes
        
          if any(ROI_cluster(c).list(s,:))
            
            %% compare to last registered neuron (closest in time)
            %% also, compare to other ones if no fit found (or to overall ROI_cluster?)
            
            if n_ref == 0   %% register new ROI as reference ROI
              %%% missing here: no merging in first session possible
              n = ROI_cluster(c).list(s,ROI_cluster(c).list(s,:)>0);
              n = n(1);
              ROI_cluster_final(c_final).list(s,1) = n;
              
              %% set reference to first ROI detected
              n_ref = n;
              s_ref = s;
              
%                %% have alternative reference to most recently detected ROI ->
%                %% always check this one as well! (merging possibilities etc should be detected better like that)
%                n_ref_alt = n;
%                s_ref_alt = s;
            else
              
              if strcmp(mode,'threshold')
                [matches_s, p_same_s] = get_matches(ROI_cluster(c),ROI_xdata,0.05,s_ref,n_ref,s);
                [p_best_s,idx_s] = max(p_same_s);
                if p_best_s > p_thr
                  best_match_s = matches_s(idx_s);
                  
                  %% check for reciprocity
                  [matches_s_ref, p_same_s_ref] = get_matches(ROI_cluster(c),ROI_xdata,0.05,s,best_match_s,s_ref);
                  [p_best_s_ref,idx_s_ref] = max(p_same_s_ref);
                  if (matches_s_ref(idx_s_ref) == n_ref) && (p_best_s_ref > p_thr)
                    ROI_cluster_final(c_final).list(s,1) = best_match_s;
                  end
                end
              
              
              %% matching due to most probable ROI (including merging etc)
              else
                
                %% check for matches with first detected ROI
                %% also, check for matches with most recently detected ROI
                [matches_s, p_same_s] = get_matches(ROI_cluster(c),ROI_xdata,0.05,s_ref,n_ref,s);
                [~,idx_s] = max(p_same_s);
                best_match_s = matches_s(idx_s);
                
                for i=1:length(matches_s)
                %%% should only first best match (matches_s_ref) be considered? or also 2nd best?
                  [matches_s_ref, p_same_s_ref] = get_matches(ROI_cluster(c),ROI_xdata,0.05,s,matches_s(i),s_ref);
                  [~,idx_s_ref] = max(p_same_s_ref);
                  
                  if matches_s(i) == best_match_s && matches_s_ref(idx_s_ref) == n_ref    %% if they are each others favorites
                    %% additionally check, whether this probability is larger than ... something?!
                    if p_same_s_ref(idx_s_ref) > 0.05
                      ROI_cluster_final(c_final).merge_status(s_ref,3) = true;
                      ROI_cluster_final(c_final).merge_status(s,1) = true;
                      
                      if ~ismember(matches_s(i),ROI_cluster_final(c_final).list(s,:))
                        idx = nnz(ROI_cluster_final(c_final).list(s,:)) + 1;
                        ROI_cluster_final(c_final).list(s,idx) = matches_s(i);
                      end
                    end
                    
                  elseif matches_s(i) == best_match_s                 %% if chosen ROI rather matches with another one
                  %% do not match!! (or rather: how much different are they? look for merging possibility?)
                  %%% here, should check for 2nd best match
                    if (p_same_s_ref(idx_s_ref) - p_same_s(idx_s) > 0.5)  %% really wants to match another one -> do not include in this ROI_cluster (very rare)
                      change_ct = change_ct + 1;
                    else                                           %% if there might be a chance of both matching -> merge?
                      ROI_cluster_final(c_final).merge_status(s_ref,2:3) = true;
                      ROI_cluster_final(c_final).merge_status(s,1) = true;
                      merge_ct = merge_ct + 1;
                      if ~ismember(matches_s_ref(idx_s_ref),ROI_cluster_final(c_final).list(s_ref,:))
                        idx = nnz(ROI_cluster_final(c_final).list(s_ref,:)) + 1;
                        ROI_cluster_final(c_final).list(s_ref,idx) = matches_s_ref(idx_s_ref);
                      end
                      if ~ismember(matches_s(i),ROI_cluster_final(c_final).list(s,:))
                        idx = nnz(ROI_cluster_final(c_final).list(s,:)) + 1;
                        ROI_cluster_final(c_final).list(s,idx) = matches_s(i);
                      end
                    end
                    
                  elseif matches_s_ref(idx_s_ref) == n_ref
                    if (p_same_s_ref(idx_s_ref) - p_same_s(idx_s) > 0.5)  %% if probabilities far exceed, change match
                      ROI_cluster_final(c_final).merge_status(s_ref,3) = true;
%                        ROI_cluster_final(c_final).list(s,:) = 0;
                      if ~ismember(matches_s(i),ROI_cluster_final(c_final).list(s,:))
                        idx = nnz(ROI_cluster_final(c_final).list(s,:)) + 1;
                        ROI_cluster_final(c_final).list(s,idx) = matches_s(i);
                        ROI_cluster_final(c_final).merge_status(s,1) = true;
                      end
                      switch_ct = switch_ct + 1;
                    else
                      
                      merge2_ct = merge2_ct + 1;
                      ROI_cluster_final(c_final).merge_status(s,2) = true;
                      if ~ismember(matches_s(i),ROI_cluster_final(c_final).list(s,:))
                        idx = nnz(ROI_cluster_final(c_final).list(s,:)) + 1;
                        ROI_cluster_final(c_final).list(s,idx) = matches_s(i);
                      end
                    end
                    best_match_s = matches_s(i);
                  end
                  
                end
                
                if any(ROI_cluster_final(c_final).list(s,:))
                  %% allow more than one neuron to go here
                  n_ref_alt = best_match_s;   %% this should include merging possibilities
                  s_ref_alt = s;
                end
              end
            end
          end
        end
        
        %% obtain and calculate values for ROI score
        score = prepare_ROI_score(ROI_cluster_final(c_final),ROI_data,ROI_xdata);
        
        ROI_cluster_final(c_final) = ROI_cleanup(ROI_cluster_final(c_final),score,c_final,ROI_data);
        ROI_cluster_final(c_final).score = get_ROI_score(ROI_cluster_final(c_final),score,0);
        
        %%% filling gaps in ROI cluster should be done in other function for all ROIs > certain score and ct
        %%% it should...
        %%% 1. crop out region that encloses all ROIs from this cluster + some margin
        %%% 2. find all closeby ROIs (also from other clusters)
        %%%   2.1. check if there is a single ROI, that might have been sorted out but belongs to this cluster- if so, remove!
        %%% 3. initiate CNMF with initial guess of closeby ROIs + region covered by this cluster (+ some margin)
        %%% 4. if new ROI is found, implement this one + its Ca-trace in data
        %%%   4.1 if no new ROI is found, remark this one as "non-active" (or just apply average ROI from neighbouring sessions and get Ca-trace from simple filter application (- background) - check, if active or not)
        
        
        %% here, implement checking for ROI score and removing/merging/splitting accordingly
        
        %% score: high average and minimum probability, bias towards large number of neurons in one ROI_cluster
        %% check: removing one ROI from ROI_cluster: does it increase or decrease the "score"?
        %% or: possible to create subset from ROI_cluster that has high average and minimum probability?
        
        %% only now, after removing "substandard matches" from the ROI_cluster, assign "matched" status to all
          
          %% assign matched status to all within ROI_cluster
          [s w] = find(ROI_cluster_final(c_final).list);
          for i = 1:length(s)
            n = ROI_cluster_final(c_final).list(s(i),w(i));
            ROI_data(s(i)).matched(n) = true;
            ROI_data(s(i)).cluster_neuron(c) = n;     %%% assign cluster-id to ROI
          end
          
          if nnz(ROI_cluster_final(c_final).list) < 1
            ROI_cluster_final(c_final) = [];
            c_final = c_final - 1;
          end
        
          ROI_cluster_final(c_final).ct = nnz(ROI_cluster_final(c_final).list(:,1));
          
          %% in the end, create new ROI_cluster from remaining ROIs and append to ROI_cluster struct
          if nnz(ROI_cluster_final(c_final).list) ~= nnz(ROI_cluster(c).list)
            
            c_idx = length(ROI_cluster)+1;
            for s = 1:nSes
              n = ROI_cluster(c).list(s,ROI_cluster(c).list(s,:)>0);
              n = n(~ROI_data(s).matched(n));
              ROI_cluster(c_idx).list(s,1:length(n)) = n;
            end
            if nnz(ROI_cluster(c_idx).list) < 1
  %              disp('purging needed?')
              ROI_cluster(c_idx) = [];
            end
          end
%          end
        
        c = c+1;
      end
    end
    
    %%% fill up cluster_neuron arrays to cover all clusters
    for s = 1:nSes
      if length(ROI_data(s).cluster_neuron) < c_final
        ROI_data(s).cluster_neuron(c_final) = 0;
      end
%        [s length(ROI_data(s).cluster_neuron)]
%        ROI_data(s).cluster_neuron = cat(1,ROI_data(s).cluster_neuron',zeros(c_final - length(ROI_data(s).cluster_neuron),1))
      [s size(ROI_data(s).cluster_neuron)]
    end
    
    nMatches = [ROI_cluster_final.ct];
    disp(sprintf('number of ROI_clusters: %d',c_final))
    disp(sprintf('merging attempts: %d',merge_ct))
    disp(sprintf('real merges to be done: %d',merge_ct_real))
%      disp(sprintf('number of session-matchings: %d',sesmatch))
%      disp(sprintf('polygamous ROIs: %d',polygamy))
%      disp('matching done')
    toc
    fig_ses = figure('position',[100 100 800 400]);
    histogram(nMatches)
    xlabel('# sessions detected')
    ylabel('# matched ROIs')
    
    hist_matches = hist(nMatches);
    text(0.6,0.9,sprintf('# stable ROIs (s>=3): %d',sum(hist_matches(3:end))),'units','normalized','FontSize',14)
    
    %      basePath = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data';
    
%      path = sprintf('%s%d/hist_s=%d.jpg',bassePath,mouse,nSes);
    
%      path = sprintf('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/hist_s=%d.jpg',nSes);
%      saveas(fig_ses,path,'jpg')
    
    savePath = sprintf('%s%d/ROI_final_p=%4.2g.mat',basePath,mouse,p_thr);    
    save(savePath,'ROI_cluster_final','ROI_data','-v7.3')
    disp(sprintf('saved under %s',savePath))
    
end



function [n, p_same] = get_matches(ROI_cluster,ROI_xdata,p_thr,s_ref,n_ref,s)
  
  %% search for all ROIs, that are not certainly rejected due to footprint or distance
  n = ROI_cluster.list(s,ROI_cluster.list(s,:)>0);
  
  model_type = 'joint';
  switch model_type
    case 'dist'
      p_same = full(ROI_xdata(s_ref,s).p_same_dist(n_ref,n));
    case 'corr'
      p_same = full(ROI_xdata(s_ref,s).p_same_corr(n_ref,n));
    case 'joint'
      p_same = full(ROI_xdata(s_ref,s).p_same_joint(n_ref,n));
  end
  
  mask = p_same>p_thr;
  n = n(mask);
  p_same = p_same(mask);
  
end


function [ROI_cluster] = ROI_cleanup(ROI_cluster,score,c_final,data)
  
  cleanup = true;
  
  while cleanup
    ROI_score_old = get_ROI_score(ROI_cluster,score,0);
    
    % get average probability and check all ROIs that are below
    s_test = find(nanmean(score.prob) < nanmean(score.prob(:))-nanvar(score.prob(:)));
    ROI_score_test = zeros(length(s_test),1);
    
    %% here, should be checked for merging, splitting, removing some or removing all
    for i = 1:length(s_test)
      s = s_test(i);
      ROI_score_test(i) = get_ROI_score(ROI_cluster,score,s);
    end
    [ROI_score_new,i_rm] = max(ROI_score_test);
    s_rm = s_test(i_rm);
    
    if (ROI_score_new - ROI_score_old) > 0.05
      n_rm = ROI_cluster.list(s_rm,1);
      
      if nnz(ROI_cluster.list(:)) > 10
        plot_ROI_cluster(ROI_cluster,data,s_rm,score)
        pause(0.1)
      end
      score.prob(s_rm,:) = NaN;
      score.prob(:,s_rm) = NaN;
      score.fp_corr_oneway(s_rm,:,:) = NaN;
      score.fp_corr_oneway(:,s_rm,:) = NaN;
      ROI_cluster.list(s_rm,:) = 0;
      
      disp(sprintf('removing ROI from cluster %d in session #%d from ROI_cluster to enhance ROI_cluster score from %6.4g to %6.4g \t dp = %6.4g',c_final,s_rm,ROI_score_old,ROI_score_new,(ROI_score_new-ROI_score_old)))
%              pause(1)
    else
      cleanup=false;
    end
    
  end
end