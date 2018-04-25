

function [data,xdata,cluster,cluster_final] = cell_registration(data,xdata,model,para,histo,p_thr,nSes)

  %%% -------------------------------------- Cell registration ---------------------------------------
    
    %% implements also check for one-sided spatial correlation
    %% include obtaining an overall spatial filter and an overall centroid (from almost surely matched cells), and match other cells according to that (or get maximum probability, given single ROIs and overall ROI)
    %% include, that cell shape and position can change over time, therefore give higher significance to cells that can be matched from neighbouring sessions!
    
    microns_per_pixel = 530.684/512;
    
    model.type = 'joint';
    
    tic
%      model = 'centroids';
%      model = 'centroids';
    
    nCluster = 0;
    sesmatch = 0;
    polygamy = 0;
    cluster = struct;
    for s = 1:nSes
      data(s).registered = false(data(s).nROI,1);
      data(s).matched = false(data(s).nROI,1);
      for sm = 1:nSes;
        xdata(s,sm).p_same_dist = sparse(data(s).nROI,data(sm).nROI);
        switch model.type
%            case 'dist'
          case 'corr'
            xdata(s,sm).p_same_corr = sparse(data(s).nROI,data(sm).nROI);
          case 'joint'
            xdata(s,sm).p_same_joint = sparse(data(s).nROI,data(sm).nROI);
        end
      end
    end
    
    %% here, all cells that are not surely different are put together in clusters
    for s = 1:nSes
      
      for sm = 1:nSes
        if s == sm
          continue
        end
        
        for n = 1:data(s).nROI
          
          if ~data(s).registered(n)   %% add new cluster if not already belonging to one
            nCluster = nCluster + 1;
            data(s).cluster(n).ID = nCluster;
            data(s).registered(n) = true;
            cluster(data(s).cluster(n).ID).list = zeros(nSes,1);
            cluster(data(s).cluster(n).ID).list(s) = n;
            cluster(data(s).cluster(n).ID).ct = 1;
          end
          
          neighbours = xdata(s,sm).neighbours(n,:) > 0;
          
          idx_dist = max(1,ceil(para.nbins*xdata(s,sm).dist(n,neighbours)/para.dist_max));
%            [~,idx_dist] = min(abs(xdata(s,sm).dist(n,neighbours) - histo.dist_x));
          xdata(s,sm).p_same_dist(n,neighbours) = model.p_same_dist(idx_dist);
          xdata(sm,s).p_same_dist(neighbours,n) = model.p_same_dist(idx_dist);
          
          idx_corr = max(1,ceil(para.nbins*xdata(s,sm).fp_corr(n,neighbours)/para.corr_max));
          if strcmp(model.type,'corr')
            xdata(s,sm).p_same_corr(n,neighbours) = model.p_same_corr(idx_corr);
            xdata(sm,s).p_same_corr(neighbours,n) = model.p_same_corr(idx_corr);
          end
          
          idx_nb = find(neighbours);
          for i=1:length(idx_dist)
            xdata(s,sm).p_same_joint(n,idx_nb(i)) = model.p_same_joint(idx_dist(i),idx_corr(i));
            xdata(sm,s).p_same_joint(idx_nb(i),n) = model.p_same_joint(idx_dist(i),idx_corr(i));
          end
          
          cluster_candidates = find(xdata(s,sm).p_same_dist(n,:)>(1-p_thr));    %% all ROIs in sm that are candidates to be same as ROI (s,n)
          
          for m = cluster_candidates
            if ~data(sm).registered(m)
              data(sm).cluster(m).ID = data(s).cluster(n).ID;
              data(sm).registered(m) = true;
            end
            
            for c = data(s).cluster(n).ID         %% in case, there are several from one session that fall into one cluster
              if ~ismember(m,cluster(c).list(sm,:))
                idx = nnz(cluster(c).list(sm,:))+1;
                cluster(c).list(sm,idx) = m;
                if idx == 1
                  cluster(c).ct = cluster(c).ct + 1;
                end
              end
            end
          end
        end
      end
    end
    
    nMatches = [cluster.ct];
    disp(sprintf('number of clusters: %d',nCluster))
    disp(sprintf('number of real clusters: %d',sum(nMatches > 1)))
%      disp(sprintf('number of session-matchings: %d',sesmatch))
%      disp(sprintf('polygamous ROIs: %d',polygamy))
    toc
    disp('pre-clustering done')
    
    figure
    histogram(nMatches)
    
    
%      mode = 'threshold';
    mode = 'other';
    
    
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
    c_final = 1;
    change_ct = 0;
    switch_ct = 0;
    p_corr_thr = 0.95;
    p_dist_thr = 0.95;
    A_thr = 400;
    
    disp('registering')
%      p_thr = 0.95;
    while c < length(cluster)
      
      %% first, remove matched ROIs from cluster
      [s w] = find(cluster(c).list);
      for i = 1:length(s)
        n = cluster(c).list(s(i),w(i));
        if data(s(i)).matched(n)
          cluster(c).list(s(i),w(i)) = 0;
        end
      end
      
      if nnz(cluster(c).list) < 2   %% only look at clusters, that actually have some matching possibilities
        cluster(c) = [];
      else
        c_final = c_final + 1;
        n_ref = 0;
        s_ref = 0;
        
        %% merge status: for every neuron in the final_list, have 3 entries: previous, current and following session match status
        %% match status does not refer to matching to a certain neuron, but rather assigning to this cluster!
        cluster_final(c_final) = struct('list',zeros(nSes,1),'fp_corr',zeros(nSes,1),'dist',zeros(nSes,1),'p_same',zeros(nSes,1),'merge_status',false(nSes,3),'ct',0);
        
        for s = 1:nSes
        
          if any(cluster(c).list(s,:))
            
            %% compare to last registered neuron (closest in time)
            %% also, compare to other ones if no fit found (or to overall cluster?)
            
            if n_ref == 0   %% register new ROI as reference ROI
              %%% missing here: no merging in first session possible
              n = cluster(c).list(s,cluster(c).list(s,:)>0);
              n = n(1);
              cluster_final(c_final).list(s,1) = n;
              cluster_final(c_final).p_same(s,1) = -1;
              
              n_ref = n;
              s_ref = s;
            else
              
              if strcmp(mode,'threshold')
                [matches_s, p_same] = get_matches(cluster(c),xdata,model.type,1-p_thr,s_ref,n_ref,s);
                [p_corr_match,idx_match] = max(p_same);
                if p_corr_match > p_thr
                  best_match_s = matches_s(idx_match);
                  
                  %% check reciprocity
                  for sm = s-1:-1:s_ref
                    [matches_s_ref, p_same_recip] = get_matches(cluster(c),xdata,model.type,1-p_thr,s,best_match_s,sm);
                    [p_same_recip,idx_recip] = max(p_same_recip);
                    
                    
                    %%% now, add possibility to find several ROIs within one session
                    if p_same_recip > p_thr
                      n_match_recip = matches_s_ref(idx_recip);
                      
                      if sm==s_ref && n_match_recip == n_ref
                        cluster_final(c_final).list(s,1) = best_match_s;
                        
%                          cluster_final(c_final).fp_corr(s,1) = xdata(s_ref,s).fp_corr(n_match_recip,best_match_s);
%                          cluster_final(c_final).dist(s,1) = xdata(s_ref,s).dist(n_match_recip,best_match_s);
%                          cluster_final(c_final).p_same(s,1) = p_same_recip;
                        
                        n_ref = best_match_s;
                        s_ref = s;
                        
                      else  %% break if ROI is found to belong to other cluster
                        break
                      end
                    end
                  end
                  %% if went all the way through here without finding any partner, go to next session
                else
                  continue
                end
              
              
              
              else
%                  disp('----------------- new ---------------')
                [matches_s, p_same_match] = get_matches(cluster(c),xdata,model.type,1-p_thr,s_ref,n_ref,s);
                [~,idx_match] = max(p_same_match);
                best_match_s = matches_s(idx_match);
                
                for i=1:length(matches_s)
                  [matches_s_ref, p_same_recip] = get_matches(cluster(c),xdata,model.type,1-p_thr,s,matches_s(i),s_ref);
                  [~,idx_recip] = max(p_same_recip);
                  
                  %% what, if 2 or more occurences happen?
                  %% -> should gather all possible matches first, then decide which one to take, based on probability
                  
                  if matches_s(i) == best_match_s && matches_s_ref(idx_recip) == n_ref    %% if they are each others favorites
                    %% additionally check, whether this probability is larger than ... something?!
                    if p_same_recip(idx_recip) > 0.05
                      cluster_final(c_final).merge_status(s_ref,3) = true;
                      cluster_final(c_final).merge_status(s,1) = true;
                      
                      if ~ismember(matches_s(i),cluster_final(c_final).list(s,:))
                        idx = nnz(cluster_final(c_final).list(s,:)) + 1;
                        cluster_final(c_final).list(s,idx) = matches_s(i);
                      end
                    end
                    
                  elseif matches_s(i) == best_match_s                 %% if chosen ROI rather matches with another one
                  %% do not match!! (or rather: how much different are they? look for merging possibility?)
                  %%% here, should check for 2nd best match
                    if (p_same_recip(idx_recip) - p_same_match(idx_match) > 0.5)  %% really wants to match another one -> do not include in this cluster
                      change_ct = change_ct + 1;
%                        disp(sprintf('get me outa here, p_recip: %6.4g, p_match: %6.4g, change_ct=%d',p_same_recip(idx_recip),p_same_match(idx_match),change_ct))
                    else                                           %% if there might be a chance of both matching -> merge?
                      cluster_final(c_final).merge_status(s_ref,2:3) = true;
                      cluster_final(c_final).merge_status(s,1) = true;
                      merge_ct = merge_ct + 1;
%                        disp(sprintf('merge me!, p_recip: %6.4g, p_match: %6.4g, merge_ct=%d',p_same_recip(idx_recip),p_same_match(idx_match),merge_ct))
                      if ~ismember(matches_s_ref(idx_recip),cluster_final(c_final).list(s_ref,:))
                        idx = nnz(cluster_final(c_final).list(s_ref,:)) + 1;
                        cluster_final(c_final).list(s_ref,idx) = matches_s_ref(idx_recip);
                      end
                      if ~ismember(matches_s(i),cluster_final(c_final).list(s,:))
                        idx = nnz(cluster_final(c_final).list(s,:)) + 1;
                        cluster_final(c_final).list(s,idx) = matches_s(i);
                      end
                    end
                    
                  elseif matches_s_ref(idx_recip) == n_ref
  %                    disp('new lover found')
                    if (p_same_recip(idx_recip) - p_same_match(idx_match) > 0.5)  %% if probabilities far exceed, change match
                      cluster_final(c_final).merge_status(s_ref,3) = true;
%                        cluster_final(c_final).list(s,:) = 0;
                      if ~ismember(matches_s(i),cluster_final(c_final).list(s,:))
                        idx = nnz(cluster_final(c_final).list(s,:)) + 1;
                        cluster_final(c_final).list(s,idx) = matches_s(i);
                        cluster_final(c_final).merge_status(s,1) = true;
                      end
                      switch_ct = switch_ct + 1;
%                        disp(sprintf('switch!, p_recip: %6.4g, p_match: %6.4g, merge_ct=%d',p_same_recip(idx_recip),p_same_match(idx_match),switch_ct))
                    else
                      
                      merge2_ct = merge2_ct + 1;
%                        disp(sprintf('merge!, p_recip: %6.4g, p_match: %6.4g, merge_ct=%d',p_same_recip(idx_recip),p_same_match(idx_match),merge2_ct))
                      cluster_final(c_final).merge_status(s,2) = true;
                      if ~ismember(matches_s(i),cluster_final(c_final).list(s,:))
                        idx = nnz(cluster_final(c_final).list(s,:)) + 1;
                        cluster_final(c_final).list(s,idx) = matches_s(i);
                      end
                    end
                    best_match_s = matches_s(i);
                  end
                  
                end
                
%                  if any(cluster_final(c_final).list(s,:))
                  %% allow more than one neuron to go here
%                    n_ref = best_match_s;   %% this should include merging possibilities
%                    s_ref = s;
%                  end
              end
            end
          end
        end
        
        %% here, implement checking for ROI score and removing/merging/splitting accordingly
        
        
        
        
        
        
        
%          if ~strcmp(mode,'threshold')
          %% now, check whether all ROIs are a good match to each other (not just from session to session)
%            %% do this before "merging", to prohibit false matches to enhance merging chance
%            for p = 1:2
%              [s w] = find(cluster_final(c_final).list);
%              skip_arr = [];
%              for i = 1:length(s)
%                n_ref = cluster_final(c_final).list(s(i),w(i));
%                prob = zeros(length(s),1);
%                for j = 1:length(s)
%                  if s(j)==s(i) || any(skip_arr==j)
%                    prob(j) = -1;
%                  else
%                    n = cluster_final(c_final).list(s(j),w(j));
%                    prob(j) = xdata(s(i),s(j)).p_same_joint(n_ref,n);
%                  end
%                end
%                prob(prob==-1) = [];
%                if any(prob < 0.05) && ~cluster_final(c_final).merge_status(s(i),2)
%      %           cluster_final(c_final).list(s(i),w(i)) = 0;
%                  skip_arr = [skip_arr i];
%                elseif any(prob < 0.05)   %% check whether merging enhances matching probability -> do this with an overall score for the cluster instead of single checks and decisions
%      %              merge_status = true(length(s));
%                  for j = 1:length(s)
%                    if s(j)~=s(i) && ~any(skip_arr==j)
%                      n = cluster_final(c_final).list(s(j),w(j));
%                      merge_status = check_merge(data,xdata,s(j),n,s(i),cluster_final(c_final).list(s(i),cluster_final(c_final).list(s(i),:)>0));
%                      if ~merge_status
%      %                    disp('remove that!')
%                        cluster_final(c_final).list(s(i),w(i)) = 0;
%                        skip_arr = [skip_arr i];
%                        break
%                      end
%                    end
%                  end
%                  if j==length(s)
%                    merge_ct_real = merge_ct_real + 1;
%                  end
%                end
%              end
%            end
        
        
        
        
        
        
          %% score: high average and minimum probability, bias towards large number of neurons in one cluster
          %% check: removing one ROI from cluster: does it increase or decrease the "score"?
          %% or: possible to create subset from cluster that has high average and minimum probability?
          
          %% now, try merging if needed
%            if any(cluster_final(c_final).merge_status(:,2))
%              merge_ct = merge_ct + sum(cluster_final(c_final).merge_status(:,2));
  %            check_merge(data,xdata,
  %            disp('try if merging enhances results')
%            end
          
        
        
        %% only now, after removing "substandard matches" from the cluster, assign "matched" status to all
          
          %% assign matched status to all within cluster
          [s w] = find(cluster_final(c_final).list);
          for i = 1:length(s)
            n = cluster_final(c_final).list(s(i),w(i));
  %            [s(i) n]
            data(s(i)).matched(n) = true;
          end
          
          if nnz(cluster_final(c_final).list) < 1
            cluster_final(c_final) = [];
            c_final = c_final - 1;
          end
        
          cluster_final(c_final).ct = nnz(cluster_final(c_final).list(:,1));
          
          %% in the end, create new cluster from remaining ROIs and append to cluster struct
          if nnz(cluster_final(c_final).list) ~= nnz(cluster(c).list)
            
            c_idx = length(cluster)+1;
            for s = 1:nSes
              n = cluster(c).list(s,cluster(c).list(s,:)>0);
              n = n(~data(s).matched(n));
              cluster(c_idx).list(s,1:length(n)) = n;
            end
            if nnz(cluster(c_idx).list) < 2
  %              disp('purging needed?')
              cluster(c_idx) = [];
            end
          end
%          end
        
        c = c+1;
      end
    end
    
    nMatches = [cluster_final.ct];
    disp(sprintf('number of clusters: %d',c_final))
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
    
    path = sprintf('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/hist_s=%d.jpg',nSes);
    saveas(fig_ses,path,'jpg')
    
    disp(sprintf('saved under %s',path))
    
end



function [n, p_same] = get_matches(cluster,xdata,model_type,p_thr,s_ref,n_ref,s)
  
  %% search for all ROIs, that are not certainly rejected due to footprint or distance
  n = cluster.list(s,cluster.list(s,:)>0);
  
  switch model_type
    case 'dist'
      p_same = full(xdata(s_ref,s).p_same_dist(n_ref,n));
    case 'corr'
      p_same = full(xdata(s_ref,s).p_same_corr(n_ref,n));
    case 'joint'
      p_same = full(xdata(s_ref,s).p_same_joint(n_ref,n));
  end
  
  mask = p_same>p_thr;
  n = n(mask);
  p_same = p_same(mask);
  
end


function [ROI_score] = get_ROI_score(cluster,res,idx_rm)
  
  nSes = size(res.prob,1);
  mask = true(nSes);
  if idx_rm>0
    mask(idx_rm,:) = false;
    mask(:,idx_rm) = false;
  end
  
  p_mean = nanmean(res.prob(mask));
  p_var = nanvar(res.prob(mask));
  p_thr = min(0.9,p_mean-2*p_var);
  disp(sprintf('mean correlation %4.2f, variance %6.4f, threshold: %4.2f',p_mean,p_var, p_thr))
  
  fp_corr_1w_1 = res.fp_corr_oneway(:,:,1);
  fp_corr_1w_2 = res.fp_corr_oneway(:,:,2);
  
  res.frac_1w = fp_corr_1w_1(mask)./fp_corr_1w_2(mask);
  res.frac_1w(res.frac_1w>1) = 1./res.frac_1w(res.frac_1w>1);
  frac_1w_mean = nanmean(res.frac_1w(:));
  frac_1w_var = nanvar(res.frac_1w(:));
  
  ROI_score_raw = 1/3*(p_mean^(1+p_var) + frac_1w_mean^(1+frac_1w_var) + nanmin(res.prob(:)));
  ROI_score = 1/3*(p_mean^(1+p_var) + frac_1w_mean^(1+frac_1w_var) + nanmin(res.prob(:)))*sum(~isnan(nanmean(res.prob)))/nSes;
  disp(sprintf('ROI score: %4.2g vs %4.2g',ROI_score_raw,ROI_score))
end