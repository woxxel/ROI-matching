

function ROI_matching_test(pathMouse)

  paths = struct;
  paths.mouse = pathMouse;
  paths.background_folder = 'Session01';
  paths.background = 'reduced_MF1_LK1.mat';
  paths.background_field = 'max_im';
  paths.footprints = 'resultsCNMF_MF1_LK1.mat';
  paths.footprints_field = 'A2';
  
  paths.xdata = pathcat(paths.mouse,'xdata.mat');
  paths.results = pathcat(paths.mouse,'matching_results.mat');
  
  footprints = match_loadSessions(paths)
  
%    if ~exist(paths.xdata,'file')
    [xdata, histo, para] = match_analyzeData(footprints,footprints.data.nSes,12);      %% calculate distances and footprint correlation
    [model,histo] = match_buildModel(xdata,histo,para,footprints.data.nSes,paths.mouse);
  %    [ROC] = estimate_model_accuracy(histo,model,para,pathMouse);
    
    %% and assigning probabilities to each (close) pair
    xdata = match_assign_prob(xdata,footprints.data,model,para,paths.xdata);
%    else
%      load(paths.xdata)
%    end
%    setappdata(0,'xdata',xdata)
  
%    bool_ld = get(h.uihandles.checkbox_load_processed_data,'Value') && exist(paths.results,'file');
%    if bool_ld
%      ld_data = load(paths.results);
%      clusters = ld_data.clusters_sv;
%      status = ld_data.status;
%    else
    real_matching(footprints.data,0.5);
%      clusters = getappdata(0,'clusters');
    status = [];
  end
  
end




function real_matching(data,p_thr)
  
  xdata = getappdata(0,'xdata');
  mode = 'threshold';
%      mode = 'other';
  
  c = 0;
  
  session(data.nSes) = struct;
  for s = 1:data.nSes
    session(s).ROI_matched = false(data.session(s).nROI,1);
  end
  
  disp('registering')
  
  for s = 1:data.nSes
    
    for n = 1:data.session(s).nROI
      
      if ~session(s).ROI_matched(n)
        
        c = c + 1;
        clusters(c) = struct('A',[],'centroid',[],'score',NaN,'ct',NaN,'session',struct('list',cell(data.nSes,1),'ROI',struct('score',[],'mean_score',[],'unsure',false)));
        
        clusters(c).session(s).list = n;
        
        s_ref = s;
        n_ref = n;
        
        if mod(c,500)==0
          disp(sprintf('%d clusters found.',c))
        end
        
        for sm = 1:data.nSes
          if s==sm
            continue
          end
          
          if strcmp(mode,'threshold')
            
%              matches_s = find(xdata(s_ref,sm).prob(n_ref,:)>p_thr);
            matches_s = find(xdata(s_ref,sm).prob(n_ref,:)>p_thr & ~session(sm).ROI_matched');
            p_same_s = xdata(s_ref,sm).prob(n_ref,matches_s);
            
            [p_best_s,idx_s] = max(p_same_s);
            if p_best_s > p_thr
              best_match_s = matches_s(idx_s);
              
              %% check for reciprocity
%                matches_s_ref = find(xdata(sm,s_ref).prob(best_match_s,:)>p_thr);
              matches_s_ref = find(xdata(sm,s_ref).prob(best_match_s,:)>p_thr & ~session(s_ref).ROI_matched');
              p_same_s_ref = xdata(sm,s_ref).prob(best_match_s,matches_s_ref);
              
              [p_best_s_ref,idx_s_ref] = max(p_same_s_ref);
              
              if (matches_s_ref(idx_s_ref) == n_ref) && (p_best_s_ref > p_thr)% && ~session(sm).ROI_matched(best_match_s)
                clusters(c).session(sm).list = best_match_s;
                
                s_ref = sm;
                n_ref = best_match_s;
              end
            end
          end
        end
        
%        %% assign matched status to all within pre_clusters
        occupancy = zeros(data.nSes,1);
        for sm = 1:data.nSes
          for i = 1:length(clusters(c).session(sm).list)
            n = clusters(c).session(sm).list(i);
            clusters(c).session(sm).ROI(i).unsure = false;
            session(sm).ROI_matched(n) = true;
          end
          occupancy(sm) = length(clusters(c).session(sm).list);
        end
        
        clusters(c).ct = nnz(occupancy);
        if clusters(c).ct < 2
          clusters(c) = [];
          c = c - 1;
        end
        
      end
    end
  end
  
  setappdata(0,'clusters',clusters)
  
end