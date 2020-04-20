  
  
function [matchState, data, idx_ref_match, idx_match_ref] = OnACID_get_matchState(pathMouse,std,thr,w,OnACID,old)
  
  nSes = 15;
  if nargin==6 && old
    pathMatching = pathcat(pathMouse,'matching_results.mat');
    load(pathMatching,'data','clusters_sv')
    nC = length(clusters_sv);
    assignments = zeros(nC,nSes);
    for c=1:nC
      for s=1:nSes
        if ~isempty(clusters_sv(c).session(s).list)
          assignments(c,s) = clusters_sv(c).session(s).list;
        else
          assignments(c,s) = NaN;
        end
      end
    end
  else
    if OnACID
      suffix = '_OnACID';
    else
      suffix = '';
    end
    pathAss = pathcat(pathMouse,sprintf('matching/results_matching_multi_std=%d_thr=%d_w=%d%s.mat',std,thr*100,round(w*100),suffix))
    load(pathAss,'assignments')   %% loads assignments, matchings, scores, shifts
    data = struct('nSes',nSes,'session',struct('nROI',NaN,'ROI',struct('cluster_ID',NaN)));
    nC = size(assignments,1);
  end
  %% clean of single entries
%    size(assignments)
%    assignments(find(sum(~isnan(assignments),2)<2),:) = [];
%    size(assignments)
  
  idx_ref_match = {};
  idx_match_ref = {};
  for s = 1:nSes
    if old
      [idx_ref_match{s} idx_match_ref{s}] = match_by_Ca_corr(pathMouse,s);
%        idx_ref_match{s} = 1:length(idx_match_ref{s});
%        idx_match_ref{s} = 1:length(idx_match_ref{s});
    elseif ~OnACID
%        pathNewCa = pathcat(pathMouse,'backup/save_final',sprintf('Session%02d',s),'CaData.mat');
%        CaData = matfile(pathNewCa);
      N = nanmax(assignments(:,s))+1;
      idx_ref_match{s} = 1:N;
      idx_match_ref{s} = 1:N;
      
%        [idx_ref_match{s} idx_match_ref{s}] = match_by_Ca_corr(pathMouse,s);
%        length(idx_match_ref{s})
%        idx_ref_match{s} = 1:length(idx_match_ref{s});
%        idx_match_ref{s} = 1:length(idx_match_ref{s});
%        [idx_ref_match{s} idx_match_ref{s}] = match_by_Ca_corr(pathMouse,s);
    else
      [idx_ref_match{s} idx_match_ref{s}] = match_old_new_IDs(pathMouse,s);
      
%        [idx_CNMFr_ACID idx_ACID_CNMFr] = match_old_new_IDs(pathMouse,s);
%        [ref_CNMFr, CNMFr_ref] = match_by_Ca_corr(pathMouse,s);
%        
%        idx_match_ref{s} = zeros(length(idx_ACID_CNMFr),1)*NaN;
%        for i = 1:length(idx_ACID_CNMFr)
%          if ~isnan(idx_ACID_CNMFr(i))
%            idx_match_ref{s}(i) = CNMFr_ref(idx_ACID_CNMFr(i));
%          end
%        end
%        
%        idx_ref_match{s} = zeros(length(ref_CNMFr),1)*NaN;
%        for i = 1:length(ref_CNMFr)
%          if ~isnan(ref_CNMFr(i))
%            idx_ref_match{s}(i) = idx_CNMFr_ACID(ref_CNMFr(i));
%          end
%        end
      
%        disp('ACID -> ref')
%        length(CNMFr_ref)
%        length(idx_ACID_CNMFr)
%        length(idx_match_ref{s})
%        
%        disp('ref -> ACID')
%        length(ref_CNMFr)
%        length(idx_CNMFr_ACID)
%        length(idx_ref_match{s})
      
    end
  end
  
  for s = 1:nSes
    data.session(s).nROI = 0;
%      data.session(s).ROI(length(idx_match_ref{s})).cluster_ID = NaN;
    for c = 1:nC
      n = assignments(c,s) + ~old;
      if ~isnan(n)
%          n
        if ~isnan(idx_match_ref{s}(n))% && sum(~isnan(assignments(c,:)),2)>2
          data.session(s).ROI(idx_match_ref{s}(n)).cluster_ID = c;
          data.session(s).nROI = data.session(s).nROI + 1;
%          else
%            data.session(s).ROI(idx_match_ref{s}(n)).cluster_ID = [];
        end
      end
    end
  end
  
  
  matchState = get_matchState(data,nC);
  
  


function [old_new, new_old] = match_old_new_IDs(pathMouse,session)
  
  pathSession = pathcat(pathMouse,sprintf('Session%02d',session));
  pathMatchResults = pathcat(pathSession,'matching_old_new.mat');
  
  load(pathMatchResults)
  
  N1 = length(matched_ROIs1) + length(non_matched1);
  N2 = length(matched_ROIs2) + length(non_matched2);
  
  old_new = zeros(N1,1)*NaN;
  new_old = zeros(N2,1)*NaN;
  
  for i = 1:length(matched_ROIs1)
    old_new(matched_ROIs1(i)+1) = matched_ROIs2(i)+1;
  end
  
  for i = 1:length(matched_ROIs2)
    new_old(matched_ROIs2(i)+1) = matched_ROIs1(i)+1;
  end
  

function [proc_raw, raw_proc] = match_by_Ca_corr(pathMouse,session)
  pathSession = pathcat(pathMouse,sprintf('backup/save_final/Session%02d',session));
  pathSv = pathcat(pathSession,'match_cleaned.mat');
  
  if exist(pathSv,'file')
    disp(sprintf('Path %s already present, loading results',pathSv))
    load(pathSv)
  else
    pathOldCa = pathcat(pathSession,'CaData.mat');
    pathNewCa = pathcat(pathMouse,'backup/save_final',sprintf('Session%02d',session),'CaData.mat');
    
    load(pathOldCa,'C2');
    oldCa = C2;
    load(pathNewCa,'C2');
    newCa = C2;
    
    raw_proc = zeros(size(oldCa,1),1)*NaN;
    proc_raw = zeros(size(newCa,1),1)*NaN;
    
    for n = 1:size(newCa,1)
      found = false;
      i=0;
      while ~found && n+i <= size(oldCa,1)
        tmp_corr = corrcoef(newCa(n,:),oldCa(n+i,:));
        if tmp_corr(1,2) > 0.99
          found = true;
        else
          i = i+1;
        end
      end
      if found
        proc_raw(n) = n+i;
        raw_proc(n+i) = n;
      else
        disp(sprintf('no more matches found @%d of %d neurons',n,size(newCa,1)))
        break
      end
      save(pathSv,'raw_proc','proc_raw')
    end
    
    
  end
  