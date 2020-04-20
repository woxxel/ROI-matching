

function [c_matches,match] = analyze_matchStats(pathMouse,std,thr,w,OnACID,old)
  
%    rmpath('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Programs/OASIS_matlab-master/optimization/cvx/lib/narginchk_:')
  if nargin < 2
    old = false;
  end
  
  match = struct;
  if OnACID
    suffix = '_OnACID';
  else
    suffix = '';
  end
  
  if ~old
    pathAss = pathcat(pathMouse,sprintf('matching/results_matching_multi_std=%d_thr=%d_w=%d%s.mat',std,thr*100,round(w*100),suffix))
    load(pathAss,'scores','assignments')
  %    load(pathcat(pathMouse{1},'xdata.mat'))
    
    for s = 2:15
      N_union = max(find(~isnan(assignments(:,s-1))));
      N_s = sum(~isnan(assignments(:,s)));
      match(s).bool = zeros(N_union,N_s)*NaN;
      match(s).d_s = zeros(N_union,N_s)*NaN;
%        end
    end
  end
  
  load(pathcat(pathMouse,'backup/save_final/matching_results.mat'))
  matchStats_arr{1} = get_matchState(data,clusters_sv);
  
  if old
%      for i = 1:length(pathMouse)    
    load(pathcat(pathMouse,'matching_results.mat'))
    [matchStats_arr{2}, data2, idx_ref_match, idx_match_ref] = OnACID_get_matchState(pathMouse,std,thr,w,OnACID,old);
%      matchStats_arr{2} = get_matchState(data,clusters_sv);
  else
    [matchStats_arr{2}, data2, idx_ref_match, idx_match_ref] = OnACID_get_matchState(pathMouse,std,thr,w,OnACID,old);
  end
    
  col = {'b','r'};
  
  nSes = length(matchStats_arr{1}.clusters(1).list);
  
  figure('position',[100 100 800 400])
  
  nROI_base = zeros(nSes,1);
  for s = 1:nSes
%      [length(idx_ref_match{s}),length(idx_match_ref{s})]
    nROI_base(s) = length(idx_ref_match{s});%data.session(s).nROI;
%      disp(sprintf('nROIs: %d , %d',sum(~isnan(idx_ref_match{s})),sum(~isnan(idx_match_ref{s}))))
    data.session(s).nROI = sum(~isnan(idx_match_ref{s}));
  end
  ROI_recurr_frac = zeros(nSes-1,2);
  
  %% clean up from non-matched neurons
  for s=1:nSes
    for n = 1:length(idx_ref_match{s})
      if isnan(idx_ref_match{s}(n)) && n <= length(matchStats_arr{1}.session(s).neuron)
        c = matchStats_arr{1}.session(s).neuron(n).cluster_ID;
        if ~isempty(c)
          matchStats_arr{1}.session(s).neuron(n).cluster_ID = [];
          matchStats_arr{1}.clusters(c).list(s) = 0;
        end
      end
    end
    
    for n = max(idx_match_ref{s})+1:length(matchStats_arr{1}.session(s).neuron)
      c = matchStats_arr{1}.session(s).neuron(n).cluster_ID;
      matchStats_arr{1}.session(s).neuron(n).cluster_ID = [];
      if ~isempty(c)
        matchStats_arr{1}.clusters(c).list(s) = 0;
      end
    end
    
%      for n = 1:length(idx_match_ref{s})
%        n
%        if isnan(idx_match_ref{s}(n))% && n <= length(matchStats_arr{2}.session(s).neuron)
%          c = matchStats_arr{2}.session(s).neuron(n).cluster_ID;
%          if ~isempty(c)
%            matchStats_arr{2}.session(s).neuron(n).cluster_ID = [];
%            matchStats_arr{2}.clusters(c).list(s) = 0;
%          end
%        end
%      end
%      
%      for n = max(idx_ref_match{s})+1:length(matchStats_arr{2}.session(s).neuron)
%        c = matchStats_arr{2}.session(s).neuron(n).cluster_ID;
%        matchStats_arr{2}.session(s).neuron(n).cluster_ID = [];
%        if ~isempty(c)
%          matchStats_arr{2}.clusters(c).list(s) = 0;
%        end
%      end
  end
  
  for i = 1:length(matchStats_arr)
    nC = length(matchStats_arr{i}.clusters);
    
    cluster_nROI = zeros(nSes,1);
    
    nROI = zeros(nSes,1);
    ROI_recurr = zeros(nSes,nSes-1);
    
    for c = 1:nC
      
      matchStats_arr{i}.clusters(c).nROI = nnz(matchStats_arr{i}.clusters(c).list);
      
      if matchStats_arr{i}.clusters(c).nROI
        cluster_nROI(matchStats_arr{i}.clusters(c).nROI) = cluster_nROI(matchStats_arr{i}.clusters(c).nROI) + 1;
      end
      
      nROI = nROI + (matchStats_arr{i}.clusters(c).list>0);
    end
    
    nROI_total = zeros(nSes-1,1);
    for s = 1:nSes
      for sm = s+1:nSes
        ds = sm-s;
%          nROI_total(ds) = nROI_total(ds) + nROI(s);
        nROI_total(ds) = nROI_total(ds) + nROI_base(s);
      end
    end
    
    
%      subplot(4,2,1)
%      hold on
    
%      bar(nROI,'FaceAlpha',0.6,'FaceColor',col{i})
%      xlabel('Session #')
%      ylabel('# ROIs')
    
%      subplot(2,2,2)
%      hold on
%      bar(nROI_total,'FaceAlpha',0.6)
    
    subplot(2,2,1)
    if i == 1
      dispName = 'pre';
    else
      dispName = 'post';
    end
    hold on
    bar(cluster_nROI,'FaceAlpha',0.6,'FaceColor',col{i},'DisplayName',dispName)
    xlabel('# ROIs / cluster')
    ylabel('# clusters')
    ylim([0,400])
    legend('location','NorthEast')
    pause(0.0001)
    
    
    for s = 1:nSes
      
      for sm = s+1:nSes
        
        ds = sm-s;
        
        for c = 1:nC
          
          if matchStats_arr{i}.clusters(c).list(s) && matchStats_arr{i}.clusters(c).list(sm)
            ROI_recurr(s,ds) = ROI_recurr(s,ds) + 1;
          end
          
        end
        
      end
      
    end
    
    ROI_recurr_frac(:,i) = sum(ROI_recurr,1)./nROI_total';
    
    subplot(2,2,2)
    hold on
    plot([0,15],[0.5,0.5],'k:','LineWidth',0.5,'HandleVisibility','off')
    plot(ROI_recurr_frac(:,i),'o--','DisplayName','p_{recurr}')
    xticks([])
%      ylim([0,1])
    ylabel('p_{rec}')  
    
  end
  
  subplot(2,2,2)
  hold on
  plot([0,15],[0,0],'k:','LineWidth',0.5,'HandleVisibility','off')
  plot(ROI_recurr_frac(:,2)-ROI_recurr_frac(:,1),'ko--','DisplayName','\Delta p')
  ylim([-0.1,1])
  xlabel('\Delta s')
  ylabel('\Delta p_{rec}')
  disp_txt = false;
  matches = zeros(nSes-1,4);
  xlim([0.5,15.5])
  legend('Location','NorthEast')
  
%    figure
  subplot(2,2,4)
  hold on
  p_matches_correct = plot(matches(:,1)./sum(matches,2),'g','DisplayName','correct match');
  p_matches_wrong = plot(matches(:,2)./sum(matches,2),'r','DisplayName','wrong match');
  p_matches_missing = plot(matches(:,3)./sum(matches,2),'r--','DisplayName','missing');
  p_matches_unnec = plot(matches(:,4)./sum(matches,2),'r:','DisplayName','unnec. match');
  
  legend('Location','northeastoutside')
  xlabel('\Delta s')
  yticks(linspace(0,1,6))
  yticklabels(linspace(0,1,6))
  ylabel('fraction')
  c_matches = zeros(length(matchStats_arr{1}.clusters),4); %% col-# = cluster-ID, cols: 1: new cluster-ID; 2: old #ROIs, 3: new #ROIs, 4: # correct ROIs
  xlim([0.5,15.5])
  ylim([0,1])
  pause(0.0001)
  
  
  nROI_match = zeros(nSes+1,4);
  
  %% now, compare both matchStats
  for s = 1:nSes
    disp(sprintf('Now processing s=%d',s))
    for n = 1:length(idx_ref_match{s})
      
      if isnan(idx_ref_match{s}(n)) || n > length(matchStats_arr{2}.session(s).neuron)
%          disp('why?')
        continue
      end
      
      
      c2 = matchStats_arr{2}.session(s).neuron(n).cluster_ID;
      c_ref = matchStats_arr{1}.session(s).neuron(n).cluster_ID;
      
      
%        matchStats_arr{1}.clusters(c_ref).list
%        matchStats_arr{2}.clusters(c2).list
      
      
      if ~isempty(c_ref)
        if isempty(c2)
          if ~c_matches(c_ref,4)
            c_matches(c_ref,1) = NaN;
            c_matches(c_ref,2) = nnz(matchStats_arr{1}.clusters(c_ref).list);
            c_matches(c_ref,3) = 0;
            c_matches(c_ref,4) = 0;
          end
        else
          init_correct = sum(matchStats_arr{1}.clusters(c_ref).list & matchStats_arr{1}.clusters(c_ref).list==matchStats_arr{2}.clusters(c2).list);
          if init_correct > c_matches(c_ref,4);
            c_matches(c_ref,1) = c2;
            c_matches(c_ref,2) = nnz(matchStats_arr{1}.clusters(c_ref).list);
            c_matches(c_ref,3) = nnz(matchStats_arr{2}.clusters(c2).list);
            c_matches(c_ref,4) = init_correct;
          end
        end
      end
      
%        if ~isempty(c_ref)% && nnz(matchStats_arr{1}.clusters(c_ref).list) < 2
%          continue
%        end
      
%        if ~isempty(c2)% && nnz(matchStats_arr{2}.clusters(c2).list) < 2
%          continue
%        end
      
%        init_correct
%        waitforbuttonpress
      
%        if ~isempty(c2) && matchStats_arr{i}.clusters(c2).nROI == 15
          
      for sm = s+1:nSes
        
        ds = sm-s;
        
        if disp_txt
          [s sm ds]
          n
          [c_ref,c2]
          matchStats_arr{1}.clusters(c_ref).list
          matchStats_arr{2}.clusters(c2).list
          
%            [matchStats_arr{1}.clusters(c_ref).list(s) matchStats_arr{2}.clusters(c2).list(s)]
%            [matchStats_arr{1}.clusters(c_ref).list(sm) matchStats_arr{2}.clusters(c2).list(sm)]
        end
        
%          m_true = [];
%          if ~isempty(matchStats_arr{2}.clusters(c2).list(sm))
%            m_true = matchStats_arr{2}.clusters(c2).list(sm);
%          end
        
        if isempty([c_ref,c2])
          if disp_txt
            disp('none assigned')
          end
        elseif isempty(c_ref)
          m2 = matchStats_arr{2}.clusters(c2).list(sm);
          
          if m2
            matches(ds,4) = matches(ds,4) + 1;
            if ~old
              match(sm).bool(c2,idx_ref_match{sm}(m2)) = 4;
              match(sm).d_s(c2,idx_ref_match{sm}(m2)) = sm-s;
              
              nROI = 16;
              nROI_match(nROI,4) = nROI_match(nROI,4) + 1;
            end
          end
          if disp_txt
            disp('ref n not assigned')
          end
        elseif isempty(c2)
          m = matchStats_arr{1}.clusters(c_ref).list(sm);
          if m && ~isnan(idx_ref_match{sm}(m))
            matches(ds,3) = matches(ds,3) + 1;
            if ~old
              match(sm).bool(c2,idx_ref_match{sm}(m)) = 3;
              
              nROI = matchStats_arr{1}.clusters(c_ref).nROI;
              nROI_match(nROI,3) = nROI_match(nROI,3) + 1;
            end
          end
          
          if disp_txt
            disp('processed n not assigned')
          end
        else
          
          m2 = matchStats_arr{2}.clusters(c2).list(sm);
          m = matchStats_arr{1}.clusters(c_ref).list(sm);
%            if m<=length(idx_ref_match{sm}) && (~m || ~isnan(idx_ref_match{sm}(m))) && m2<=length(idx_ref_match{sm}) && (~m2 || ~isnan(idx_ref_match{sm}(m2)))
          if ~m && ~m2
            if disp_txt
              disp('nothing, both are empty!')
            end
          elseif ~m
            matches(ds,4) = matches(ds,4) + 1;
            if ~old
              match(sm).bool(c2,idx_ref_match{sm}(m2)) = 4;
              match(sm).d_s(c2,idx_ref_match{sm}(m2)) = sm-s;
              
              nROI = matchStats_arr{1}.clusters(c_ref).nROI;
              nROI_match(nROI,4) = nROI_match(nROI,4) + 1;
            end
            if disp_txt
              disp('ref sm not assigned')
            end
          elseif ~m2
            matches(ds,3) = matches(ds,3) + 1;
            if ~old
              match(sm).bool(c2,idx_ref_match{sm}(m)) = 3;
              
              nROI = matchStats_arr{1}.clusters(c_ref).nROI;
              nROI_match(nROI,3) = nROI_match(nROI,3) + 1;
            end
            
            if disp_txt
              disp('processed sm not assigned')
            end
          elseif m2 ~= m
            matches(ds,2) = matches(ds,2) + 1;
            if ~old
              match(sm).bool(c2,idx_ref_match{sm}(m2)) = 2;
              match(sm).bool(c2,idx_ref_match{sm}(m)) = 3;
              match(sm).d_s(c2,idx_ref_match{sm}(m2)) = sm-s;
              
              nROI = matchStats_arr{1}.clusters(c_ref).nROI;
              nROI_match(nROI,2) = nROI_match(nROI,2) + 1;
              
            end
            if disp_txt
              disp('wrongly assigned')
            end
          else
            matches(ds,1) = matches(ds,1) + 1;
            if ~old
              match(sm).bool(c2,idx_ref_match{sm}(m2)) = 1;
              match(sm).d_s(c2,idx_ref_match{sm}(m2)) = sm-s;
              
              nROI = matchStats_arr{1}.clusters(c_ref).nROI;
              nROI_match(nROI,1) = nROI_match(nROI,1) + 1;
            end
            if disp_txt
              disp('properly assigned')
            end
          end
%            end
        end
%          if disp_txt
%            waitforbuttonpress
%          end
      end
      
      if ~mod(n,200)
        norm_matches = sum(matches(:,1:3),2);
        set(p_matches_correct,'YData',matches(:,1)./norm_matches)
        set(p_matches_wrong,'YData',matches(:,2)./norm_matches)
        set(p_matches_missing,'YData',matches(:,3)./norm_matches)
        set(p_matches_unnec,'YData',matches(:,4)./norm_matches)
        drawnow
%        pause(0.0001)
%          waitforbuttonpress
      end
    end
  end
  
  subplot(2,2,3)
  scatter(c_matches(:,3)+0.5*rand(size(c_matches,1),1)-0.25,c_matches(:,2)+0.5*rand(size(c_matches,1),1)-0.25,5,c_matches(:,4)./c_matches(:,3));
  colormap('jet');
  
  xlim([1.5,15.5])
  ylim([1.5,15.5])
  xlabel('#ROIs pre')
  ylabel('#ROIs post')
  cbar = colorbar;
  
  ylabel(cbar, 'ratio correct matches')
  
  
  pathSv = pathcat(pathMouse,'Figures',sprintf('results_matching_multi_std=%d_thr=%d_w=%d%s.png',std,thr*100,round(w*100),suffix))
%    path = pathcat(pathFigures,'results_manual_matching.png');
  print(pathSv,'-dpng','-r300')
  
  figure
  c_idx = c_matches(:,2)>1;
%    d_matches = c_matches(c_idx,2)-c_matches(c_idx,4);
  d_matches = c_matches(c_idx,2)-c_matches(c_idx,3);
  histogram2(c_matches(c_idx,2),d_matches,linspace(-0.5,15.5,16),linspace(-15.5,15.5,31))
    
    
    
  if ~old
    disp('now, find avg and STD of corr. and distance of wrongly assigned matches / missing matches, as well as same distr from proper matches (able to characterize those based on those values?). also, find the shifted fp_corr of those groups (better match with this?). can i find a better value, which describes the matching yet better?')
    
    disp('get something, saying how many clusters are actually estimated right, vs wrong, + how many errors are made in the clusters (per nROI?)')
    
    
    %%% gather scores for each entry in cluster, compute average or something and plot vs # of correct matches
    
    
    
    
%      
    col = {'g','r','b','m'}
%      col = {[0,1,0],[1,0,0],[0,0,1],[1,0,1]}
    dispName = {'correct','wrong','missing','unnecessary'}
    figure
    
    matches_ecdf = cell(4,1);
    
    nROI_match
    subplot(1,3,2)
    h = bar(nROI_match,'stacked');
    for i=1:4
      set(h(i),'FaceColor',col{i})
    end
%      set(h,{'FaceColor'},col)
    hold on
    
    subplot(1,3,3)
    h = bar(transpose(nROI_match),'stacked');
%      for i=1:4
%        set(h(i),'FaceColor',col{i})
%      end
    hold on
    
    for i = 1:4
      
      
      for s = 2:nSes
  %        for sm = s+1:nSes
  %        disp(sprintf('session %d',s))
        mask = logical(match(s).bool==i); %(1:size(scores{s,1},1),1:size(scores{s,1},2))
        
        matches_ecdf{i} = horzcat(matches_ecdf{i},transpose(scores{s-1,1}(mask)));
        
%          subplot(2,2,1)
%          hold on
%          scatter(scores{s-1,1}(mask),scores{s-1,3}(mask),'Color',col{i})
  %        end
      end
      subplot(1,3,1)
      [f,x,flo,fup] = ecdf(matches_ecdf{i});
      stairs(x,f,'Color',col{i},'DisplayName',dispName{i});
      hold on
      stairs(x,flo,':','Color',col{i},'HandleVisibility','off');
      stairs(x,fup,':','Color',col{i},'HandleVisibility','off');
      
%        subplot(1,2,2)
%        hold on
%        hist(matches_ecdf{i})%,'FaceAlpha',0.5,'FaceColor',col{i})
    end
    
    legend('Location','NorthWest')
    xlim([0,1])
    
    
    pathSv = pathcat(pathMouse,'Figures',sprintf('results_matching2_multi_std=%d_thr=%d_w=%d%s.png',std,thr*100,round(w*100),suffix))
  %    path = pathcat(pathFigures,'results_manual_matching.png');
    print(pathSv,'-dpng','-r300')
  end
  
%    for i = -1:1
%      for s = 1:nSes
%        for sm = s+1:nSes
%          mask = match(s,sm).bool(1:size(xdata(1,2).corr,1),1:size(xdata(1,2).corr,2))==i;
%          scatter(xdata(s,sm).dist(mask),xdata(s,sm).corr(mask),col{i})
%        end
%      end
%    end
  
  
  
  
end