

function [arrays,occupancy,ROI_recurr,N_pairs,N_norm,pop_overlap] = display_matching_stats(pathMouse,PCs,perc_thr,ct_thr,N_bs,s_offset,sv,sv_suffix,sv_ext,arrays,occupancy,ROI_recurr,N_pairs,N_norm,pop_overlap)%,ROI_rec2,ROI_tot2)%pathBase,mouse)
  
  nSes = size(PCs.status,2)
  
  %%% think about: doing all the plots in an UI, able to change parameters such as gate & reward location, sessions considered, etc
  %%% -> which plots are needed, which ones are not that important?
  
  
  %%% fitting functions and options
  F_est = @(x,data)...
              x(4)/(sqrt(2*pi)*x(2))*exp(-((data-x(1)).^2./(2*x(2)^2))) + x(3)/length(data);     %% gaussian + linear offset
  exp_dist = @(a,x_data) a(2)*exp(-x_data/a(1));
  log_exp_dist = @(a,x_data) log(a(2)) - x_data/a(1);
  
  options = statset('MaxIter',100, 'MaxFunEvals',200,'Display','off');
  
  
  %%% plotting options
  plt_presi = false;
  plot_pop = false;
  
  pathFigures = pathcat(pathMouse,'Figures');
  
  plot_arr = {'NRNG','GT','RW'};
  col = {'b','g','r'};
  col_fill = {[0.5,0.5,1],[0.5,1,0.5],[1,0.5,0.5]};
  
  s_start = 1;
  s_end = nSes;
  
  [~,mouse,~] = fileparts(pathMouse);
%% absorb this into get_t_measures!!!
%%% ----------------------------------------------
  time_real = false;
  if time_real
    t_measures = get_t_measures(mouse);
    t_mask_m = false(1,t_measures(nSes));
    for s = 1:nSes-1
      for sm = s+1:nSes
        dt = t_measures(sm)-t_measures(s);
        t_mask_m(dt) = true;
      end
    end
    t_data_m = find(t_mask_m);
    t_ses = t_measures;
    t_mask = t_mask_m;
    t_data = t_data_m;
    nT = length(t_data);
  else
    t_measures = get_t_measures(mouse);
    t_measures = t_measures(s_offset:s_offset+nSes-1)
%      t_measures = 1:nSes;    %% remove!
%      t_measures
    t_ses = linspace(1,nSes,nSes);
    t_data = linspace(1,nSes,nSes);
    t_mask = true(nSes,1);
    nT = nSes;
  end
%%% ----------------------------------------------


  %%% get this from paras
  para = set_paras([],[],mouse);
  
  PC_bar = zeros(para.nbin,1);
  PC_bar(para.idx_zone{3}) = 1;
  
  gt_bar = zeros(para.nbin,1);
  
  if length(para.idx_zone{4})>1
    gt_bar(para.idx_zone{4}) = 1;
  end
  
  rw_bar = zeros(para.nbin,1);
  rw_bar(para.idx_zone{5}) = 1;
  
  tolerance = 6./(100/para.nbin);
  
  
  %% remove clusters, consisting of single neurons, only
  idx_delete = sum(PCs.status(:,:,2),2)<2;
  
  PCs.clusterID(idx_delete,:) = [];
  PCs.status(idx_delete,:,:) = [];
  PCs.trials(idx_delete,:) = [];
  PCs.firingrate(idx_delete,:) = [];
  PCs.shifts(idx_delete) = [];
  
  PCs.fields.center(idx_delete,:,:) = [];
  PCs.fields.width(idx_delete,:,:) = [];
  PCs.fields.status(idx_delete,:,:) = [];
  PCs.fields.firingmap(idx_delete,:,:) = [];
  
  PCs.MI.value(idx_delete,:) = [];
  PCs.MI.p_value(idx_delete,:) = [];
%    PCs.MI.rand_distr(idx_delete,:,:) = [];
  
  nC = size(PCs.status,1);
  nC
  
  ROI_ct = sum(PCs.status(:,:,2),2);
  
  %% build blue-red colormap
  n = 51;   %% must be an even number
  cm = ones(3,n);
  cm(1,1:ceil(n/2)) = linspace(0,1,ceil(n/2));      %% red
  cm(2,1:ceil(n/2)) = linspace(0,1,ceil(n/2));      %% green
  cm(2,ceil(n/2):n) = linspace(1,0,floor(n/2)+1);   %% green
  cm(3,ceil(n/2):n) = linspace(1,0,floor(n/2)+1);   %% blue
  
  
  plot_fig = false(100,1);%[false,true,false,false,false,false,false,false,false,false,false];
%    plot_fig = [true,true,true,true,true,true,true,true,true,true];
%    plot_fig = [false,false,true,false,false,false,false,false,false,false];
  
%    plot_fig(1) = true;       %% basic cluster stats (sessions active / cluster, cluster score)
%    plot_fig(2) = true;       %% neuron numbers (active, silent, PCs (NRNG,GT,RW))
%  %  %  %  %    plot_fig(3) = true;       %% (old)
%    plot_fig(4) = true;       %% number of fields
%    plot_fig(5) = true;       %% positions of stable fields
%    plot_fig(6) = true;       %% remapping
%    plot_fig(7) = true;       %% trial threshold
%    plot_fig(8) = true;       %% MI per session
%    plot_fig(9) = true;       %% PC coverage
%    plot_fig(10) = true;      %% PC remapping (presi)
%    plot_fig(11) = true;      %% overrepresentation //correlation of activity / PCs
%    plot_fig(12) = true;      %% plot pie charts of neuron numbers / session
%    plot_fig(13) = true;      %% firing rates
%    plot_fig(14) = true;      %% (no plot) pie charts: silent in- and out- transitions
%    plot_fig(15) = true;      %% (old) shift-distributions
%    plot_fig(16) = false;      %% (check) traces of PCs
%    plot_fig(17) = true;      %% ROI recurrence

%    plot_fig(18) = true;      %% population test
%    plot_fig(19) = true;      %% (shift to 10) PC/active numbers histogram
%    plot_fig(20) = true;      %% (not yet done) correlated population changes
%    plot_fig(21) = true;      %% firing rate / population
%    plot_fig(22) = true;      %% coding / session
%    plot_fig(23) = true;      %% firing rate correlation
%    plot_fig(24) = true;      %% population fraction reliability
%    plot_fig(25) = true;      %% number of fields
%    plot_fig(26) = true;      %% stability in populations (hard threshold)
%    plot_fig(27) = true;      %% population sizes
%    plot_fig(28) = true;      %% formation / deformation ratio
%    plot_fig(29) = true;      %% intercoding state
%    plot_fig(30) = true;      %% stability in populations (fit to distribution)
%    plot_fig(31) = true;      %% waiting time to next coding
  plot_fig(32) = true;      %% shift distributions dependence on time passed
%    plot_fig(33) = true;
%    plot_fig(34) = true;
%    plot_fig(:) = true;
  
%    if nargin < 8
    %% get neurons & PCs per session
    
    
%      PC_recurr = zeros(nSes,t_ses(end),3);  %% 1 = NRNG, 2 = GT, 3 = RW
%      PC_stable = zeros(nSes,t_ses(end),3);  %% 1 = NRNG, 2 = GT, 3 = RW
%      stable_pos = zeros(nSes,t_ses(end),para.nbin);
    
%      remap_pos = zeros(nSes,nSes-1,para.nbin+2,para.nbin+2);
    
%      remap_df = struct('hist',zeros(nSes-1,para.nbin,3),'ds',struct('rw',cell(nSes-1,1),'gt',cell(nSes-1,1),'nrng',cell(nSes-1,1)));
    
%      remap_df = zeros(nC,nSes-1,t_ses(end))*NaN;    %% session, delta session
    
    
    %%% ---------- post-process status: only take into account fields that pass the threshold ----------%%%
    
    %% remove fields, that dont pass the trial_frac-threshold
    [PCs,frac_active] = find_PCs(PCs,0.1,perc_thr,3);
    
    nROI = zeros(nSes,5); %% 1 = silent, 2 = active (nPC), 3 = NRNG, 4 = GT, 5 = RW
    nROI(:,1) = nansum(PCs.status(:,:,1),1);
    nROI(:,2) = sum(PCs.status(:,:,2)>0 & ~any(PCs.status(:,:,3:5),3),1);
    
    field_ct = nansum(PCs.status(:,:,3:5),3);
    nROI(:,3:5) = nansum(PCs.status(:,:,3:5)./field_ct,1);    %% count neurons with multiple fields as half (third, etc) only
    
%      nROI
    
    if nargin < 14 || isempty(arrays)
      
      ROI_recurr = zeros(nSes,t_ses(end));
      occupancy = zeros(nSes,para.nbin+2);
      
      N_pairs = zeros(t_ses(end),1);
      N_norm = zeros(nT,sum(t_mask),3);
      
      arrays = struct;
      arrays.ICI_s = {[],[],[]};
      arrays.ICI_t = {[],[],[]};
      
      arrays.s_ref = {[],[],[]};
      arrays.shift = {[],[],[]};
      arrays.r_coding = {[],[],[]};
      arrays.r_silent = {[],[],[]};
      
      arrays.ICI_discont = {[],[],[]};
      arrays.shift_discont = {[],[],[]};
      
      arrays.ICI_cont = {[],[],[]};
      arrays.shift_cont = {[],[],[]};
      
      tic
      for c = 1:nC
        
        if ROI_ct(c) >= ct_thr
          for s = s_start:nSes
            %% total number of ROIs active (and in cluster) per session
            
            if any(PCs.status(c,s,3:5))
              idx = round(PCs.fields.center(c,s,~isnan(PCs.fields.status(c,s,:))));
              occupancy(s,idx) = occupancy(s,idx) + 1/field_ct(c,s);
            elseif PCs.status(c,s,2)
              occupancy(s,end-1) = occupancy(s,end-1) + 1;
            else
              occupancy(s,end) = occupancy(s,end) + 1;
            end
            
            if any(PCs.status(c,s,3:5),3)
            
              for sm = s+1:nSes
                ds = sm-s;
                if time_real
                  dt = t_measures(sm) - t_measures(s);
                else
                  dt = ds;
                end
                for p = 1:3
                  if any(PCs.fields.status(c,s,:)==p+2,3)
                    
                    if time_real
                      idx_T = find(t_data_m==dt);
                    else
                      idx_T = find(t_data==dt);
                    end
                    N_norm(ds,idx_T,p) = N_norm(ds,idx_T,p)+1;  %% for each field, gather possible stable pairs
                
                    if any(PCs.status(c,sm,3:5),3)
                      
                      N_pairs(ds) = N_pairs(ds) + 1;
                      
                      if ds > 1
                        r_tmp = sum(any(PCs.status(c,s+1:sm-1,3:5),3))/(ds-1);
  %                        r_tmp = sum(any(PCs.status(c,s+1:sm-1,2:5),3))/(ds-1);
                        r_silent_tmp = sum(any(PCs.status(c,s+1:sm-1,1),3))/(ds-1);
                      else
                        r_tmp = 1;
                        r_silent_tmp = 1;
                      end
                      
                      field_ref = squeeze(PCs.fields.center(c,s,find(~isnan(PCs.fields.status(c,s,:)))));
                      field_1 = squeeze(PCs.fields.center(c,sm,find(~isnan(PCs.fields.status(c,sm,:)))));
                      
                      shift_tmp0 = mod(field_1'-field_ref+para.nbin/2,para.nbin)-para.nbin/2;
                      
                      if length(field_ref) == length(field_1)
                        [~,idx] = min(abs(shift_tmp0),[],2);
                        for j = 1:length(field_ref)
                          shift_tmp = shift_tmp0(j,idx(j));
                          
                          arrays.shift{p} = [arrays.shift{p},shift_tmp];
                          arrays.ICI_s{p} = [arrays.ICI_s{p},ds];
                          arrays.ICI_t{p} = [arrays.ICI_t{p},dt];
                          arrays.s_ref{p} = [arrays.s_ref{p},s];
                          arrays.r_coding{p} = [arrays.r_coding{p},r_tmp];
                          arrays.r_silent{p} = [arrays.r_silent{p},r_silent_tmp];
                        end
                      else
                        [~,idx] = min(abs(shift_tmp0(:)));
                        shift_tmp = shift_tmp0(idx);
                        
                        arrays.shift{p} = [arrays.shift{p},shift_tmp];
                        arrays.ICI_s{p} = [arrays.ICI_s{p},ds];
                        arrays.ICI_t{p} = [arrays.ICI_t{p},dt];
                        arrays.s_ref{p} = [arrays.s_ref{p},s];
                        arrays.r_coding{p} = [arrays.r_coding{p},r_tmp];
                        arrays.r_silent{p} = [arrays.r_silent{p},r_silent_tmp];
                      end
                    end
                  end
                end
              end
            end
            
              for sm = s+1:nSes
  %                ds_0 = sm-s;
                ds = t_ses(sm)-t_ses(s);
                if PCs.status(c,s,2) && PCs.status(c,sm,2)    %%% if active in session s and sm
  %                
                  ROI_recurr(s,ds) = ROI_recurr(s,ds) + 1;
                end
              end
  %                  if any(PCs.status(c,s,3:5)) && any(PCs.status(c,sm,3:5))    %%% if PC in session s and sm
  %                    
  %                    centers_tmp = PCs.fields.center(c,s,~isnan(PCs.fields.center(c,s,:)));
  %                    
  %                    PC_id = PCs.fields.status(c,s,~isnan(PCs.fields.status(c,s,:)))-2;
  %                    for ID = PC_id
  %                      PC_recurr(s,ds,ID) = PC_recurr(s,ds,ID) + 1;
  %                    end
  %                    
  %                    for fm = find(~isnan(PCs.fields.status(c,sm,:)))'
  %                    
  %                      pos = PCs.fields.center(c,sm,fm);
  %                      PC_id = PCs.fields.status(c,sm,fm)-2;
  %                      
  %                      disc = 1/field_ct(c,sm) * 1/field_ct(c,s);
  %                      
  %                      for f = find(~isnan(PCs.fields.status(c,s,:)))'
  %                          
  %                          pos_ref = PCs.fields.center(c,s,f);
  %                          remap_pos(s,ds_0,pos_ref,pos) = remap_pos(s,ds_0,pos_ref,pos) + disc;
  %                          
  %                      end
  %                      
  %                      if any(abs(mod(pos - centers_tmp + para.nbin/2,para.nbin)-para.nbin/2) < tolerance)
  %                        PC_stable(s,ds,PC_id) = PC_stable(s,ds,PC_id) + 1;
  %                        stable_pos(s,ds_0,pos) = stable_pos(s,ds_0,pos) + 1;
  %                      end
  %                    
  %                    end
  %                  elseif any(PCs.status(c,s,3:5))
  %  %                    disc = 1/field_ct(c,s);
  %  %                    for f = find(~isnan(PCs.fields.status(c,s,:)))'
  %                        pos_ref = PCs.fields.center(c,s,~isnan(PCs.fields.status(c,s,:)));
  %                        remap_pos(s,ds_0,pos_ref,para.nbin+1) = remap_pos(s,ds_0,pos_ref,para.nbin+1) + 1/field_ct(c,s);
  %  %                    end
  %                  elseif any(PCs.status(c,sm,3:5))
  %  %                    disc = 1/field_ct(c,sm);
  %  %                    for fm = find(~isnan(PCs.fields.status(c,sm,:)))'
  %                        pos = PCs.fields.center(c,sm,~isnan(PCs.fields.status(c,sm,:)));
  %                        remap_pos(s,ds_0,para.nbin+1,pos) = remap_pos(s,ds_0,para.nbin+1,pos) + 1/field_ct(c,sm);
  %  %                    end
  %                  else
  %                    remap_pos(s,ds_0,para.nbin+1,para.nbin+1) = remap_pos(s,ds_0,para.nbin+1,para.nbin+1) + 1;
  %                  end
  %                  
  %                elseif PCs.status(c,sm,2)
  %                  if any(PCs.status(c,sm,3:5))
  %  %                    disc = 1/field_ct(c,sm);
  %  %                    for fm = find(~isnan(PCs.fields.status(c,sm,:)))'
  %                        pos = PCs.fields.center(c,sm,~isnan(PCs.fields.status(c,sm,:)));
  %                        remap_pos(s,ds_0,para.nbin+2,pos) = remap_pos(s,ds_0,para.nbin+2,pos) + 1/field_ct(c,sm);
  %  %                    end
  %                  else
  %                    remap_pos(s,ds_0,para.nbin+2,para.nbin+1) = remap_pos(s,ds_0,para.nbin+2,para.nbin+1) + 1;
  %                  end
  %                elseif PCs.status(c,s,2)
  %                  if any(PCs.status(c,s,3:5))
  %  %                    disc = 1/field_ct(c,s);
  %  %                    for f = find(~isnan(PCs.fields.status(c,s,:)))'
  %                        pos_ref = PCs.fields.center(c,s,~isnan(PCs.fields.status(c,s,:)));
  %                        remap_pos(s,ds_0,pos_ref,para.nbin+2) = remap_pos(s,ds_0,pos_ref,para.nbin+2) + 1/field_ct(c,s);
  %  %                    end
  %                  else
  %                    remap_pos(s,ds_0,para.nbin+1,para.nbin+2) = remap_pos(s,ds_0,para.nbin+1,para.nbin+2) + 1;
  %                  end
  %                else
  %                  remap_pos(s,ds_0,para.nbin+2,para.nbin+2) = remap_pos(s,ds_0,para.nbin+2,para.nbin+2) + 1;
  %                end
  %              end
          end
        end
      end
    
      toc
        
  %      ROI_recurr(~ROI_recurr) = NaN;
  %      PC_recurr(~PC_recurr) = NaN;
  %      stable_pos(~stable_pos) = NaN;
  %      occupancy(~occupancy) = NaN;
  %      remap_pos(~remap_pos) = NaN;
      
  %      remap_df.hist(~remap_df.hist) = NaN;
    
    end
    
    ROI_total = zeros(nSes,t_ses(end),5);
    for s = 1:nSes
      for sm = s+1:nSes
        ds = t_ses(sm)-t_ses(s);
        for i = 1:5
          ROI_total(s,ds,i) = ROI_total(s,ds,i) + nROI(s,i);
        end
      end
    end
    ROI_total(~ROI_total) = NaN;
      
  
  
%%% ------------ get information about populations ------------- %%%
  
  pop_thr = 3;
  nOcc = squeeze(sum(PCs.status,2));
  nOcc(:,2) = nOcc(:,2) - squeeze(sum(any(PCs.status(:,:,3:5),3),2));
%    nOcc
  nPop = sum(nOcc>=pop_thr);
  
  idx_pure = sum(nOcc(:,3:5)>=pop_thr,2)==1;
  idx_none = sum(nOcc(:,3:5)>=pop_thr,2)==0;
  idx_mixed = sum(nOcc(:,3:5)>=pop_thr,2)>1;
  
  idx_pop = nOcc>=pop_thr;
%    idx_pop_pure = idx_pop & idx_pure;
  
  nPC_only = idx_none & idx_pop(:,2);
  NRNG_only = idx_pure & idx_pop(:,3);
  GT_only = idx_pure & idx_pop(:,4);
  RW_only = idx_pure & idx_pop(:,5);
  
  nPop_nPC = sum(nPC_only);
  nPop_GT = sum(GT_only);
  nPop_RW = sum(RW_only);
  nPop_NRNG = sum(NRNG_only);
  nPop_multi = sum(idx_mixed);
  
  pTotal = sum(nOcc) / (nC*nSes)
  
  
%%% ------------------------- plot basic cluster stats ----------------------- %%%
  
  if plot_fig(1)
    
    PC_number = nSes - sum(nOcc(:,1:2),2);
    
    h_edges = linspace(-0.5,nSes+0.5,nSes+2);
    n_edges = linspace(1,nSes,nSes);
    
    figure('position',[100 100 400 200])
    axes('position',[0.2,0.2,0.75,0.75])
%      subplot(2,1,1)
    
    hold on
%      histo_ct = bar(linspace(1,nSes,nSes),dat,'stacked');
%      y_lims = ylim;
%      plot([ct_thr ct_thr]+0.5,[0,y_lims(2)],'r--')
    
    histogram(ROI_ct,h_edges,'FaceColor','k','DisplayName','# sessions active');
    
    xlabel('# sessions','FontSize',14)
    ylabel('# neurons','FontSize',14)
    legend('Location','NorthEast')
    xlim([-0.5,nSes+0.5])
    ylim([0,250])
    if plt_presi
      path = pathcat(pathFigures,sprintf('Nactive%s.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
%      set(histo_ct(2),'YData',hist(ROI_ct,n_edges))
    
%      histo_ct(1).FaceColor = 'b';
%      histo_ct(2).FaceColor = 'r';
    
%      subplot(2,1,2)
    histogram(PC_number,h_edges,'FaceColor',[0.8,0.8,0.8],'DisplayName','# sessions coding')
    
    
%      ylim([0,1000])
%      xticks(linspace(0,nSes,4))
    
%      figure()
%      scatter(ROI_ct+0.5*rand(nC,1),PC_number+0.5*rand(nC,1),10,'.')
%      fit line through here
    
    
    if sv
      save_fig(pathFigures,mouse,'Nactive_NPC',sv_suffix,sv_ext,fig_pos)
%        path = pathcat(pathFigures,sprintf('matching_stats_ct_score%s.%s',sv_suffix,sv_ext));
%        path = pathcat(pathFigures,sprintf('Nactive_NPC%s.%s',sv_suffix,sv_ext));
%        print(path,sprintf('-d%s',sv_ext),'-r300')
%        disp(sprintf('Figure saved as %s',path))
    end
    
  end
  
  
%%% ------------------------------- plot matching results ------------------------------- %%%
  
  if plot_fig(2)
    
    
    fig_pos = [200 200 400 200];
    figure()
%      figure('position',[200 200 400 300])
    
%      ax2 = axes('position',[0.15 0.2 0.75 0.3]);
%  %      subplot(2,1,2)
%      active_time = zeros(nSes,1);
%      for s = 1:s_end
%        s_eff = s+s_offset-1;
%        pathSession = pathcat(pathMouse,sprintf('Session%02d',s_eff));
%        pathBH = dir(pathcat(pathSession,'*aligned.mat'));
%        pathBH = pathcat(pathSession,pathBH.name);
%        
%        bh = load(pathBH);
%        bh = bh.alignedData.resampled;
%        
%        active_time(s) = sum(bh.longrunperiod)/length(bh.longrunperiod);
%      end
%      
%      plot(ax2,t_measures(1:s_end),active_time,'k')
%      xlim(ax2,[0,t_measures(nSes)])
%      ylim(ax2,[0,1])
%      xlabel(ax2,'t [h]','FontSize',14)
%      ylabel(ax2,'active time')
    
    
    mask_plot = sum(nROI(:,2:5),2) > 500;
    ax1 = axes('position',[0.2 0.3 0.75 0.65]);
%      ax1 = subplot(1,1,1)
%      set(ax1,'position',[0.15 0.55 0.75 0.4],'XTick',[])
    hold(ax1,'on')
%      plot(linspace(1,nSes,nSes),nROI_raw,'k--')
%      plot(linspace(1,nSes,nSes),sum(nROI,2),'ko-')
    ylim([0,nC*1.5])
    xlabel(ax1,'session s','FontSize',12)
%      ylim([0,6000])
    
    ylabel('# neurons','FontSize',12)
%      legend('Location',[0.4 0.77 0.3 0.15])
%      if plt_presi
%        path = pathcat(pathFigures,sprintf('ROInum%s_1.png',sv_suffix));
%        print(path,'-dpng','-r300')
%      end
    plot(t_ses(mask_plot),ones(sum(mask_plot),1)*nC,'k:','DisplayName','# neurons')
    scatter(t_ses(mask_plot),sum(nROI(mask_plot,2:5),2),20,'ko','DisplayName','# active neurons')
%      ylim([0,nC*1.2])
    xlim([0,t_ses(nSes)])
    
%      scatter(t_ses,ones(nSes,1)*nC,'ko','DisplayName','# Neurons')
    legend('Location',[0.4 0.79 0.3 0.15])
    if plt_presi
      path = pathcat(pathFigures,sprintf('ROInum%s_2.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
%      scatter(t_measures(mask_plot),sum(nROI(mask_plot,3:5),2),10,'ko','filled','DisplayName','# place cells')
    legend('Location',[0.4 0.79 0.3 0.15])
    if plt_presi
      path = pathcat(pathFigures,sprintf('ROInum%s_3.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
%      plot(t_ses,nROI(:,3),'b-','DisplayName','# NRNG')
%      plot(t_ses,nROI(:,4),'r-','DisplayName','# GT')
%      plot(t_ses,nROI(:,5),'g-','DisplayName','# RW')
    lgd = legend('Location',[0.4 0.79 0.3 0.15]);
    lgd.FontSize = 10;
    if plt_presi
      path = pathcat(pathFigures,sprintf('ROInum%s_4.%s',sv_suffix,sv_ext));
      print(path,sprintf('-d%s',sv_ext),'-r300')
      disp(sprintf('Figure saved as %s',path))
    end
    hold off
    
%      set(gcf,'position',[200 200 300 250])
%      xlabel(ax1,'')
    
    
    if sv
      save_fig(pathFigures,mouse,'neuron_numbers',sv_suffix,sv_ext,fig_pos)
%        path = pathcat(pathFigures,sprintf('neuron_numbers%s.%s',sv_suffix,sv_ext));
%        print(path,sprintf('-d%s',sv_ext),'-r300')
%        disp(sprintf('Figure saved as %s',path))
    end
    
    
    idx_t = find(any(ROI_recurr,1));
    
    ROI_tot = squeeze(nansum(nansum(ROI_total(s_start:s_end,idx_t,2:5),3),1));
    
%      ROI_tot
    
    fig_pos = [600 200 400 200];
    fig = figure('Position',fig_pos);
    ax1 = axes('position',[0.2 0.3 0.75 0.65]);
%      subplot(2,1,1)
%      ROI_recurr
    ROI_recurr_0 = squeeze(nansum(ROI_recurr(s_start:s_end,idx_t),1));
    hold on
    plot(ax1,[0,180],[0.75,0.75],'k:','HandleVisibility','off')
    
%      h = plot(idx_t,ROI_recurr_0'./nansum(ROI_tot(:,2:5),2),'ko-','DisplayName','active neurons');
    scatter(ax1,idx_t,ROI_recurr_0./ROI_tot,20,'ko','filled','DisplayName','tracked neuron fraction');
%      set(h, {'MarkerFaceColor'}, get(h,'Color')); 
%      plot(idx_t,ROI_recurr_0(:,3)./ROI_tot(:,5),'go-','DisplayName','reward vPC');
%      plot(idx_t,ROI_recurr_0(:,2)./ROI_tot(:,4),'ro-','DisplayName','gate vPC')
%      plot(idx_t,ROI_recurr_0(:,1)./ROI_tot(:,3),'ko-','DisplayName','NRNG vPC')
    hold off

    y_lims = ylim;
    ylim(ax1,[0,1])
    xlim(ax1,[0,t_ses(end)])
    set(ax1,'YTick',linspace(0,1,5))
%      set(gca,'XTick',linspace(0,200,5))
%      yticks(linspace(0,1,11))
    xlabel(ax1,'session difference \Delta s','FontSize',12)
    ylabel(ax1,'fraction','FontSize',12)
    lgd = legend('Location',[0.3,0.35,0.35,0.13]);
    lgd.FontSize = 10;
%      set(gcf,'OuterPosition',[100 100 250 200],'InnerPosition',[100 100 250 200])
    
    if sv
      figure(fig)
      save_fig(pathFigures,mouse,'ROI_stability',sv_suffix,sv_ext,fig_pos)
%        path = pathcat(pathFigures,sprintf('ROI_stability%s.%s',sv_suffix,sv_ext));
%        print(path,sprintf('-d%s',sv_ext),'-r300')
%        disp(sprintf('Figure saved as %s',path))
    end
  
    
%      yl = ylim;
%      ylim([0,yl(2)])
    
    
    
    
%      path = pathcat(pathFigures,sprintf('/ROInum%s_1.png',sv_suffix));
%      print(path,'-dpng','-r300')
    
%      ax = subplot(1,2,2);
%      hold on
%      plot(ROI_recurr./ROI_total(:,2),'ko-')%,'DisplayName','total')
%  %      plot(sum(ROI_recurr,2)./ROI_total_raw,'k--')%,'DisplayName','total')
%      plot(ROI_recurr./ROI_total(:,2),'ro-')%,'DisplayName','processed part')
%  %      plot(ROI_recurr(:,2)./ROI_total(:,2),'bo-')%,'DisplayName','unprocessed part')
%      
%  %      plot(ROI_rec2(:,2)./ROI_tot2(:,2),'mo-')%,'DisplayName','raw')
%      
%      ROI_arr = ROI_recurr./ROI_total(:,2);
%      
%      for s = 1:3:nSes-1
%        text(s-0.5,ROI_arr(s)+0.01,sprintf('%4.2g %% (n = %d)',(100*ROI_arr(s)),ROI_total(s,2)))
%      end
%      
%      hold off
%  %      legend(ax,{'total','processed','unprocessed','raw'})
%      ylim([0.6,0.8])
%      yticks(linspace(0,1,6))
%      xlabel('\Delta s')
%      ylabel('prob. reoccurence')
%      xlim([0 t_ses(end)])
    
    
  end
  
%%% ----------------------------- plot PC stability ------------------------------------ %%%
%    PC_stable_0 = squeeze(nansum(PC_stable(s_start:s_end,:,:),1));
%    idx_t = find(any(PC_stable_0,2));
%    ROI_tot = squeeze(nansum(ROI_total(s_start:s_end,idx_t,:),1));
  if plot_fig(3)
    
    
    %% should be merged into the above graph
    figure('position',[200 200 400 250])
%      subplot(2,2,1)
%      hold on
%  %      plot(linspace(1,nSes,nSes),sum(nPC,2),'ko-')
%      plot(linspace(1,nSes,nSes),sum(nPC,2),'ko-')
%  %      plot(linspace(1,nSes,nSes),nPC(:,2),'bo-')
%      hold off
%      yl = ylim;
%      ylim([0,yl(2)])
%      xlabel('session')
%      ylabel('# PCs')
    
    
%      subplot(2,2,2)
    hold on
%      sum(squeeze(sum(PC_stable(s_start:s_end,:,:),1)),2)
%      sum(squeeze(sum(PC_stable(s_start:s_end,:,:),1)),2)./sum(PC_total,2)
    
    PC_stable_0 = PC_stable_0(idx_t,:);
    
%      size(PC_stable_0)
    ROI_tot = squeeze(nansum(ROI_total(s_start:s_end,idx_t,:),1));
    
    plot(idx_t,nansum(PC_stable_0,2)./nansum(ROI_tot(:,3:6),2),'k:','DisplayName','all vPC')
    plot(idx_t,PC_stable_0(:,3)./ROI_tot(:,6),'go-','DisplayName','reward vPC');
    plot(idx_t,PC_stable_0(:,2)./ROI_tot(:,5),'ro-','DisplayName','wall vPC')
    plot(idx_t,PC_stable_0(:,2)./ROI_tot(:,4),'ro-','DisplayName','gate vPC')
    plot(idx_t,PC_stable_0(:,1)./ROI_tot(:,3),'ko-','DisplayName','NRNG vPC')
    
    hold off
    y_lims = ylim;
    ylim([0,0.5])
    xlim([0,180])
    yticks(linspace(0,1,11))
    xlabel('\Delta time [h]','Fontsize',14)
    ylabel('stability','Fontsize',14)
    legend()
    
    
%      subplot(2,2,4)
%      hold on
%      plot(idx_t,nansum(ROI_tot(:,3:5),2),'k:','DisplayName','all vPC')
%      plot(idx_t,ROI_tot(:,5),'go-','DisplayName','reward vPC');
%      plot(idx_t,ROI_tot(:,4),'ro-','DisplayName','gate vPC')
%      plot(idx_t,ROI_tot(:,3),'ko-','DisplayName','NRNG vPC')
%      
%      hold off
%      y_lims = ylim;
%      xlim([0,180])
%      xlabel('\Delta time [h]','Fontsize',14)
%      ylabel('# of possible pairs','Fontsize',14)
%      set(gca,'yscale','log')
%      
%      subplot(2,2,1)
%      PC_recurr_0 = squeeze(nansum(PC_recurr(s_start:s_end,idx_t,:),1));
%  %      PC_recurr_0
%      hold on
%      plot(idx_t,sum(PC_recurr_0,2)./nansum(ROI_tot(:,3:5),2),'k:','DisplayName','all vPC')
%      plot(idx_t,PC_recurr_0(:,3)./ROI_tot(:,5),'go-','DisplayName','reward vPC');
%      plot(idx_t,PC_recurr_0(:,2)./ROI_tot(:,4),'ro-','DisplayName','gate vPC')
%      plot(idx_t,PC_recurr_0(:,1)./ROI_tot(:,3),'ko-','DisplayName','NRNG vPC')
%      hold off
%  
%      y_lims = ylim;
%      ylim([0,1])
%      xlim([0,180])
%      yticks(linspace(0,1,11))
%      xlabel('\Delta time [h]')
%      ylabel('prob. reoccurence (vPC)')
    
    if sv
      path = pathcat(pathFigures,sprintf('PC_stability%s.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
    
  end
  
  
  if plot_fig(4)
    
    %% # of fields
    mean_fc = zeros(nSes,1);
    for s = 1:nSes
      mean_fc(s) = mean(field_ct(field_ct(:,s)>0,s));
    end
    figure('position',[500 400 600 400])
    plot(mean_fc,'ko-')
    ylim([0,2])
    xlabel('Session')
    ylabel('avg. #fields')
    
    if sv
      path = pathcat(pathFigures,sprintf('field_num%s.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
    
  end
  
  if plot_fig(5)
    
    
    early_idx = 1:3;
    mid_idx = 4:10;
    end_idx = 11:15;
    %%% stable positions
%      stable_pos_proc = squeeze(stable_pos(:,1,:));
%      stable_pos_raw = squeeze(stable_pos(:,1,:,2));
    stable_pos_total = squeeze(stable_pos(:,1,:));
    
%      figure('position',[100 100 1200 300])
%      subplot(2,2,1)
%      bar(sum(stable_pos_proc,1))
    
%      subplot(2,2,3)
%      bar(sum(stable_pos_raw,1))
    
%      subplot(2,2,2)
%      bar(sum(stable_pos_proc,2))
%      
%      subplot(2,2,4)
%      bar(sum(sum(stable_pos(:,:,:,1),1),3))
    
    figure('position',[500,100,400,300])
    
    ax1 = subplot(3,1,1);
%      axes('position',[0.1 0.7])
    hold on
    bar(gt_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none')
%      bar(wc_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none')
    bar(rw_bar,1,'FaceColor',[0.8,1,0.8],'EdgeColor','none')
    bar(PC_bar,1,'FaceColor',[0.7,0.7,1],'EdgeColor','none')
    
    bar(nansum(stable_pos_total(early_idx,:),1)/sum(nansum(stable_pos_total(early_idx,:))),'k','EdgeColor','none')
%      bar(sum(occupancy(1:3,1:end-2)/sum(sum(occupancy(1:3,1:end-2))),1),0.9,'k','EdgeColor','none')
    
    ylims = [0,0.15]
    title('early (1-3)')
    ylim(ylims)
    xticks([])
    set(gca,'YTick',linspace(0,0.1,3),'YTickLabel',linspace(0,10,3))
    ylabel('% of PC')
    
    ax2 = subplot(3,1,2);
    hold on
    bar(gt_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none')
    bar(rw_bar,1,'FaceColor',[0.8,1,0.8],'EdgeColor','none')
    bar(PC_bar,1,'FaceColor',[0.7,0.7,1],'EdgeColor','none')
    bar(nansum(stable_pos_total(mid_idx,:),1)/sum(nansum(stable_pos_total(mid_idx,:))),'k','EdgeColor','none')
%      bar(sum(occupancy(4:10,1:end-2)/sum(sum(occupancy(4:10,1:end-2))),1),0.9,'k','EdgeColor','none')
    title('middle (4-10)')
    ylim(ylims)
    xticks([])
    set(gca,'YTick',linspace(0,0.1,3),'YTickLabel',linspace(0,10,3))
    ylabel('% of PC')
    
    ax3 = subplot(3,1,3);
    hold on
    bar(gt_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none')
    bar(rw_bar,1,'FaceColor',[0.8,1,0.8],'EdgeColor','none')
    bar(PC_bar,1,'FaceColor',[0.7,0.7,1],'EdgeColor','none')
    bar(nansum(stable_pos_total(end_idx,:),1)/sum(nansum(stable_pos_total(end_idx,:))),'k','EdgeColor','none')
%      bar(sum(occupancy(11:15,1:end-2)/sum(sum(occupancy(11:15,1:end-2))),1),0.9,'k','EdgeColor','none')
    title('late (11-15)')
    set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5))
%      xlim([0,81])
    xlabel('Location [cm]')
    ylim(ylims)
    set(gca,'YTick',linspace(0,0.1,3),'YTickLabel',linspace(0,10,3))
    ylabel('% of PC')
    
    
%      subplot(1,3,1)
%      
%      title('early')
%      ylim([0,40])
%      set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5))
%      
%      subplot(1,3,2)
%      bar(nansum(stable_pos_total(mid_idx,:),1))
%      title('middle')
%      ylim([0,40])
%      set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5))
%      
%      subplot(1,3,3)
%      bar(nansum(stable_pos_total(end_idx,:),1),'k','EdgeColor','none')
%      title('late')
%      ylim([0,40])
%      set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5))
    
    if sv
      path = pathcat(pathFigures,sprintf('PC_stability%s.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
  end
  
  
  
  if plot_fig(6)
    
    
%      ds_arr = 1:10
    ds_arr = 1;
    for ds = ds_arr
      remap_ds = squeeze(remap_pos(:,ds,:,:));
      
%        remap_PC = squeeze(sum(remap_ds(:,para.idx_zone{3},:),2));
%        remap_gt = squeeze(sum(remap_ds(:,para.idx_zone{4},:),2));
%        remap_rw = squeeze(sum(remap_ds(:,para.idx_zone{5},:),2));
      
%        csvFile = pathcat(pathMouse,sprintf('ROI_transitions_gt.csv'));
%        dlmwrite(csvFile,remap_gt,'delimiter',';')
%        
%        csvFile = pathcat(pathMouse,sprintf('ROI_transitions_rw.csv'));
%        dlmwrite(csvFile,remap_rw,'delimiter',';')
%        
%        csvFile = pathcat(pathMouse,sprintf('ROI_transitions_PC.csv'));
%        dlmwrite(csvFile,remap_PC,'delimiter',';')
%        
%        x_margin = 0.03;
%        y_margin = 0.0;
      
      
%        figure('position',[500 100 1200 900])
      
      figure('position',[500 100 1200 300])
      
      ylims = [0.005,0.01,0.02,0.06,0.06,0.06];
      
      bar_zone{3} = PC_bar*ylims(3);
      bar_zone{4} = gt_bar*ylims(4);
%        bar_zone{5} = wc_bar*ylims(5);
      bar_zone{5} = rw_bar*ylims(6);
      
      ylabels = {'silent','nPC','PC_{NRNG}','PC_{test}','PC_{reward}'};
      
      parts = 3;
      for j = 1:parts
        if j == 1
          idx_state = 1:3;
          title_str = 'early (1-3)';
        elseif j == 2
          idx_state = 4:10;
          title_str = 'middle (4-10)';
        else
          idx_state = 11:15;
          title_str = 'late (11-15)';
        end
      
%        for j = 1:parts
%          idx_state = (2*j-1):min(2*j,nSes);
%          title_str = '';
        idx_state
        idx_state_in = max(1,(idx_state(1)-1)-(ds-1)):((idx_state(end)-1)-(ds-1));
        idx_state_out = idx_state(1):min(idx_state(end),nSes-ds);
        idx_state_norm_in = max(idx_state(1),1+ds):idx_state(end);
        idx_state_norm_out = idx_state(1):min(idx_state(end),(nSes-ds));
        
  %        remap_hist_tot = squeeze(sum(remap_pos(idx_state,1,:,:),1));
        
  %        occupancy_phase = sum(occupancy(idx_state,:),1);
        remap_ds_in = squeeze(sum(remap_ds(idx_state_in,1:para.nbin,:),1));
        remap_ds_out = squeeze(sum(remap_ds(idx_state_out,:,1:para.nbin),1));
        
        for i = 5%1:6
          any_in(i,:) = sum(remap_ds_in(:,para.idx_zone{i}),2)/sum(sum(occupancy(idx_state_norm_in,para.idx_zone{i})));
          any_out(i,:) = sum(remap_ds_out(para.idx_zone{i},:),1)/sum(sum(occupancy(idx_state_norm_out,para.idx_zone{i})));
        end
        
        
        for i = 5%1:6
%            subplot(5,parts,(i-1)*parts+1+(j-1))
          subplot(1,3,j)
          hold on
          if i >= 3
            bar(bar_zone{i},1,'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
            bar(-bar_zone{i},1,'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
          end
          bar(any_in(i,:),'g','EdgeColor','none')
          bar(-any_out(i,:),'r','EdgeColor','none')
          bar(any_in(i,:)-any_out(i,:),'b','EdgeColor','none')
          ylim([-ylims(i),ylims(i)])
          set(gca,'XTick',[])
          
          title(title_str,'FontSize',10)
          if j == 1
            ylabel(ylabels{i},'FontSize',16)
            
          end
          if j == 2 && i == 1
            text(50,3/4*ylims(1),'recruiting from')
            text(50,-3/4*ylims(1),'get recruited to')
          end
          
          if i >= 4
            text(para.idx_zone{i}(end-1)-25,3/4*ylims(i),sprintf('%4.2g %%',(100*sum(any_in(i,para.idx_zone{i})))),'FontSize',14)
            text(para.idx_zone{i}(end-1)-25,-3/4*ylims(i),sprintf('%4.2g %%',(100*sum(any_out(i,para.idx_zone{i})))),'FontSize',14)
          end
%            if i == 5
            set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5),'FontSize',14)
            xlabel('Location [cm]','FontSize',14)
%            end
        end
        
      end
      if sv
        if plt_presi
          path = pathcat(pathFigures,sprintf('stable/PC_recruitment%d%s.png',ds,sv_suffix));
        else
          path = pathcat(pathFigures,sprintf('PC_recruitment%s_test.png',sv_suffix));
        end
        print(path,'-dpng','-r300')
      end
      
    end
  end
  
  
  if plot_fig(7)
    
    perc = linspace(0.1,1,10);
    nPC_perc = zeros(length(perc),nSes);
    
    for c = 1:nC
      if ROI_ct(c) >= ct_thr
        for s = 1:nSes
          if any(PCs.status(:,s,3:5))
            for i = 1:length(perc)
              for f = 1:sum(~isnan(PCs.fields.center(c,s,:)))
                if frac_active(c,s,f) >= perc(i)
                  nPC_perc(i,s) = nPC_perc(i,s) + 1;
                  break
                end
              end
            end
          end
        end
      end
    end
    
    nPC_perc_norm = nPC_perc./nPC_perc(1,:);%sum(nROI(:,3:5),2)'
    
    figure('position',[500 100 400 250])
    hold on
    plot([perc_thr,perc_thr],[0,1],'r--')
    for i = 1:length(perc)
      scatter(perc(i)*ones(nSes,1),nPC_perc_norm(i,:),10,[0.8, 0.8, 0.8],'o','filled')
      nPC_mean = mean(nPC_perc_norm(i,:));
      nPC_var = var(nPC_perc_norm(i,:));
      errorbar(perc(i),nPC_mean,nPC_var,'kx','MarkerSize',5)
    end
    hold off
    xlabel('active trial threshold [%]')
    set(gca,'XTick',linspace(0,1,6),'XTickLabels',linspace(0,100,6))
    ylabel('vPCs [fraction]')
    ylim([0,1])
    xlim([-0.1,1.1])
    
%      figure('position',[500 100 600 400])
%      hold on
%      
%      idx = 5;
%      for i = 1:length(perc)
%        frac = nPC_perc(i,:)./sum(nROI(:,3:5),2)';
%        plot(frac,'Color',[i i i]./double(length(perc)))
%        if ~mod(i+1,2)
%  %          [max_val, max_pos] = max(frac);
%          idx = idx + 1;
%          text(idx,frac(idx)+0.02,sprintf('%02d %%',round(100*perc(i))))
%        end
%      end
%        
%      hold off
%      xlabel('Session')
%      ylabel('fraction vPCs')
%      xlim([0.5,15.5])
    if sv
      path = pathcat(pathFigures,sprintf('field_thr%s.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
  end
  
  
  if plot_fig(8)
    
    pop_mean = mean(PCs.MI.rand_distr,3);
    pop_var = std(PCs.MI.rand_distr,0,3);
    
    zs = (PCs.MI.value-pop_mean)./pop_var;
    
    MI_mean = nanmean(PCs.MI.value,1);
    MI_var = sqrt(nanvar(PCs.MI.value,0,1));
    
    
    figure('position',[500 100 900 600])
    subplot(2,1,1)
    errorbar(linspace(1,nSes,nSes),MI_mean,MI_var)
    xlim([0,16])
%      ylim([0,5])
    ylabel('MI')
    
    subplot(2,1,2)
    errorbar(linspace(1,nSes,nSes),nanmean(zs,1),sqrt(nanvar(zs,0,1)))
    
    xlim([0,16])
%      ylim([1,2])
    xlabel('Session')
    ylabel('MI fraction of 95% shuffle')
    
    
%      figure
%      for s = 1:nSes
%        subplot(3,5,s)
%        hold on
%        histogram(zs(:,s),linspace(0,30,31))
%        histogram(zs(any(PCs.status(:,s,3:5),3),s),linspace(0,30,31),'FaceColor','r')
%      end
    
    zs_hist = histc(zs(:),linspace(0,20,51));
    zs_PC_hist = histc(zs(any(PCs.status(:,:,3:5),3)),linspace(0,20,51));
    
    MI_hist = histc(PCs.MI.value(:),linspace(0,0.3,51));
    MI_PC_hist = histc(PCs.MI.value(any(PCs.status(:,:,3:5),3)),linspace(0,0.3,51));
    
    figure
    subplot(2,2,1)
    hold on
    histogram(zs(:),linspace(0,20,51))
    histogram(zs(any(PCs.status(:,:,3:5),3)),linspace(0,20,51),'FaceColor','r')
    xlim([0,20])
    
    subplot(2,2,3)
    plot(linspace(0,20,51),zs_PC_hist./zs_hist,'r')
    xlim([0,20])
    ylim([0,1])
    
    subplot(2,2,2)
    hold on
    histogram(PCs.MI.value(:),linspace(0,0.3,51))
    histogram(PCs.MI.value(any(PCs.status(:,:,3:5),3)),linspace(0,0.3,51),'FaceColor','r')
    xlim([0,0.3])
    
    subplot(2,2,4)
    plot(linspace(0,0.3,51),MI_PC_hist./MI_hist,'r')
    xlim([0,0.3])
    ylim([0,1])
    
    
    MI_stats = zeros(nSes,5,2);
    for s = 1:nSes
      p=1; %% all
      MI_stats(s,p,1) = nanmean(zs(:,s));
      MI_stats(s,p,2) = sqrt(nanvar(zs(:,s)));
      
      p=2; %% all
      MI_stats(s,p,1) = nanmean(zs(any(PCs.status(:,s,3:5),3),s));
      MI_stats(s,p,2) = sqrt(nanvar(zs(any(PCs.status(:,s,3:5),3),s)));
      
      for p=3:5
        MI_stats(s,p,1) = nanmean(zs(PCs.status(:,s,p),s));
        MI_stats(s,p,2) = sqrt(nanvar(zs(PCs.status(:,s,p),s)));
      end
    end
    
    figure
    hold on
    
    col_arr = {'k','k:','b','r','m','g'};
    for p=1:5
      errorbar(linspace(1,nSes,nSes),MI_stats(:,p,1),MI_stats(:,p,2),col_arr{p})
    end
    
    
    if sv
      path = pathcat(pathFigures,sprintf('MI%s.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
  end
  
  
  if plot_fig(9)
    
    %% get centers of all coding cells in session s
    s_ref = 15;
    n_plots = 7;
    ordered = true;
    
    if ordered
      centers = zeros(0,2);
      i=1;
      for c=1:nC
        if any(PCs.status(c,s_ref,3:5))
          centers(i,1) = nanmin(PCs.fields.center(c,s_ref,:));
          centers(i,2) = c;
          i=i+1;
        end
      end
      centers = sortrows(centers,1);
      nID = size(centers,1);
    end
    
%      figure('position',[200 200 1800 300])
%      for s=1:nSes
    for ds = -(n_plots-1)/2:(n_plots-1)/2
%        ScrSz = get(0,'ScreenSize')
%        f = figure();%'position',[200 200 200 300])
%        set(0,'DefaultAxesColor','none')
        s = s_ref+ds
        if ~ordered
          centers = zeros(0,2);
          i=1;
          for c=1:nC
            if any(PCs.status(c,s,3:5))
              centers(i,1) = nanmin(PCs.fields.center(c,s,:));
              centers(i,2) = c;
              i=i+1;
            end
          end
          centers = sortrows(centers,1);
          nID = size(centers,1);
        end
      
%        subplot(3,5,s)
      
        f = figure();
        ds = 0;
        firingstacks = zeros(nID,para.nbin);
        i=1;
        for c=centers(:,2)'
          firingstacks(i,:) = PCs.fields.firingmap(c,s+ds,:) - min(PCs.fields.firingmap(c,s+ds,:));
          firingstacks(i,:) = firingstacks(i,:)./max(firingstacks(i,:));   % normalize to maximum
          
          i = i+1;
        end
%          ax = subplot(1,n_plots+1,ds+(n_plots-1)/2+1);
        ax = axes('Position',[0.25,0.15,0.65,0.65]);
        hold on
        imagesc(ax,firingstacks)
        colormap('jet')
        ax.YDir = 'reverse';
        
%          if ds == -(n_plots-1)/2
%            plot([gt_zone(1) gt_zone(1)],[1,nID],'r:','LineWidth',3)
%            plot([gt_zone(2) gt_zone(2)],[1,nID],'r:','LineWidth',3)
%            
%            plot([rw_zone(1) rw_zone(1)],[1,nID],'g:','LineWidth',3)
%            plot([rw_zone(2) rw_zone(2)],[1,nID],'g:','LineWidth',3)
%            
%            ylabel('Neuron ID')
%          end
        
        xticks(ax,linspace(1,para.nbin,5))
        xticklabels(ax,linspace(0,100,5))
        xlim(ax,[1,para.nbin])
        ylim(ax,[1,nID])
        xlabel(ax,'Location [cm]')
        
        ax_hist = axes('Position',[0.25,0.8,0.65,0.15]);
        bar(ax_hist,linspace(1,para.nbin,para.nbin),occupancy(s,1:end-2),'FaceColor','k')
        set(ax_hist,'Xtick',[],'YTick',[],'visible','off')
        suptitle(sprintf('s=%d',s))
%          set()
        
        set(f,'OuterPosition',[200 200 200 300],'InnerPosition',[200 210 200 290])
        
%          set(gcf,'Color','None')
%          set(ax,'Color','None')
%          set(ax_hist,'Color','None')
        
        
%          histcounts(centers,linspace(1,para.nbin,para.nbin));
%        end
        
%        ax = subplot(1,n_plots+1,n_plots+1);
%        colormap('jet')
%        cb = colorbar();%'Location','West');
%        ylabel(cb,'normalized Ca^{2+} activity above baseline')
%        set(ax,'visible','off')
        
        
        if sv
          if ordered
            path = pathcat(pathFigures,sprintf('PC_coverage_session_%d_aligned_%d%s.%s',s,s_ref,sv_suffix,sv_ext));
          else
            path = pathcat(pathFigures,sprintf('PC_coverage_session_%d_nonaligned_%d%s.%s',s,s_ref,sv_suffix,sv_ext));
          end
          export_fig(path,sprintf('-%s',sv_ext),'-transparent','-r300')
%            print(path,sprintf('-d%s',sv_ext),'-r300')
          disp(sprintf('Figure saved as %s',path))
        end
      end

%      if sv
%        path = pathcat(pathFigures,sprintf('PC_coverage_multiple%s.png',sv_suffix));
%        print(path,sprintf('-d%s',sv_ext),'-r300')
%        disp(sprintf('Figure saved as %s',path))
%      end
    
    plt_single = true;
    if plt_single
      figure('position',[500,100,400,150])
      ax1 = axes('position',[0.15,0.3,0.75,0.65]);
  %      axes('position',[0.1 0.7])
      hold on
      bar(gt_bar,1,'FaceColor',[0.8,1,0.8],'EdgeColor','none')
      bar(rw_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none')
      bar(PC_bar,1,'FaceColor',[0.7,0.7,1],'EdgeColor','none')
      
      bar(sum(occupancy(4:s_end,1:end-2)/sum(sum(occupancy(4:s_end,1:end-2))),1),0.9,'k','EdgeColor','none')
      
      ylims = [-0.02,0.06];
%        title('early (1-3)')
      ylim(ylims)
%        xticks([])
      set(gca,'YTick',linspace(0,0.1,3),'YTickLabel',linspace(0,10,3))
      ylabel('% of PC','FontSize',14)
      
      axes(ax1)
      yyaxis right
      hold on
      plot(sum(occupancy(4:s_end,1:para.nbin))/(sum(sum(nROI(4:s_end,3:5)))/para.nbin),'r:')
      plot(imgaussfilt(sum(occupancy(4:s_end,1:para.nbin))/(sum(sum(nROI(4:s_end,3:5)))/para.nbin),2,'Padding','circular'),'r','LineWidth',2)
      ylim([0,4])
      ylabel('ratio','FontSize',14)
      set(gca,'ycolor','r')
      xlabel('Position [cm]','FontSize',14)
      set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5))
      
      yyaxis left
%        text(19,0.055,'GT','FontSize',14)
%        text(67,0.055,'RW','FontSize',14)
      
      if sv
        path = pathcat(pathFigures,sprintf('PC_coverage_single%s.%s',sv_suffix,sv_ext));
        print(path,sprintf('-d%s',sv_ext),'-r300')
        disp(sprintf('Figure saved as %s',path))
      end
    else
    
    
    
      figure('position',[500,100,400,300])
      
      ax1 = subplot(3,1,1);
  %      axes('position',[0.1 0.7])
      hold on
      bar(gt_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none')
      bar(rw_bar,1,'FaceColor',[0.8,1,0.8],'EdgeColor','none')
      bar(PC_bar,1,'FaceColor',[0.7,0.7,1],'EdgeColor','none')
      
      bar(sum(occupancy(1:3,1:end-2)/sum(sum(occupancy(1:3,1:end-2))),1),0.9,'k','EdgeColor','none')
      
      ylims = [-0.035,0.07]
      title('early (1-3)')
      ylim(ylims)
      xticks([])
      set(gca,'YTick',linspace(0,0.1,3),'YTickLabel',linspace(0,10,3))
      ylabel('% of PC')
      
      ax2 = subplot(3,1,2);
      hold on
      bar(gt_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none')
      bar(rw_bar,1,'FaceColor',[0.8,1,0.8],'EdgeColor','none')
      bar(PC_bar,1,'FaceColor',[0.7,0.7,1],'EdgeColor','none')
      bar(sum(occupancy(4:10,1:end-2)/sum(sum(occupancy(4:10,1:end-2))),1),0.9,'k','EdgeColor','none')
      title('middle (4-10)')
      ylim(ylims)
      xticks([])
      set(gca,'YTick',linspace(0,0.1,3),'YTickLabel',linspace(0,10,3))
      ylabel('% of PC')
      
      ax3 = subplot(3,1,3);
      hold on
      bar(gt_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none')
      bar(rw_bar,1,'FaceColor',[0.8,1,0.8],'EdgeColor','none')
      bar(PC_bar,1,'FaceColor',[0.7,0.7,1],'EdgeColor','none')
      bar(sum(occupancy(11:15,1:end-2)/sum(sum(occupancy(11:15,1:end-2))),1),0.9,'k','EdgeColor','none')
      title('late (11-15)')
      set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5))
  %      xlim([0,81])
      xlabel('Location [cm]')
      ylim(ylims)
      set(gca,'YTick',linspace(0,0.1,3),'YTickLabel',linspace(0,10,3))
      ylabel('% of PC')
      
      if plt_presi
        path = pathcat(pathFigures,sprintf('PC_coverage_presi%s.png',sv_suffix));
        print(path,'-dpng','-r300')
      end
      
      axes(ax1)
      yyaxis right
      hold on
      plot(sum(occupancy(1:3,1:para.nbin))/(sum(sum(nROI(1:3,3:6)))/para.nbin),'r:')
      plot(imgaussfilt(sum(occupancy(1:3,1:para.nbin))/(sum(sum(nROI(1:3,3:6)))/para.nbin),2,'Padding','circular'),'r','LineWidth',2)
      ylim([0,3])
      ylabel('ratio')
      set(gca,'ycolor','r')
      
      axes(ax2)
      yyaxis right
      hold on
      plot(sum(occupancy(4:10,1:para.nbin))/(sum(sum(nROI(4:10,3:6)))/para.nbin),'r:')
      plot(imgaussfilt(sum(occupancy(4:10,1:para.nbin))/(sum(sum(nROI(4:10,3:6)))/para.nbin),2,'Padding','circular'),'r','LineWidth',2)
      ylim([0,3])
      ylabel('ratio')
      set(gca,'ycolor','r')
      
      axes(ax3)
      yyaxis right
      plot(sum(occupancy(11:15,1:para.nbin))/(sum(sum(nROI(11:15,3:6)))/para.nbin),'r:')
      plot(imgaussfilt(sum(occupancy(11:15,1:para.nbin))/(sum(sum(nROI(11:15,3:6)))/para.nbin),2,'Padding','circular'),'r','LineWidth',2)
      ylim([0,3])
      ylabel('ratio')
      set(gca,'ycolor','r')
      
      if sv
        path = pathcat(pathFigures,sprintf('PC_coverage%s.png',sv_suffix));
        print(path,'-dpng','-r300')
        disp(sprintf('Figure saved as %s',path))
      end
    end
%      overrepr = occupancy(:,1:para.nbin)./(sum(nROI(:,3:5),2)/para.nbin);
    
    
    
    
    
  end
  
  
  
  if plot_fig(10)
    
%      s = ;
    for s = 10
      figure('position',[500 500 400 300])
%        suptitle(sprintf('session %d',s))
      
      ax1 = axes('position',[0.15 0.15 0.65 0.8]);
      
      set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5))
      set(gca,'YTick',linspace(-10,10,5),'YTickLabel',[10,5,0,5,10])
      xlabel('Location [cm]')
      ylim([-12,12])
      xticks([])
%        set(gca,'YTick',linspace(0,0.1,3),'YTickLabel',linspace(0,10,3))
      ylabel('# PC')
      hold on
      
      ax2 = axes('position',[0.82 0.15 0.1 0.8]);
      set(gca,'XTick',[1,2],'XTickLabel',{'nPC','silent'},'YAxisLocation','right')
      set(gca,'YTick',linspace(-50,50,3),'YTickLabel',[50,0,50])
      xtickangle(60)
      ylim([-50 50])
      xlim([0,3])
      hold on
      
      axes(ax1)
      h_bar_gt_u = bar(12*gt_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none');
      h_bar_gt_d = bar(-12*gt_bar,1,'FaceColor',[1,0.8,0.8],'EdgeColor','none');
      h_bar_occ = bar(occupancy(s,1:para.nbin),0.9,'k','EdgeColor','none');
%        bar(ax2,occupancy(s,para.nbin+1:para.nbin+2),0.5,'k','EdgeColor','none');
%        waitforbuttonpress
      if plt_presi
        path = pathcat(pathFigures,sprintf('Presi/remapping1%s.png',sv_suffix));
        print(path,'-dpng','-r300')
      end
        
      delete(h_bar_occ)
      h_bar_occ_gt = bar(para.idx_zone{4},occupancy(s,para.idx_zone{4}),0.9,'k','EdgeColor','none');
      set(h_bar_gt_u,'FaceColor',[0.9 0.9 0.9])
      set(h_bar_gt_d,'FaceColor',[0.9 0.9 0.9])
      remap_ds = squeeze(remap_pos(:,1,:,:));
      remap_gt_pre = sum(remap_ds(s-1,:,para.idx_zone{4}),3);
      remap_gt_post = squeeze(sum(remap_ds(s,para.idx_zone{4},:),2));
%        waitforbuttonpress
      if plt_presi
        path = pathcat(pathFigures,sprintf('Presi/remapping2%s.png',sv_suffix));
        print(path,'-dpng','-r300')
      end
      
      bar(remap_gt_pre(1:para.nbin),0.9,'g','EdgeColor','none');
      bar(ax2,remap_gt_pre(para.nbin+1:para.nbin+2),0.5,'g','EdgeColor','none');
      text(40,8,'recruited from')
%        waitforbuttonpress
      if plt_presi
        path = pathcat(pathFigures,sprintf('Presi/remapping3%s.png',sv_suffix));
        print(path,'-dpng','-r300')
      end
      
      bar(-remap_gt_post(1:para.nbin,:),0.9,'r','EdgeColor','none');
      bar(ax2,-remap_gt_post(para.nbin+1:para.nbin+2),0.5,'r','EdgeColor','none');
      text(40,-10,'recruited to')
%        waitforbuttonpress
      if plt_presi
        path = pathcat(pathFigures,sprintf('Presi/remapping4%s.png',sv_suffix));
        print(path,'-dpng','-r300')
      end
      
      bar(remap_gt_pre(1:para.nbin)-remap_gt_post(1:para.nbin)',0.9,'b','EdgeColor','none')
      bar(ax2,remap_gt_pre(para.nbin+1:para.nbin+2)-remap_gt_post(para.nbin+1:para.nbin+2)',0.5,'b','EdgeColor','none');
%        waitforbuttonpress
      if plt_presi
        path = pathcat(pathFigures,sprintf('Presi/remapping5%s.png',sv_suffix));
        print(path,'-dpng','-r300')
      end
      
      delete(h_bar_occ_gt)
      if plt_presi
        path = pathcat(pathFigures,sprintf('Presi/remapping6%s.png',sv_suffix));
        print(path,'-dpng','-r300')
      end
      
    end
    
  end
  
  overlap_calc = false;
  if ~exist('pop_overlap','var')
    disp('pop_overlap non-existant')
    pop_overlap = struct;
    overlap_calc = true;
  end
  
  
  
  if plot_fig(11)
    %% plot # of neurons active in session x and y
    
    if overlap_calc
      pop_overlap.overrepr_act = zeros(nSes,nSes);
      pop_overlap.frac_act = zeros(nSes,nSes);
      
      pop_overlap.overrepr_PC = zeros(nSes,nSes);
      pop_overlap.frac_PC = zeros(nSes,nSes);
      
      L = 100;
      for s = 1:s_end
        
        mask_act = PCs.status(:,s,2);
        mask_PC = any(PCs.status(:,s,3:5),3);
        disp(sprintf('s=%d, active neurons: %d, PCs: %d',s,sum(mask_act),sum(mask_PC)))
        for sm = 1:s_end
          
          if s==sm
            pop_overlap.overrepr_act(s,sm) = 0;
            pop_overlap.frac_act(s,sm) = 1;
            
            pop_overlap.overrepr_PC(s,sm) = 0;
            pop_overlap.frac_PC(s,sm) = NaN;
          else
            
            overlap_act = nansum(PCs.status(mask_act,sm,2));
            overlap_PC = nansum(any(PCs.status(mask_PC,sm,3:5),3));
            
            rand_pull_act = zeros(L,1);
            rand_pull_PC = zeros(L,1);
            
            K1_act = round(sum(nROI(s,2:5)));
            K2_act = round(sum(nROI(sm,2:5)));
            
            K1_PC = round(sum(nROI(s,3:5)));
            K2_PC = round(sum(nROI(sm,3:5)));
            
            offset = K1_act - overlap_act;
            
            for l = 1:L
              rand_pull_act(l) = nnz(ismember(datasample(1:nC,K1_act,'Replace',false),1:K2_act));
              rand_pull_PC(l) = nnz(ismember(datasample(1:K1_act,K1_PC,'Replace',false),datasample(offset+1:offset+K2_act,K2_PC,'Replace',false)));
            end
            
            pop_overlap.overrepr_act(s,sm) = overlap_act/mean(rand_pull_act);%(overlap_act-mean(rand_pull_act))/sqrt(var(rand_pull_act));
            pop_overlap.overrepr_PC(s,sm) = overlap_PC/mean(rand_pull_PC);%(overlap_PC-mean(rand_pull_PC))/sqrt(var(rand_pull_PC));
            
            pop_overlap.overrepr_act_var(s,sm) = std(rand_pull_act)/mean(rand_pull_act);
            pop_overlap.overrepr_PC_var(s,sm) = std(rand_pull_PC)/mean(rand_pull_PC);
            
            pop_overlap.frac_act(s,sm) = overlap_act/sum(mask_act);
            pop_overlap.frac_PC(s,sm) = overlap_PC/sum(mask_PC);
          end
        end
      end
    end
    fig_pos = [100 100 600 400]
    figure('position',fig_pos)
    subplot(2,2,1)
    imagesc(pop_overlap.overrepr_act)
    cb = colorbar;
    ylabel(cb,'Overrepresentation')
    xlabel('Session')
    ylabel('Session')
    title('active ROIs')
    caxis([1,2])
    
    subplot(2,2,2)
    imagesc(pop_overlap.frac_act)
    cb = colorbar;
    ylabel(cb,'Overlap')
    xlabel('Session')
    ylabel('Session')
    title('active ROIs')
    caxis([0,1])
    
    subplot(2,2,3)
    imagesc(pop_overlap.overrepr_PC)
    colormap('jet')
    cb = colorbar;
    ylabel(cb,'Overrepresentation')
    xlabel('Session')
    ylabel('Session')
    title('vPCs')
    caxis([1,2])
    
    subplot(2,2,4)
    imagesc(pop_overlap.frac_PC)
    cb = colorbar;
    ylabel(cb,'Overlap')
    xlabel('Session')
    ylabel('Session')
    title('vPCs')
    caxis([0,1])
    
    if sv
      save_fig(pathFigures,mouse,'Overlap',sv_suffix,sv_ext,fig_pos)
    end
    
%      N=round(mean(sum(nROI(:,2:5),2)))
%      L=500;
%      rand_pull_act=zeros(L,1);
%      for l=1:L
%        rand_pull_act(l) = nnz(ismember(datasample(1:nC,N,'Replace',false),1:N));
%      end
%      histogram(rand_pull_act)
    
    overlap = zeros(nSes,nSes-1)*NaN;
    frac_act = zeros(nSes,nSes-1)*NaN;
    
    fig_pos = [100 100 400 300];
    overlap_var = zeros(nSes,nSes-1)*NaN;
    figure('Position',fig_pos)
    ax1 = axes('Position',[0.2,0.175,0.75,0.375]);
    hold(ax1,'on')
    ax2 = axes('Position',[0.2,0.6,0.75,0.375]);
    hold(ax2,'on')
    
    for ds = 1:nSes-1
      for s = 1:nSes-ds
        if sum(nROI(s,2:5)) > 500 && sum(nROI(s+ds,2:5)) > 500
          overlap(s,ds) = pop_overlap.overrepr_act(s,s+ds);
          overlap_var(s,ds) = pop_overlap.overrepr_act_var(s,s+ds);
          
          frac_act(s,ds) = pop_overlap.frac_act(s,s+ds);
          
          plot(ax1,ds,pop_overlap.overrepr_act(s,s+ds),'.','Color',[0.7,0.7,0.7],'MarkerSize',5,'HandleVisibility','off')
          plot(ax2,ds,pop_overlap.frac_act(s,s+ds),'.','Color',[0.7,0.7,0.7],'MarkerSize',5,'HandleVisibility','off')
          
        end
      end
    end
    
%      for s = 1:nSes
%        plot(ax1,1:nSes-1,overlap(s,:),'-','Color',[0.9,0.9,0.9],'LineWidth',0.5,'HandleVisibility','off')
%        plot(ax2,1:nSes-1,frac_act(s,:),'-','Color',[0.9,0.9,0.9],'LineWidth',0.5,'HandleVisibility','off')
%      end
    
    [ol_median, ol_CI] = bootstrapping(overlap,1000,'mean');
    plot_CI_as_fill(ol_median,ol_CI,t_data,ax1,{[0,0,0],[0.5,0.5,0.5]},'data')
    plot_CI_as_fill(ones(nSes-1,1),nanmean(overlap_var),t_data,ax1,{[1,0,0],[1,0.5,0.5]},'random')
    
    xlabel(ax1,'\Delta s','FontSize',14)
    ylabel(ax1,'overrepr.','FontSize',14)
    
    xlim(ax1,[0,nSes])
    ylim(ax1,[0.8,2])
    legend(ax1,'Location',[0.65,0.4,0.3,0.1])
    
    [frac_median, frac_CI] = bootstrapping(frac_act,1000,'mean');
    plot_CI_as_fill(frac_median,frac_CI,t_data,ax2,{[0,0,0],[0.5,0.5,0.5]},'data')
    ylabel(ax2,'% overlap','FontSize',14)
    set(ax2,'XTick',[],'Xlim',[0,nSes],'Ylim',[0,1])
    
    suptitle('Active population')
    
    if sv
      save_fig(pathFigures,mouse,'Drift_active',sv_suffix,sv_ext,fig_pos)
    end
    
    overlap_PC = zeros(nSes,nSes-1)*NaN;
    frac_PC = zeros(nSes,nSes-1)*NaN;
    overlap_PC_var = zeros(nSes,nSes-1)*NaN;
    figure('Position',fig_pos)
    ax1 = axes('Position',[0.2,0.175,0.75,0.375]);
    hold(ax1,'on')
    ax2 = axes('Position',[0.2,0.6,0.75,0.375]);
    hold(ax2,'on')
    for ds = 1:nSes-1
      for s = 1:nSes-ds
        if sum(nROI(s,3:5)) > 50 && sum(nROI(s+ds,3:5)) > 50
          overlap_PC(s,ds) = pop_overlap.overrepr_PC(s,s+ds);
          overlap_PC_var(s,ds) = pop_overlap.overrepr_PC_var(s,s+ds);
          
          frac_PC(s,ds) = pop_overlap.frac_PC(s,s+ds);
          
          plot(ax1,ds,pop_overlap.overrepr_PC(s,s+ds),'.','Color',[0.7,0.7,0.7],'MarkerSize',5,'HandleVisibility','off')
          plot(ax2,ds,pop_overlap.frac_PC(s,s+ds),'.','Color',[0.7,0.7,0.7],'MarkerSize',5,'HandleVisibility','off')
        end
      end
    end
    
%      for s = 1:nSes
%        plot(ax1,1:nSes-1,overlap_PC(s,:),'-','Color',[0.9,0.9,0.9],'LineWidth',0.5,'HandleVisibility','off')
%        plot(ax2,1:nSes-1,frac_PC(s,:),'-','Color',[0.9,0.9,0.9],'LineWidth',0.5,'HandleVisibility','off')
%      end
    
    [ol_median, ol_CI] = bootstrapping(overlap_PC,1000,'mean');
    plot_CI_as_fill(ol_median,ol_CI,t_data,ax1,{[0,0,0],[0.5,0.5,0.5]},'data')
    plot_CI_as_fill(ones(nSes-1,1),nanmean(overlap_PC_var),t_data,ax1,{[1,0,0],[1,0.5,0.5]},'random')
    
    xlabel(ax1,'\Delta s','FontSize',14)
    ylabel(ax1,'overrepr.','FontSize',14)
    
    ylabel(ax2,'% overlap','FontSize',14)
    
    xlim(ax1,[0,nSes])
    ylim(ax1,[0.8,2])
    legend(ax1,'Location',[0.65,0.4,0.3,0.1])
    
    [frac_median, frac_CI] = bootstrapping(frac_PC,1000,'mean');
    plot_CI_as_fill(frac_median,frac_CI,t_data,ax2,{[0,0,0],[0.5,0.5,0.5]},'data')
    ylabel(ax2,'% overlap','FontSize',14)
    set(ax2,'XTick',[],'Xlim',[0,nSes],'Ylim',[0.0,1])
    
    suptitle('Place cell population')
    
    %% plot # of neurons coding in session x and y
    if sv
      save_fig(pathFigures,mouse,'Drift_PC',sv_suffix,sv_ext,fig_pos)
    end
  end
  
  
  
  if plot_fig(12)
    
    figure('position',[100 100 1500 1200])
    explode = [0,0,1,1,1];
    for s = 1:s_end
      subplot(3,5,s)
      pie(nROI(s,:)./nC,explode,{'silent','nvPC','NRNG','GT','RW'})
      title(gca,sprintf('Session %02d',s))
    end
    csvFile = pathcat(pathMouse,sprintf('ROInumbers.csv'));
    dlmwrite(csvFile,nROI,'delimiter',';')
    
  end
  
  
  if plot_fig(13)
  
    mean_nu = zeros(nSes,2);
    
    figure('position',[100 100 400 200])
    for s = 1:s_end
%        pathSession = pathcat(pathMouse,sprintf('Session%02d',s));
%        pathBH = dir(pathcat(pathSession,'deconvData.mat'));
%        pathCa = pathcat(pathSession,'results_OnACID.mat');
      
%        act = load(pathCa,'S');
%        S = act.S;
      
      
%        nu = zeros(size(S,1),3);
%        for n = 1:size(S,1)
%          nu(n,1) = sum(floor(sqrt(S(n,:)/(act.baseline10(n)*3))));
%          nu(n,2) = sum(floor(sqrt(S(n,:)/(act.baseline20(n)*3))));
%          nu(n,3) = sum(floor(sqrt(S(n,:)/(act.baseline30(n)*3))));
%        end
%        nu = nu./(8989/15.);
%        mean_nu(s,:) = nanmean(nu);
      activeCells = any(PCs.status(:,s,2:5),3);
      mean_nu(s,1) = mean(PCs.firingrate(:,s));
      mean_nu(s,2) = mean(PCs.firingrate(activeCells,s));
      
%        for c = find(activeCells)'
%          disp('-----')
%          n = PCs.neuronID(c,s,2);
%          [c,n]
%          [spikeNr,md,sd_r] = get_spikeNr(S(S(:,n)>0,n),linspace(1/15,600,8989));
%          spikeNr
%          [spikeNr/(8989/15), PCs.firingrate(c,s)]
%          
%          if isnan(PCs.firingrate(c,s)) || ~PCs.firingrate(c,s)
%            p1 = plot(S(:,n),'r');
%          end
%          if ~isnan(PCs.firingrate(c,s)) && PCs.firingrate(c,s)
%            p1 = plot(S(:,n),'b');
%          end
%          hold on
%          p2 = plot([0,9000],[md+2*sd_r md+2*sd_r],'k:','LineWidth',2);
%          waitforbuttonpress
%  %          close('all')
%          hold off
%          delete(p1)
%          delete(p2)
%        end
      
      
      if s == 5
%          subplot(1,2,1)
        hold on
        histogram(PCs.firingrate(activeCells,s),linspace(0,10,51),'Normalization','probability','HandleVisibility','off')
%          bar(0,sum(~activeCells),'FaceColor','r','HandleVisibility','off')
%          ylim([0,35])
        y_lims = ylim;
        
        plot([mean_nu(s,1) mean_nu(s,1)],[0,y_lims(2)],'r:','LineWidth',2,'DisplayName','$\bar{\nu}$ (all)')
        plot([mean_nu(s,2) mean_nu(s,2)],[0,y_lims(2)],'r--','LineWidth',2,'DisplayName','$\bar{\nu}$ (active only)')
%          xlim([0,10])
        xlabel('\nu')
        ylabel('p(\nu)')
        leg = legend();
        set(leg,'interpreter','latex')
      end
    end
    ax1 = axes('position',[0.65,0.45,0.3,0.25]);
    hold on
%      plot(mean_nu(:,1),'r:')
%      plot(mean_nu(:,2),'r-.')
    plot(mean_nu(:,1),'r:')
    plot(mean_nu(:,2),'r--')
    ylim([0,7])
    xlabel('session')
    ylabel('$\bar{\nu}$','interpreter','latex')
    
    if sv
      path = pathcat(pathFigures,sprintf('activityRate%s.%s',sv_suffix,sv_ext));
      print(path,sprintf('-d%s',sv_ext),'-r300')
      disp(sprintf('Figure saved as %s',path))
    end
     
  end
  
  
  
%    if plot_fig(14)
%      
%  %      f1 = figure('position',[100 100 1800 1500]);
%  %      f2 = figure('position',[100 100 1800 1500]);
%  %      f3 = figure('position',[100 100 1800 1500]);
%      
%      explode = [0,1,1,1];
%      
%      xls_in_dat = zeros(nSes,5);
%      xls_out_dat = zeros(nSes,5);
%      
%      for s = 1:s_end
%        
%        title_str = sprintf('Session %02d',s);
%        
%        remap_hist_tot = squeeze(remap_pos(:,1,:,:));
%        
%  %        occupancy_phase = occupancy(s,:);
%        
%        
%        silent_any = zeros(para.nbin,1);
%        any_silent = zeros(1,para.nbin);
%        
%        silent_out_trans = zeros(5,1);
%        
%        if s < 15
%          
%  %          silent_any = squeeze(remap_hist_tot(s,para.idx_zone{1},1:para.nbin)/sum(occupancy(s,para.idx_zone{1})));
%          silent_any = squeeze(remap_hist_tot(s,para.idx_zone{1},1:para.nbin));
%          
%          silent_out_trans(1) = remap_hist_tot(s,para.idx_zone{1},para.idx_zone{1});%/sum(occupancy(s,para.idx_zone{1}));
%          silent_out_trans(2) = remap_hist_tot(s,para.idx_zone{1},para.idx_zone{2});%/sum(occupancy(s,para.idx_zone{1}));
%          silent_out_trans(3) = sum(silent_any(para.idx_zone{3}));
%          silent_out_trans(4) = sum(silent_any(para.idx_zone{4}));
%          silent_out_trans(5) = sum(silent_any(para.idx_zone{5}));
%          
%          xls_out_dat(s,:) = silent_out_trans;
%          
%  %          figure(f2)
%  %          ax_pie = subplot(3,5,s);%axes('Position',[ax_pos(1)+0.15 ax_pos(2)+0.1 0.075 0.05]);
%  %    %        hold on
%  %    %        pie(silent_trans,explode,{'silent','nvPC','NRNG','GT','RW'})
%  %          perc = 100*silent_out_trans(2:5)/sum(silent_out_trans(2:5));
%  %          pie(silent_out_trans(2:5)/sum(silent_out_trans(2:5)),explode,{sprintf('nPC %3.1f%%',perc(1)),sprintf('NGNR %3.1f%%',perc(2)),sprintf('GT %3.1f%%',perc(3)),sprintf('RW %3.1f%%',perc(4))})
%  %  %          title(sprintf('Session %02d',s))
%        end
%        
%        silent_in_trans = zeros(5,1);
%        
%        if s > 1
%  %          any_silent = squeeze(remap_hist_tot(s-1,1:para.nbin,para.idx_zone{1})/sum(occupancy(s,para.idx_zone{1})));
%          any_silent = squeeze(remap_hist_tot(s-1,1:para.nbin,para.idx_zone{1}));
%          
%          silent_in_trans(1) = remap_hist_tot(s-1,para.idx_zone{1},para.idx_zone{1});%/sum(occupancy(s,para.idx_zone{1}));
%          silent_in_trans(2) = remap_hist_tot(s-1,para.idx_zone{2},para.idx_zone{1});%/sum(occupancy(s,para.idx_zone{1}));
%          silent_in_trans(3) = sum(any_silent(para.idx_zone{3}));
%          silent_in_trans(4) = sum(any_silent(para.idx_zone{4}));
%          silent_in_trans(5) = sum(any_silent(para.idx_zone{5}));
%          
%          xls_in_dat(s,:) = silent_in_trans;
%  %          figure(f3)
%  %          perc = 100*silent_in_trans(2:5)/sum(silent_in_trans(2:5));
%  %          ax_pie = subplot(3,5,s);%axes('Position',[ax_pos(1)+0.15 ax_pos(2)+0.1 0.075 0.05]);
%  %          pie(silent_in_trans(2:5)/sum(silent_in_trans(2:5)),explode,{sprintf('nPC %3.1f%%',perc(1)),sprintf('NGNR %3.1f%%',perc(2)),sprintf('GT %3.1f%%',perc(3)),sprintf('RW %3.1f%%',perc(4))})
%  %  %          title(sprintf('Session %02d',s))
%        end
%        
%        x_margin = 0.03;
%        y_margin = 0.0;
%        
%  %        ylims = [40,20,30,30,30];
%        ylims = [0.05,0.02,0.02,0.01,0.01];
%        
%        
%        
%  %        if j == 1
%  %          ylabel('PC_{gate}')
%  %          text(50,3/4*ylims(1),'recruiting from')
%  %          text(50,-3/4*ylims(1),'get recruited to')
%  %        end
%        
%  %        figure(f1)
%  %        ax = subplot(3,5,s);
%  %        ax_pos = get(ax,'Position');
%  %        hold on
%  %        bar(any_silent,'g','EdgeColor','none')
%  %        bar(-silent_any,'r','EdgeColor','none')
%  %        bar(any_silent-silent_any','b','EdgeColor','none')
%  %        hold off
%  %        set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5))
%  %  %        if j == 1
%  %  %          ylabel('silent')
%  %  %        end
%  %        ylim([-ylims(5),ylims(5)])
%  %        xlabel('Location at s [cm]')
%  %        title(sprintf('Session %02d',s))
%      end
%      
%      csvFile = pathcat(pathMouse,sprintf('ROI_transitions_out.csv'));
%      dlmwrite(csvFile,xls_out_dat,'delimiter',';')
%      
%      csvFile = pathcat(pathMouse,sprintf('ROI_transitions_in.csv'));
%      dlmwrite(csvFile,xls_in_dat,'delimiter',';')
%      
%  %      if sv
%  %        figure(f1)
%  %        path = pathcat(pathFigures,sprintf('silence_bars.png'));
%  %        print(path,'-dpng','-r300')
%  %        
%  %        figure(f2)
%  %        suptitle('recruiting to...')
%  %        path = pathcat(pathFigures,sprintf('silence_pie_out.png'));
%  %        print(path,'-dpng','-r300')
%  %        
%  %        figure(f3)
%  %        suptitle('recruiting from...')
%  %        path = pathcat(pathFigures,sprintf('silence_pie_in.png'));
%  %        print(path,'-dpng','-r300')
%  %      end
%    end
  
  if plot_fig(16)
    
    figure('position',[100 100 600 400])
    hold on
    ylim([-para.nbin/2,para.nbin/2])
    plot([0,16],[0,0],'k--')
    set(gca,'YTick',linspace(-para.nbin/2,para.nbin/2,5),'YTickLabel',linspace(-50,50,5))
    xlabel('Session')
    ylabel('dx [cm]')
    
    idx = find(sum(PCs.status >= 3,2)>=10);
    
    ref_pos = 0;
    
    j = 0;
    trace = zeros(1,nSes)*NaN;
    for i_0 = idx'
      i = PCs.clusterID(i_0,2);
      
%        field_ct = max([PC_clusters(i).session.field_ct]);
      ref_pos = 0;
      breakit = false;
      
      for s = nSes:-1:1
        if any(PCs.status(i,s,3:5))
          
          [max_val, max_idx] = max([frac_active(i,s,:)]);
          if ~ref_pos
            if max_val>=perc_thr
              ref_pos = PCs.fields.center(i,s,max_idx);
              if ~(ref_pos >= gt_zone(1) && ref_pos <= gt_zone(2))
%                if ~(ref_pos >= rw_zone(1) && ref_pos <= rw_zone(2))
%                if (ref_pos >= gt_zone(1) && ref_pos <= gt_zone(2)) || (ref_pos >= rw_zone(1) && ref_pos <= rw_zone(2))
                breakit = true;
                break
              else
                j = j+1;
                trace(j,1) = 0;
                trace(j,:) = trace(j,:)*NaN;
              end
            end
          end
          df = mod(PCs.fields.center(i,s,max_idx)-ref_pos + para.nbin/2,para.nbin)-para.nbin/2;
          trace(j,s) = df;
        end
      end
      if ~breakit
        if j > 1
          set(plt_trace,'Color',[0.8,0.8,0.8],'LineWidth',0.2)
        end
        plt_trace = plot(trace(j,:),'-','Color','k');
%          title(sprintf('Neuron #%d',i))
%          path = pathcat(pathFigures,sprintf('traces/PC_trajectories_GT10_%d%s.png',j,sv_suffix));
%          print(path,'-dpng','-r300')
%          waitforbuttonpress
      end
    end
    
    errorbar(linspace(1,nSes,nSes),nanmean(trace,1),sqrt(nanvar(trace,1)),'k-','LineWidth',3)
    sv
    if sv
      path = pathcat(pathFigures,sprintf('PC_trajectories_GT10%s.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
    
  end
  
  
  
  
  if plot_fig(18)
    
    pPop = zeros(5,5);
%      nOcc
    for p = 1:5
%        sum(nOcc(idx_pop(:,i),:),1);
%        pPop(i,:) = sum(nOcc(idx_pop(:,i),:),1) / (sum(idx_pop(:,i))*nSes); %% overall probability
      nOcc_tmp = nOcc;
%        nOcc_tmp(idx_pop(:,p),:)
      nOcc_tmp(:,p) = nOcc_tmp(:,p) - pop_thr;
      pPop(p,:) = sum(nOcc_tmp(idx_pop(:,p),:),1) / (sum(idx_pop(:,p))*(nSes-pop_thr)); %% overall probability
%          pPop_corr(:,i) = sum(nOcc(idx,:),1) / sum((sum(nOcc(idx,:),2)-nOcc(idx,i)));
    end
    
    figure('position',[500 500 1200 400])
    
    axes('position',[0.05 0.12 0.2 0.85])
    explode = [0,0,1,1,1];
    h_pie = pie(pTotal,explode,{'silent','nvPC','NRNG','GT','RW'});
    
    set(h_pie(1),'FaceColor',[1 1 1])
    set(h_pie(3),'FaceColor',[0 0 0])
    set(h_pie(5),'FaceColor',[0 0 1])
    set(h_pie(7),'FaceColor',[0 1 0])
    set(h_pie(9),'FaceColor',[1 0 0])
%      set(h_pie(11),'FaceColor',[0 1 0])
    
          
    axes('position',[0.35 0.175 0.25 0.8])
    imagesc(log10(pPop ./pTotal))
%        imagesc(log10(pPop_corr ./pTotal'))
    colormap(cm')
    h_cb = colorbar;
    set(h_cb,'YTick',linspace(-1,1,3),'YTickLabel',{'0.1','1','10'})
    ylabel(h_cb,'times average occurence','FontSize',14)
    caxis([-1,1])
    
    set(gca,'XTick',linspace(1,5,5),'XTickLabels',{'silent','nPC','NRNG','GT','RW'},...
            'YTick',linspace(1,5,5),'YTickLabels',{'silent','nPC','NRNG','GT','RW'})
    xlabel('coding for','FontSize',14)
    ylabel(sprintf('population (# coding >= %d)',pop_thr),'FontSize',14)
    xtickangle(60)
    
    
    axes('position',[0.7 0.175 0.25 0.8])
    hold on
    barh(nPop,'FaceColor',[0.4 0.4 0.4],'HandleVisibility','off')
    plot([nC,nC],[0.5,5.5],'r--','HandleVisibility','off')
    
    bar_stack = [sum(idx_pop(:,1)&idx_pop(:,3)) sum(idx_pop(:,1)&idx_pop(:,4)) sum(idx_pop(:,1)&idx_pop(:,5));...
                 sum(idx_pop(:,2)&idx_pop(:,3)) sum(idx_pop(:,2)&idx_pop(:,4)) sum(idx_pop(:,2)&idx_pop(:,5));...
                 sum(idx_pop(:,3)&idx_pure) sum(idx_pop(:,3)&idx_pop(:,4)) sum(idx_pop(:,3)&idx_pop(:,5));...
                 sum(idx_pop(:,4)&idx_pop(:,3)) sum(idx_pop(:,4)&idx_pure) sum(idx_pop(:,4)&idx_pop(:,5));...
                 sum(idx_pop(:,5)&idx_pop(:,3)) sum(idx_pop(:,5)&idx_pop(:,4)) sum(idx_pop(:,5)&idx_pure)];
                  
    h_bar = barh([1;2;3;4;5],bar_stack,'stacked');
    
    set(h_bar,{'FaceColor'},{'b';'r';'g'})
    set(h_bar,{'DisplayName'},{'NRNG';'GT';'RW'})
    ylim([0.5,5.5])
%      xlim([0,500])
    set(gca,'YTick',linspace(1,5,5),'YTickLabels',{'silent','nPC','NRNG','GT','RW'},'YDir','reverse')
    
    leg = legend('location','SouthEast');
    title(leg,'other populations')
%      XScale('log')
    xlim([1,1000])
    
    if sv
      path = pathcat(pathFigures,sprintf('populations%d%s.png',pop_thr,sv_suffix));
      print(path,'-dpng','-r300')
    end
    
  end
  
  
  
  if plot_fig(19)
    
    act_number = nSes - nOcc(:,1);
    PC_number = nSes - sum(nOcc(:,1:2),2);
    
    
    figure('position',[100 100 1000 700])
%      subplot(2,2,1)
    axes('position',[0.1 0.5 0.4 0.4])
%      title('PC vs active')
    hold on
    plot([0,nSes],[0,nSes],'k--')
%      jitter = 0.5*(0.5-rand(nC,2));
    sz = 20;
    scatter(act_number(idx_none)+0.5*(0.5-rand(sum(idx_none),1)),PC_number(idx_none)+0.5*(0.5-rand(sum(idx_none),1)),sz/2,'ko','filled')
    scatter(act_number(GT_only)+0.5*(0.5-rand(sum(GT_only),1)),PC_number(GT_only)+0.5*(0.5-rand(sum(GT_only),1)),sz,'go','filled')
    scatter(act_number(RW_only)+0.5*(0.5-rand(sum(RW_only),1)),PC_number(RW_only)+0.5*(0.5-rand(sum(RW_only),1)),sz,'ro','filled')
    scatter(act_number(NRNG_only)+0.5*(0.5-rand(sum(NRNG_only),1)),PC_number(NRNG_only)+0.5*(0.5-rand(sum(NRNG_only),1)),sz,'bo','filled')
    scatter(act_number(idx_mixed)+0.5*(0.5-rand(sum(idx_mixed),1)),PC_number(idx_mixed)+0.5*(0.5-rand(sum(idx_mixed),1)),sz,'mo','filled')
%      scatter(act_number+jitter(:,1),PC_number+jitter(:,2),2,'x')
    xlim([-0.5,nSes+0.5])
    ylim([-0.5,nSes+0.5])
    
    set(gca,'XAxisLocation','top','FontSize',12)
    ylabel('# Session vPC','FontSize',14)
    xlabel('# Session active','FontSize',14)
    
    axes('position',[0.16 0.74 0.13 0.1])
    PC_ratio = PC_number./act_number;
    hold on
    scatter(2*ones(sum(idx_none),1),PC_ratio(idx_none),act_number(idx_none),'ko','filled')
    scatter(3*ones(sum(NRNG_only),1),PC_ratio(NRNG_only),act_number(NRNG_only),'bo','filled')
    scatter(4*ones(sum(GT_only),1),PC_ratio(GT_only),act_number(GT_only),'go','filled')
    scatter(5*ones(sum(RW_only),1),PC_ratio(RW_only),act_number(RW_only),'ro','filled')
    scatter(6*ones(sum(idx_mixed),1),PC_ratio(idx_mixed),act_number(idx_mixed),'mo','filled')
    
    errorbar(2,mean(PC_ratio(idx_none)),std(PC_ratio(idx_none)),'kx','LineWidth',2)
    errorbar(3,mean(PC_ratio(NRNG_only)),std(PC_ratio(NRNG_only)),'kx','LineWidth',2)
    errorbar(4,mean(PC_ratio(GT_only)),std(PC_ratio(GT_only)),'kx','LineWidth',2)
    errorbar(5,mean(PC_ratio(RW_only)),std(PC_ratio(RW_only)),'kx','LineWidth',2)
    errorbar(6,mean(PC_ratio(idx_mixed)),std(PC_ratio(idx_mixed)),'kx','LineWidth',2)
    xlim([1.5,6.5])
    ylim([0,1])
    set(gca,'XTick',linspace(2,6,5),'XTickLabel',{'nPC','NRNG','GT','RW','mixed'},'FontSize',10)
    xtickangle(60)
    ylabel('#PC / #active')
    
%      suptitle('add fraction of GT/RW/NRNG/multi-population per session active / session PC')
%      subplot(2,2,3)
    axes('position',[0.1 0.3 0.4 0.2])
    hold on
    
    h_edges = linspace(-0.5,nSes+0.5,17);
    h_centers = linspace(0,nSes,16);
    act_hist = histcounts(act_number,h_edges)';
    histogram(act_number,h_edges,'FaceColor','k')
    
    GT_hist = histcounts(nSes-nOcc(GT_only,1),h_edges)';%sum(nOcc(GT_only,2:5),2),h_edges)';
    RW_hist = histcounts(nSes-nOcc(RW_only,1),h_edges)';%sum(nOcc(RW_only,2:5),2),h_edges)';
    NRNG_hist = histcounts(nSes-nOcc(NRNG_only,1),h_edges)';%sum(nOcc(NRNG_only,2:5),2),h_edges)';
    multi_hist = histcounts(nSes-nOcc(idx_mixed,1),h_edges)';%sum(nOcc(idx_mixed,2:5),2),h_edges)';
    
    h_bar = bar(h_centers,[NRNG_hist,GT_hist,RW_hist,multi_hist],1,'stack','FaceAlpha',0.8);
    set(h_bar,{'FaceColor'},{'b';'g';'r';'m'})
    
    
    set(gca,'YDir','reverse','YAxisLocation','right','XTick',[],'FontSize',12)
    xlim([-0.5,nSes+0.5])
    ylabel('# cells','FontSize',14)
    
    axes('position',[0.1 0.175 0.4 0.1])
    hold on
%      plot(nPC_hist/nPop_GT,'r')
%      plot(h_centers,GT_hist/nPop_GT,'r')
%      plot(h_centers,RW_hist/nPop_RW,'g')
%      plot(h_centers,NRNG_hist/nPop_NRNG,'b')
%      plot(h_centers,multi_hist/nPop_multi,'m')
    
    plot(h_centers,GT_hist./act_hist,'g')
    plot(h_centers,RW_hist./act_hist,'r')
    plot(h_centers,NRNG_hist./act_hist,'b')
    plot(h_centers,multi_hist./act_hist,'m')
    xlim([-0.5,nSes+0.5])
    ylim([0,0.4])
%      xticks([])
    
    
%      subplot(2,2,2)
    axes('position',[0.5 0.5 0.2 0.4])
    hold on
    GT_hist = histcounts(nSes-sum(nOcc(GT_only,1:2),2),h_edges)';%sum(nOcc(GT_only,3:5),2),h_edges)';
    RW_hist = histcounts(nSes-sum(nOcc(RW_only,1:2),2),h_edges)';%sum(nOcc(RW_only,3:5),2),h_edges)';
    NRNG_hist = histcounts(nSes-sum(nOcc(NRNG_only,1:2),2),h_edges)';%sum(nOcc(NRNG_only,3:5),2),h_edges)';
    multi_hist = histcounts(nSes-sum(nOcc(idx_mixed,1:2),2),h_edges)';%sum(nOcc(idx_mixed,3:5),2),h_edges)';
    
    PC_num_hist = histcounts(PC_number,h_edges)';
    histogram(PC_number,h_edges,'FaceColor','k','orientation','horizontal')
    h_bar = barh(h_centers,[NRNG_hist,GT_hist,RW_hist,multi_hist],1,'stack','FaceAlpha',0.8);
    set(h_bar,{'FaceColor'},{'b';'g';'r';'m'})
    
    xlim([0,400])
    ylim([-0.5,nSes+0.5])
    set(gca,'YTick',[],'FontSize',12)
    xlabel('# cells','FontSize',14)
    
    axes('position',[0.725 0.5 0.1 0.4])
    hold on
    
    plot(GT_hist./PC_num_hist,h_centers,'g')
    plot(RW_hist./PC_num_hist,h_centers,'r')
    plot(NRNG_hist./PC_num_hist,h_centers,'b')
    plot(multi_hist./PC_num_hist,h_centers,'m')
%      yticks([])
    ylim([-0.5,nSes+0.5])
    xlim([0,1])
    set(gca,'YAxisLocation','right','FontSize',12)
    
    
%      jitter = zeros(nC,2);
%      axes('position',[0.8 0.1 0.15 0.15])
%      hold on
%      plot([0,15],[0,15],'k--')
%      jitter = 0.5*(0.5-rand(nC,2));
%      scatter(nOcc(:,4)+jitter(:,1),nOcc(:,5)+jitter(:,2),2,'x')
%      xlim([-0.5,15.5])
%      ylim([-0.5,15.5])
%      xlabel('GT')
%      ylabel('RW')
%      
%      axes('position',[0.8 0.3 0.15 0.15])
%      hold on
%      plot([0,15],[0,15],'k--')
%      jitter = 0.5*(0.5-rand(nC,2));
%      scatter(nOcc(:,4)+jitter(:,1),nOcc(:,3)+jitter(:,2),2,'x')
%      xlim([-0.5,15.5])
%      ylim([-0.5,15.5])
%      xlabel('GT')
%      ylabel('NRNG')
%      
%      axes('position',[0.6 0.1 0.15 0.15])
%      hold on
%      plot([0,15],[0,15],'k--')
%      jitter = 0.5*(0.5-rand(nC,2));
%      scatter(nOcc(:,3)+jitter(:,1),nOcc(:,5)+jitter(:,2),2,'x')
%      xlim([-0.5,15.5])
%      ylim([-0.5,15.5])
%      xlabel('NRNG')
%      ylabel('RW')
    
    if sv
      path = pathcat(pathFigures,sprintf('activity%d%s.png',pop_thr,sv_suffix));
      print(path,'-dpng','-r300')
    end
  end
  
%    if plot_fig(20)
%      for c = cell_idx(GT_only)'
%        plot_clusters('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/M879',PC_clusters,c,0.3)
%      end
%    end
  
  if plot_fig(21)
  
    %%% firing rates and MI of "special cells":
    
    %% test, whether "coding" cells have higher / lower firing rate
    %% -> doesn't seem like it
    %% also, do simulation, whether random spiking cells (later, cells with field) are detected to be PCs, depending on rate!
    
    %% test, whether GT/RW/NRNG cells have specific firing rate
    %% -> doesn't seem like it
    
    %% test, whether firing rate is modulated by coding / non-coding
    
    %% test, whether firing rate is modulated in correlated fashion within one group!
    
    h_bins = linspace(0,prctile(PCs.firingrate(:),99),41);
    
    figure('position',[100 100 600 400])
    axes('position',[0.15 0.4 0.8 0.5])
    hold on
    
    offset = (h_bins(2)-h_bins(1))/2;
    act_hist = histc(PCs.firingrate(PCs.status(:,:,2)),h_bins-offset);
    NRNG_hist = histc(PCs.firingrate(PCs.status(:,:,3)),h_bins-offset);
    GT_hist = histc(PCs.firingrate(PCs.status(:,:,4)),h_bins-offset);
    RW_hist = histc(PCs.firingrate(PCs.status(:,:,5)),h_bins-offset);
    
%      GT_hist = histc(PCs.firingrate(GT_only),h_bins);
%      RW_hist = histc(PCs.firingrate(RW_only),h_bins);
%      NRNG_hist = histc(PCs.firingrate(NRNG_only),h_bins);
    
    histogram(PCs.firingrate(PCs.status(:,:,2)),h_bins-offset,'FaceColor',[0.3 0.3 0.3],'DisplayName','active')
    
    h_bar = bar(h_bins,[GT_hist,RW_hist,NRNG_hist],1,'stack');
    set(h_bar,{'FaceColor'},{'r';'g';'b'})
%      set(gca,'XTick',[])
    set(h_bar,{'DisplayName'},{'GT';'RW';'NRNG'})
%      xlim([0,400])
    ylabel('number of ROIs')
    legend()
    
    
    axes('position',[0.15 0.1 0.8 0.25])
    hold on
    plot(h_bins,GT_hist./act_hist,'r')
    plot(h_bins,RW_hist./act_hist,'g')
    plot(h_bins,NRNG_hist./act_hist,'b')
%      xlim([0,400])
    xlabel('activity rate')
    ylabel('fraction')
    
    if sv
      path = pathcat(pathFigures,sprintf('populations_firingrate%d%s.png',pop_thr,sv_suffix));
      print(path,'-dpng','-r300')
    end
  end
  
  
  if plot_fig(22)
    
%      population = 4;   %% 2=active, 3=NRNG, 4=GT, 5=RW
    for population = [3,4,5]
      figure('position',[100 100 400 600])
      
      c_idx = find(sum(PCs.status(:,:,population),2)>=pop_thr);
      
      color_arr = zeros(length(c_idx),nSes);
      for i = 1:length(c_idx)
        c = c_idx(i);
        
        for s = 1:s_end
          color_arr(i,s) = sum(PCs.status(c,s,2:5));
          if color_arr(i,s) > 2  %% coding for more than one thing
            color_arr(i,s) = 5;
          else
            if PCs.status(c,s,4)
              color_arr(i,s) = 3;
            elseif PCs.status(c,s,5)
              color_arr(i,s) = 4;
            end
          end
        end
      end
      
      cm = ones(3,6);
      cm(:,1) = [1 1 1];  %% silent
      cm(:,2) = [0 0 0];  %% active nPC
      cm(:,3) = [0 0 1];  %% NRNG PC
      cm(:,4) = [1 0 0];  %% GT PC
      cm(:,5) = [0 1 0];  %% RW PC
      cm(:,6) = [1 0 1];  %% multi PC
      
      imagesc(color_arr)
      ylim([0.5,min(50,length(c_idx))+0.5])
      ylabel('ID')
      xlabel('Session')
      switch population
        case 3
          title('NRNG vPC population')
          sv_str = 'NRNG';
        case 4
          title('GT vPC population')
          sv_str = 'GT';
        case 5
          title('RW vPC population')
          sv_str = 'RW';
      end
      
      colormap(cm')
      h_cbar = colorbar('Location','SouthOutside');
      caxis([-0.5,5.5])
      
      set(h_cbar,'YTick',linspace(0,5,6),'YTickLabel',{'silent','nPC','NRNG','GT','RW','multi PC'})
      
      if sv
        path = pathcat(pathFigures,sprintf('raster%d_%s%s.png',pop_thr,sv_str,sv_suffix));
        print(path,'-dpng','-r300')
      end
    end
  end
  
  if plot_fig(23)
    cm = ones(3,21);
    cm(1,1:11) = linspace(0,1,11);    %% red
    cm(2,1:11) = linspace(0,1,11);  %% green
    
    cm(3,11:21) = linspace(1,0,11);   %% blue
    cm(2,11:21) = linspace(1,0,11);   %% green
    
    %% firing rate correlation overall, in populations, in "pure" populations
    population = 5;
    c_idx = find(sum(PCs.status(:,:,population),2)>=pop_thr);
    
    fr_real = PCs.firingrate(c_idx,:)';
    fr_corr = corr(fr_real);
%      fr_corr(logical(eye(length(c_idx)))) = 0;
%      mean(nanmean(fr_corr))
    
    figure
    subplot(1,2,1)
    imagesc(fr_corr)
    colormap(cm')
    colorbar
    caxis([-1,1])
    
    subplot(1,2,2)
%      c_idx = find(sum(PCs.status(:,:,4),2)>=4);
%      fr_real = PCs.firingrate(c_idx,:)';
%      fr_corr= corr(fr_real);
    
    n  = size(fr_corr,1);             % number of nodes
    M  = 1:n;                   % initial community affiliations
    Q0 = -1; Q1 = 0;            % initialize modularity values
    while Q1-Q0>1e-5;           % while modularity increases
      Q0 = Q1;                % perform community detection
      [M, Q1] = community_louvain(fr_corr, 5, M,'negative_sym');
    end
    [~,ind] = sort(M);
    imagesc(fr_corr(ind,ind))
    colormap(cm')
    colorbar
    caxis([-1,1])
    
    figure
    histogram(fr_corr)
    
  end
  
  if plot_fig(24)
    
    xlin = linspace(0,1,51);
    
    figure('position',[500 500 400 300])
    
    
    axes('position',[0.05 0.7 0.7 0.25])
    mask_NRNG = PCs.fields.status == 3;
    hold on
    histogram(frac_active(mask_NRNG),xlin,'FaceColor','b','Normalization','probability')
    errorbar(mean(frac_active(mask_NRNG)),0,std(frac_active(mask_NRNG)),'k','horizontal','LineWidth',2)
    set(gca,'YTick',[],'XTick',[])
    xlim([0.4,1])
    ylim([0,0.15])
    mean_frac(1) = mean(frac_active(mask_NRNG));
    var_frac(1) = std(frac_active(mask_NRNG));
    legend({'NRNG'})
    
    axes('position',[0.05 0.425 0.7 0.25])
    mask_GT = PCs.fields.status == 4;
    hold on
    histogram(frac_active(mask_GT),xlin,'FaceColor','r','Normalization','probability')
    errorbar(mean(frac_active(mask_GT)),0,std(frac_active(mask_GT)),'k','horizontal','LineWidth',2)
    set(gca,'YTick',[],'XTick',[])
    xlim([0.4,1])
    ylim([0,0.15])
    mean_frac(2) = mean(frac_active(mask_GT));
    var_frac(2) = std(frac_active(mask_GT));
    legend({'GT'})
    
    
    axes('position',[0.05 0.15 0.7 0.25])
    mask_RW = PCs.fields.status == 5;
    hold on
    histogram(frac_active(mask_RW),xlin,'FaceColor','g','Normalization','probability')
    errorbar(mean(frac_active(mask_RW)),0,std(frac_active(mask_RW)),'k','horizontal','LineWidth',2)
    set(gca,'YTick',[])
    xlim([0.4,1])
    ylim([0,0.15])
    xlabel('reliability')
    mean_frac(3) = mean(frac_active(mask_RW));
    var_frac(3) = std(frac_active(mask_RW));
    legend({'RW'})
    
    
    axes('position',[0.8 0.15 0.15 0.8])
    hold on
    errorbar(linspace(1,3,3),mean_frac,var_frac,'ko')
    
    plot(1,mean_frac(1),'b.','MarkerSize',30)
    plot(2,mean_frac(2),'r.','MarkerSize',30)
    plot(3,mean_frac(3),'g.','MarkerSize',30)
    set(gca,'XTick',linspace(1,3,3),'XTickLabel',{'NRNG','GT','RW'},'YTick',linspace(0,1,3))
    ylim([0.4,1])
    xlim([0.5,3.5])
    xtickangle(60)
    
    if sv
      path = pathcat(pathFigures,sprintf('population_reliability%d_%s.png',pop_thr,sv_suffix));
      print(path,'-dpng','-r300')
    end
    
%      figure
%      errorbar(linspace(1,3,3),mean_frac,var_frac,'kx')
%      
%      xlin = linspace(0,30,31);
%      figure;
%      subplot(3,1,1)
%      mask_GT = PCs.fields.status == 4 & frac_active > 0.5;
%      hold on
%      histogram(PCs.fields.width(mask_GT),xlin,'FaceColor','r')
%      errorbar(mean(PCs.fields.width(mask_GT)),110,sqrt(var(PCs.fields.width(mask_GT))),'kx','horizontal','LineWidth',2)
%      xlim([0,30])
%      mean(PCs.fields.width(mask_GT))
%      
%      subplot(3,1,2)
%      mask_RW = PCs.fields.status == 5 & frac_active > 0.5;
%      hold on
%      histogram(PCs.fields.width(mask_RW),xlin,'FaceColor','g')
%      errorbar(mean(PCs.fields.width(mask_RW)),110,sqrt(var(PCs.fields.width(mask_RW))),'kx','horizontal','LineWidth',2)
%      xlim([0,30])
%      mean(PCs.fields.width(mask_RW))
%      
%      subplot(3,1,3)
%      mask_NRNG = PCs.fields.status == 3 & frac_active > 0.5;
%      hold on
%      histogram(PCs.fields.width(mask_NRNG),xlin,'FaceColor','b')
%      errorbar(mean(PCs.fields.width(mask_NRNG)),110,sqrt(var(PCs.fields.width(mask_NRNG))),'kx','horizontal','LineWidth',2)
%      xlim([0,30])
%      mean(PCs.fields.width(mask_NRNG))

  end
  
  if plot_fig(25)
    %% 2nd peaks
    nFields = sum(~isnan(PCs.fields.center),3);
    NRNG_idx = any(PCs.fields.status==3,3);
    GT_idx = any(PCs.fields.status==4,3);
    RW_idx = any(PCs.fields.status==5,3);
    
    mean(nFields(NRNG_idx))
    mean(nFields(GT_idx))
    mean(nFields(RW_idx))
    
%      PCs.fields.bins
    
%      GT_only
    
%      idx_pop(4)
    GTs = PCs.fields.center(idx_pop(:,4),:,:);
    GTs = GTs(~isnan(GTs));
%      GTs
    
    all(NRNG_only == idx_pop(:,3))
    sum(NRNG_only == idx_pop(:,3))
    sum(idx_pop(:,3))
    
    figure
    subplot(3,1,1)
    hold on
    histogram(PCs.fields.center(idx_pop(:,3),:,:),linspace(0,80,81))
    histogram(PCs.fields.center(NRNG_only,:,:),linspace(0,80,81),'FaceColor','r')
    subplot(3,1,2)
    hold on
    histogram(PCs.fields.center(idx_pop(:,4),:,:),linspace(0,80,81))
    histogram(PCs.fields.center(GT_only,:,:),linspace(0,80,81),'FaceColor','r')
    subplot(3,1,3)
    hold on
    histogram(PCs.fields.center(idx_pop(:,5),:,:),linspace(0,80,81))
    histogram(PCs.fields.center(RW_only,:,:),linspace(0,80,81),'FaceColor','r')
    
    %% histogram of 2nd fields, when coding for gt etc. (in population, or outside?!)
    %% histogram of fields, when not coding for population bins
    
  end
  
  if plot_fig(26)
    %%% plot stability of different populations (with tolerance = hard threshold)
    
    figure('position',[100 100 400 200])
%      ax1 = subplot(2,1,1);
    ax1 = subplot(1,1,1);
    hold on
%      ax2 = subplot(2,1,2);
%      hold on
    lbl_arr = {'other cell (NRNG)','gate cell (GT)','reward cell (RW)'};
    stable_dt = zeros(nT,N_bs,3);
    for p = 1:3
      
      N_data = length(arrays.ICI_s{p});
      if N_data > 0
        bs_idx = randi(N_data,N_bs,N_data);
          
        tic
        for L = 1:N_bs
          ICI_L = arrays.ICI_s{p}(bs_idx(L,:));
          shift_L = arrays.shift{p}(bs_idx(L,:));
          for ds = 1:nT
            idxes_ds = ICI_L==t_data(ds);
            stable_dt(ds,L,p) = sum(abs(shift_L(idxes_ds)) < tolerance)/sum(N_norm(ds,:,p),2);
          end
        end
        
        stable_mean = mean(stable_dt(:,:,p),2);
        stable_CI = prctile(stable_dt(:,:,p),[5,95],2);
        
        mask = ~isnan(stable_mean);
        x2 = [t_data(mask) fliplr(t_data(mask))];
        inBetween = [stable_CI(mask,1)', flipud(stable_CI(mask,2))'];
        fill(ax1,x2,inBetween,col_fill{p},'FaceAlpha',0.5,'EdgeColor','None','HandleVisibility','off');
        plot(ax1,t_data,stable_mean,'o','Color',col{p},'DisplayName',lbl_arr{p})
      end
%        plot(ax2,t_data,N_norm(:,p),'o-','Color',col{p})
    end
    set(ax1,'FontSize',12)
    xlim(ax1,[0,t_ses(end)])
    set(gca,'yscale','log')
    ylim([10^-2,2*10^-1])
%      ylim(ax1,[0,0.55])
    xlabel(ax1,'\Delta s','FontSize',14)
    ylabel(ax1,'r^*_{stable}','FontSize',14)
    legend(ax1,'Location',[0.52,0.68,0.4,0.3])
    
%      xlim(ax2,[0,t_ses(end)])
%      ylim(ax2,[0,3000])
    
    if sv
      path = pathcat(pathFigures,sprintf('stability_thr%s.%s',sv_suffix,sv_ext));
      print(path,sprintf('-d%s',sv_ext),'-r300')
      disp(sprintf('Figure saved as %s',path))
    end
  end
  
  
  if plot_fig(27)
  %% population sizes
    
    pop_thr_arr = 1:7;
    for p = pop_thr_arr
      
      idx_pure = sum(nOcc(:,3:5)>=p,2)==1;
      idx_none = sum(nOcc(:,3:5)>=p,2)==0;
      idx_mixed = sum(nOcc(:,3:5)>=p,2)>1;
      
      idx_pop = nOcc>=p;
      
      nPC_only = idx_none & idx_pop(:,2);
      NRNG_only = idx_pure & idx_pop(:,3);
      GT_only = idx_pure & idx_pop(:,4);
      RW_only = idx_pure & idx_pop(:,5);
      
%        nPop_nPC(p) = sum(nPC_only);
%        nPop_GT(p) = sum(GT_only);
%        nPop_RW(p) = sum(RW_only);
%        nPop_NRNG(p) = sum(NRNG_only);
%        nPop_multi(p) = sum(idx_mixed);
      
      nPop_nPC(p) = sum(idx_pop(:,2));
      nPop_GT(p) = sum(idx_pop(:,4));
      nPop_RW(p) = sum(idx_pop(:,5));
      nPop_NRNG(p) = sum(idx_pop(:,3));
      nPop_multi(p) = sum(idx_mixed);
      
      nPop_rnd(p) = nchoosek_sum(nSes,p,pTotal(3));
    end
    
    figure('position',[500 500 400 200])
    
    hold on
%      plot(nPop_nPC/nPop_nPC(1),'k')
    plot(nPop_NRNG/nPop_NRNG(1),'b','DisplayName',sprintf('NRNG (p = %3.1f%%)',100*pTotal(3)))
    plot(nPop_GT/nPop_GT(1),'r','DisplayName',sprintf('GT (p = %3.1f%%)',100*pTotal(4)))
    plot(nPop_RW/nPop_RW(1),'g','DisplayName',sprintf('RW (p = %3.1f%%)',100*pTotal(5)))
%      plot(nPop_multi/nPop_multi(1),'m')
    plot(nPop_rnd/nPop_rnd(1),'k--','DisplayName',sprintf('expectation (p=%3.1f%%)',100*pTotal(3)))
    ylim([0,1.1])
    xlabel('threshold to population \Theta')
    ylabel('N_p(\Theta) / N_p(\Theta=1)')
    legend()
    
    if sv
      path = pathcat(pathFigures,sprintf('population_overoccurence_%s.png',sv_suffix));
      print(path,'-dpng','-r300')
    end
    
%      figure
%      plot((nPop_GT/nPop_GT(1))./(nPop_rnd/nPop_rnd(1)))
%      set(gca,'yscale','log')
  end
  
  
  disp('plot rates of: formation, stabilization, disrecruitment for gate cells & reward cells (from silent, nPC, PC, each)')
  
  if plot_fig(28)
    
    ds = 1;
    
    remap_ds = squeeze(remap_pos(:,ds,:,:));
%      remap_gt = zeros(nSes,para.nbin+2);
%      remap_rw = zeros(nSes,para.nbin+2);
%      remap_PC = zeros(nSes,para.nbin+2);
    
%      for s = 1:s_end
%        remap_gt(s,:) = sum(remap_ds(s,para.idx_zone{4},:),2);
%        remap_rw(s,:) = sum(remap_ds(s,para.idx_zone{5},:),2);
%        remap_PC(s,:) = sum(remap_ds(s,para.idx_zone{3},:),2);
%      end
    
    
%      figure('position',[500 100 1200 900])
    
    ylims = [0.06,0.06,0.02,0.01,0.005];
    
    gt_bar1 = gt_bar*ylims(1);
    rw_bar1 = rw_bar*ylims(2);
    PC_bar1 = PC_bar*ylims(3);
    
    parts = 14;
    
    p_in = zeros(parts,5,5);
    p_out = zeros(parts,5,5);
    
    p_in_r = zeros(parts,5,5);
    p_out_r = zeros(parts,5,5);
    
    
    for j = 1:parts
      idx_state = j;%(2*j-1):min(2*j,nSes);
      title_str = '';
%        if j > 1
        idx_state_in = j;
%        else
%          idx_state_in = [];
%        end
      idx_state_out = j;
      idx_state_norm_in = idx_state_in+1;
      idx_state_norm_out = idx_state_out;
      
      remap_ds_in = squeeze(sum(remap_ds(idx_state_in,:,:),1));
      remap_ds_out = squeeze(sum(remap_ds(idx_state_out,:,:),1));
      
      for i = 1:5
        for k = 1:5
          p_in(j,i,k) = sum(sum(remap_ds_in(para.idx_zone{i},para.idx_zone{k}),2)/sum(sum(occupancy(idx_state_norm_in,para.idx_zone{k}))));
          p_out(j,k,i) = sum(sum(remap_ds_out(para.idx_zone{k},para.idx_zone{i}),1)/sum(sum(occupancy(idx_state_norm_out,para.idx_zone{k}))));
          
          p_in_r(j,i,k) = sum(sum(remap_ds_in(para.idx_zone{i},para.idx_zone{k}),2)/sum(sum(occupancy(idx_state_in,para.idx_zone{i}))));
%            p_out_r(j,k,i) = sum(sum(remap_ds_out(para.idx_zone{k},para.idx_zone{i}),1)/sum(sum(occupancy(idx_state_norm_out,para.idx_zone{k}))));
        end
      end
    end
    occupancy
    col = {'k--','k','b','g','r'};
    
    
    figure('position',[100 100 900 350])
    
    t_ses
    title_arr = {'NRNG','GT','RW'}
    for p = 3:5
      
      idxes = 1:5;
      idxes(p) = [];
      
%        suptitle(sprintf('population %d',p))
      subplot(2,3,p-2)
      title(title_arr{p-2})
      hold on
%        for i = 1:5
%          plot(p_in(:,i,p),col{i})
%        end
      plot(t_ses(2:end),p_in(:,p,p),col{p},'LineWidth',3)
%        plot(sum(p_in(:,idxes,p),2),'m','LineWidth',3)
      
%        plot(sum(p_in(:,1:2,p),2),'k','LineWidth',3)
      size(t_ses(2:end))
      plot(t_ses(2:end),p_in(:,idxes(3),p),col{idxes(3)},'LineWidth',3)
      plot(t_ses(2:end),p_in(:,idxes(4),p),col{idxes(4)},'LineWidth',3)
      xlim([0,t_ses(end)])
      ylim([0,0.5])
      
      if p==3
        ylabel('N_{i\leftarrow j} / N_i','fontsize',14)
      end
      
      subplot(2,3,p-2+3)
%        title('in (r population)')
      hold on
%        for i = 1:6
%          plot(p_in(:,i,p),col{i})
%        end
      plot(t_ses(2:end),p_in_r(:,p,p),col{p},'LineWidth',3)
%        plot(sum(p_in(:,idxes,p),2),'m','LineWidth',3)
      
      plot(t_ses(2:end),sum(p_in_r(:,1:2,p),2),'k','LineWidth',3)
      plot(t_ses(2:end),p_in_r(:,idxes(3),p),col{idxes(3)},'LineWidth',3)
      plot(t_ses(2:end),p_in_r(:,idxes(4),p),col{idxes(4)},'LineWidth',3)
      ylim([0,0.5])
      xlim([0,t_ses(end)])
      xlabel('t [h]')
      if p==3
        ylabel('N_{i\leftarrow j} / N_j','fontsize',14)
      end
      
      
      
      
%        subplot(3,1,2)
%        title('out')
%        hold on
%  %        for i = 1:5
%  %          plot(p_out(:,p,i),col{i})
%  %        end
%        plot(p_out(:,p,p),col{p},'LineWidth',3)
%        plot(sum(p_in(:,1:2,p),2),'k','LineWidth',3)
%        plot(p_in(:,idxes(3),p),col{idxes(3)},'LineWidth',3)
%        plot(p_in(:,idxes(4),p),col{idxes(4)},'LineWidth',3)
%  %        plot(sum(p_out(:,p,idxes),3),'m','LineWidth',3)
%        xlim([1,14])
      
%        subplot(3,1,3)
%        title('diff')
%        hold on
%        plot([1,parts],[0 0],'k:')
%  %        for i = 1:5
%  %          plot(p_in(:,i,p) - p_out(:,p,i),col{i})
%  %        end
%        plot(p_in(:,p,p)-p_out(:,p,p),col{p},'LineWidth',3)
%        plot(sum(p_in(:,idxes,p),2)-sum(p_out(:,p,idxes),3),'m','LineWidth',3)
%        xlim([1,14])
    end
    
    
    if sv
      path = pathcat(pathFigures,sprintf('recruitment_rates%s.png',sv_suffix));
      print(path,'-dpng','-r300')
      disp(sprintf('Figure saved as %s',path))
    end
    
  end
  
  if plot_fig(29)
    disp('what is a GT/RW/NRNG cell coding for after non-coding / silence?')
    disp('doesnt seem to be much of a difference, if non-coding state is in between...')
    
    figure
    for p = 1:3
      recode_hist = zeros(para.nbin,2);
      for c = find(idx_pop(:,p+2))'
        for s=1:s_end-2
          if PCs.status(c,s,p+2) && any(PCs.status(c,s+1,3:5)) && any(PCs.status(c,s+2,3:5))
            for f = find(~isnan(PCs.fields.center(c,s+2,:)))'
              pos = PCs.fields.center(c,s+2,f);
              recode_hist(pos,1) = recode_hist(pos,1) + 1;
            end
          end
          if PCs.status(c,s,p+2) && ~any(PCs.status(c,s+1,3:5))
            
            s_recode = find(any(PCs.status(c,s+1:end,3:5),3),1,'first') + s;
            for f = find(~isnan(PCs.fields.center(c,s_recode,:)))'
              pos = PCs.fields.center(c,s_recode,f);
              recode_hist(pos,2) = recode_hist(pos,2) + 1;
            end
          end
        end
      end
      
      subplot(3,1,p)
      hold on
      bar(recode_hist(:,1),'b','FaceAlpha',0.5)
      bar(-recode_hist(:,2),'r','FaceAlpha',0.5)
%        bar(recode_hist(:,1)-recode_hist(:,2),'k')
    end
    
  end
  
  if plot_fig(30)
    
    %% plot stability with estimate from gaussian-constant fit
    x_data = linspace(-40,40,para.nbin+1);
    
    mu_0 = 0;
    std_0 = 5;
    offset_0 = 0;
    ampl_0 = 1;
    
    fit_shifts = struct;
    fit_shifts.init = [mu_0 std_0 offset_0 ampl_0];
    fit_shifts.lb = [-para.nbin/2 0 0 0];                     %% lower bounds
    fit_shifts.ub = [para.nbin/2, para.nbin/2, 1, 1];                  %% upper bounds
    
    if plot_pop
      fig = figure('position',[100 100 400 200]);
      ax1 = axes('position',[0.2,0.3,0.4,0.65]);
      hold on
      ax2 = axes('position',[0.65,0.3,0.175,0.65]);
      hold on
      for p = 1:3
        
        N_data = length(arrays.ICI_s{p});
        if N_data > 0
          bs_idx = randi(N_data,N_bs,N_data);
          
          p_shifts = zeros(nT,N_bs,4);
          tic
          for L=1:N_bs
              
            ICI_L = arrays.ICI_s{p}(bs_idx(L,:));
            shift_L = arrays.shift{p}(bs_idx(L,:));
            N_ds = zeros(nT,1);
            for ds = 1:nT
              idxes_ds = ICI_L==t_data(ds);
              
              if sum(idxes_ds) > 20
                N_ds(ds) = sum(idxes_ds);
                shift_hist = histcounts(shift_L(idxes_ds),x_data-0.5,'Normalization','probability');
                [p_shifts(ds,L,:),resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F_est,fit_shifts.init,x_data(1:end-1),shift_hist,fit_shifts.lb,fit_shifts.ub,options);    %% fit data
              end
            end
            ampl = p_shifts(:,L,4);
            mask = ampl>0;
            
            X = [ones(size(ampl,1),1),t_data'];
            Y = log(ampl);
            
            %% jackknifing the whole thing
            N_data_tmp = sum(mask);
            p_decay_tmp = zeros(N_data_tmp,2);
            j = 1;
            for i=find(mask)'
              
              mask_ampl = mask;
              mask_ampl(i) = false;
              
              [p_tmp,STDX] = lscov(X(mask_ampl,:),Y(mask_ampl),N_ds(mask_ampl));
              p_decay_tmp(j,1) = -1/p_tmp(2);
              p_decay_tmp(j,2) = exp(p_tmp(1));
              j=j+1;
            end
            p_decay(L,:) = median(p_decay_tmp,1);
          end
          toc
          
          ampl = median(p_shifts(:,:,4),2)';
          ampl_CI = squeeze(prctile(p_shifts(:,:,4),[5,95],2));
          
          tau_r = median(p_decay(:,1))
          tau_CI = prctile(p_decay(:,1),[5,95])
          
          mask = ampl>0;
          
          ampl_CI(ampl_CI(mask,1) == 0,1) = 0.1;
          x2 = [t_data(mask) fliplr(t_data(mask))];
          inBetween = [ampl_CI(mask,1)', flipud(ampl_CI(mask,2))'];
          fill(ax1,x2,inBetween,col_fill{p},'FaceAlpha',0.5,'EdgeColor','None','HandleVisibility','off');
          
          plot(ax1,t_data(mask),ampl(mask),'o','Color',col{p})
          plot(ax1,t_data,exp_dist(median(p_decay,1),t_data),'--','Color',col{p})
          
          errorbar(ax2,p,tau_r,abs(tau_CI(1)-tau_r),abs(tau_CI(2)-tau_r),'o-','Color',col{p})
        end
      end
      set(ax2,'YAxisLocation','right')
      xlim(ax2,[0.5,3.5])
      ylim(ax2,[0,20])
      set(ax2,'XTick',linspace(1,3,3),'XTickLabels',{'NRNG','GT','RW'},'FontSize',12)
      xtickangle(ax2,60)
      ylabel(ax2,'\tau_r [sessions]','FontSize',14)
    else
      fig_pos = [100 100 400 300];
      fig = figure('position',fig_pos);
%        ax1 = axes('position',[0.25,0.3,0.7,0.65]);
      ax1 = subplot(2,1,2);
      hold on
      ax2 = subplot(2,1,1);
      hold on
      
      ICI_data = [arrays.ICI_s{:}];
      shift_data = [arrays.shift{:}];
      N_data = length(ICI_data);
      bs_idx = randi(N_data,N_bs,N_data);
      X = [ones(nT,1),t_data'];
      p_shifts = zeros(nT,N_bs,4);
      tic
      for L=1:N_bs
        
        ICI_L = ICI_data(bs_idx(L,:));
        shift_L = shift_data(bs_idx(L,:));
        N_ds = zeros(nT,1);
        for ds = 1:nT
          idxes_ds = ICI_L==t_data(ds);
          if sum(idxes_ds) > 20
            N_ds(ds) = sum(idxes_ds);
            shift_hist = histcounts(shift_L(idxes_ds),x_data-0.5,'Normalization','probability');
            [p_shifts(ds,L,:),resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F_est,fit_shifts.init,x_data(1:end-1),shift_hist,fit_shifts.lb,fit_shifts.ub,options);    %% fit data
          end
        end
        ampl = p_shifts(:,L,4);
        mask = ampl>0 & (1:nSes <= 10)';
        Y = log(ampl);
        
        %% jackknifeing the whole thing
        N_data_tmp = length(ampl);
        p_decay_tmp = zeros(N_data_tmp,2)*NaN;
        j = 1;
        for i=find(mask)'
          mask_ampl = mask;
          mask_ampl(i) = false;
          
%            tbl = table(t_data', ampl);
%            tbl
%            modelfun = @(beta,x) beta(1) * x(:, 1) .^ + beta(2) + beta(3);  
%            beta_1 = 1;
%            beta_2 = -1;
%            beta_3 = 0;
%            beta0 = [beta_1, beta_2, beta_3]; % Guess values to start with.  Just make your best guess.
%            % Now the next line is where the actual model computation is done.
%            mdl = fitnlm(tbl, modelfun, beta0);
%            coefficients = mdl.Coefficients{:, 'Estimate'}
%            modelfun(coefficients,[t_data',ampl])
%            figure
%            hold on
%            plot(t_data,ampl,'ok')
%            plot(t_data,modelfun(coefficients,[t_data',ampl]),'r')
%            
%            waitforbuttonpress
          
          p_tmp = lscov(X(mask_ampl,:),Y(mask_ampl),N_ds(mask_ampl));
          p_decay_tmp(j,1) = -1/p_tmp(2);
          p_decay_tmp(j,2) = exp(p_tmp(1));
          j=j+1;
        end
        p_decay(L,:) = nanmedian(p_decay_tmp,1);
      end
      toc
      
      fig1 = figure('position',[500 500 600 200]);
      k = 1;
      
      for ds = 1:nT
        idxes_ds = ICI_data==t_data(ds);
        if sum(idxes_ds) > 20
          shift_hist = histcounts(shift_data(idxes_ds),x_data-0.5,'Normalization','probability');
          p_shifts_data(1,ds,:) = lsqcurvefit(F_est,fit_shifts.init,x_data(1:end-1),shift_hist,fit_shifts.lb,fit_shifts.ub,options);    %% fit data
        
        
          if ismember(t_data(ds),[t_data(1),t_data(3),t_data(5)])
            subplot(1,3,k)
            hold on
            bar(x_data(1:end-1),shift_hist,'FaceColor','k','DisplayName','field shifts \Delta x')
            plot(x_data,F_est(p_shifts_data(1,ds,:),x_data),'r:','LineWidth',2,'HandleVisibility','off')
            ylim([0,0.14])
            set(gca,'FontSize',10)
            set(gca,'XTick',linspace(1,para.nbin,5),'XTickLabel',linspace(0,100,5))
            xlabel('\Delta x [cm]','FontSize',14)
            title(['\Delta ',sprintf('s = %d',t_data(ds))])
            if k==1
              ylabel('p(\Delta x)','FontSize',14)
            end
            if k == 3
              legend('Location','NorthWest');
            end
            k = k+1;
          end
        end
      end
      
      if sv
        figure(fig1)
        save_fig(pathFigures,mouse,'shift_distr',sv_suffix,sv_ext,[500 500 600 200])
      end
      
      ampl = median(p_shifts(:,:,4),2)
      ampl_CI = squeeze(prctile(p_shifts(:,:,4),[2.5,97.5],2))
      ampl_CI(ampl_CI <=0) = 0.001
      ampl(ampl<=0) = NaN;
      ampl_CI(ampl<=0) = NaN;
      
      
      sig = median(p_shifts(:,:,2),2)*100/para.nbin;
      sig_CI = squeeze(prctile(p_shifts(:,:,2),[2.5,97.5],2))*100/para.nbin;
      
      tau_r = median(p_decay(:,1));
      tau_CI = prctile(p_decay(:,1),[5,95]);
      
      plot_CI_as_fill(ampl,ampl_CI,t_data,ax1,{'k',[0.5,0.5,0.5]},[])
      plot(ax1,t_data,exp_dist(median(p_decay,1),t_data),'--','Color','k','DisplayName',['\tau_r',sprintf('= %4.2g, [%4.2g,%4.2g]',tau_r,tau_CI(1),tau_CI(2))])%'A_0 exp(-\Delta s / \tau_r)')
      legend(ax1,'Location',[0.52,0.42,0.35,0.1])
      
      plot(ax2,[0,t_ses(end)],[3.5,3.5],'k:')
      plot_CI_as_fill(sig,sig_CI,t_data,ax2,{'k',[0.5,0.5,0.5]},[])
      
      ylabel(ax2,'\sigma [cm]')
      set(ax2,'FontSize',12)
      xlim(ax2,[0,nSes])
      ylim(ax2,[0,10])
      
    end
    set(ax1,'Yscale','log')
    set(ax1,'FontSize',12)
%      xlim(ax1,[0,t_ses(end)])
    xlim(ax1,[0,nSes])
    ylim(ax1,[0.05,1])
    xlabel(ax1,'\Delta s','FontSize',14)
    ylabel(ax1,'r_{stable}','FontSize',14)    
    
%      
%      figure()
%      theta = median(p_shifts(:,:,1),2)*100/para.nbin;
%      theta_CI = squeeze(prctile(p_shifts(:,:,1),[2.5,97.5],2))*100/para.nbin;
%  %      theta(theta==0) = NaN
%      
%      plot_CI_as_fill(theta,theta_CI,t_data,gca,{'k',[0.5,0.5,0.5]},[])
%      
%      ylim([-50,50])
    
    if sv
      figure(fig)
      save_fig(pathFigures,mouse,'stability_fit',sv_suffix,sv_ext,fig_pos)
    end
  end
  
  
  
  if plot_fig(31)
    %% plot waiting time distribution between coding sessions
    
    X = [ones(nT,1),t_data'];
    
    if plot_pop
      
      tau_w = zeros(3,1);
      tau_CI = zeros(3,2);
      
      figure('position',[100 100 400 200]);
      ax1 = axes('position',[0.2,0.3,0.4,0.65]);
      hold on
      ax2 = axes('position',[0.65,0.3,0.175,0.65]);
      hold on
      
      for p=1:3
        
        idx_discont = arrays.r_coding{p} == 0;% | arrays.ICI_s{p}==1;
        ICI_data = arrays.ICI_s{p}(idx_discont);
        
        N_data = length(ICI_data);
        
        if N_data > 0
          bs_idx = randi(N_data,N_bs,N_data);
          
          tic
          p_decay = zeros(N_bs,2);
          hist_ICI = zeros(N_bs,nT);
          for L = 1:N_bs
            data_L = ICI_data(bs_idx(L,:));
            hist_ICI_tmp = histcounts(data_L,[t_data-0.5,t_data(end)+0.5],'Normalization','probability');
            hist_ICI(L,:) = hist_ICI_tmp(t_mask);
            
            p_tmp = expfit_jackknife(X,hist_ICI(L,:)');
            p_decay(L,1) = median(-1./p_tmp(:,2));
            p_decay(L,2) = median(exp(p_tmp(:,1)));
          end
          toc
          
          ICI_med = median(hist_ICI,1);
          ICI_CI = prctile(hist_ICI,[5,95],1)';
          ICI_CI(ICI_CI==0) = 10^(-3);
          
          tau_w(p) = median(p_decay(:,1));
          tau_CI(p,:) = prctile(p_decay(:,1),[5,95]);
          mask = ICI_med > 0;
          
          x2 = [t_data(mask) fliplr(t_data(mask))];
          inBetween = [ICI_CI(mask,1)', flipud(ICI_CI(mask,2))'];
          fill(ax1,x2,inBetween,col_fill{p},'FaceAlpha',0.5,'EdgeColor','None','HandleVisibility','off');
        
          plot(ax1,t_data,exp(log_exp_dist(median(p_decay,1),t_data)),'--','Color',col{p},'LineWidth',1,'HandleVisibility','off')
          plot(ax1,t_data,ICI_med,'o','Color',col{p},'DisplayName',plot_arr{p})
          
          errorbar(ax2,p,tau_w(p),abs(tau_CI(p,1)-tau_w(p)),abs(tau_CI(p,2)-tau_w(p)),'o-','Color',col{p})
        end
      end
      set(ax2,'YAxisLocation','right')
      ylabel(ax2,'\tau_w [sessions]','FontSize',14)
      xlim(ax2,[0.5,3.5])
      ylim(ax2,[0,3.5])
      set(ax2,'XTick',linspace(1,3,3),'XTickLabels',{'NRNG','GT','RW'},'FontSize',12)
      xtickangle(ax2,60)
    
    else
      figure('position',[100 100 300 150])
      ax1 = axes('position',[0.25,0.32,0.7,0.6]);
      hold on
%        [p_out,resnorm,resid,~,~,~,J] = lsqcurvefit(exp_dist,p_0,t_data,sum(ICI,2)/sum(ICI(:)),lb,ub,options);
      
      idx_discont = [arrays.r_coding{:}] == 0;% | [arrays.ICI_s{:}] == 1;
      ICI_data = [arrays.ICI_s{:}];
      ICI_data = ICI_data(idx_discont);
      
      N_data = length(ICI_data);
      bs_idx = randi(N_data,N_bs,N_data);
      
      p_decay = zeros(N_bs,2);
      hist_ICI = zeros(N_bs,nT);
      for L = 1:N_bs
        hist_ICI_tmp = histcounts(ICI_data(bs_idx(L,:)),[t_data-0.5,t_data(end)+0.5],'Normalization','probability');
        hist_ICI(L,:) = hist_ICI_tmp(t_mask);
        mask = hist_ICI > 0;
        
        p_tmp = expfit_jackknife(X,hist_ICI(L,:)');
        p_decay(L,1) = median(-1./p_tmp(:,2));
        p_decay(L,2) = median(exp(p_tmp(:,1)));
        
      end
      
      ICI_med = median(hist_ICI,1);
      ICI_CI = prctile(hist_ICI,[5,95],1)';
      ICI_CI(ICI_CI==0) = 10^(-3);
      
      mask = ICI_med > 0;
      tau_w = median(p_decay(:,1));
      tau_CI = prctile(p_decay(:,1),[5,95]);
        
      
      x2 = [t_data(mask) fliplr(t_data(mask))];
      inBetween = [ICI_CI(mask,1)', flipud(ICI_CI(mask,2))'];
      fill(ax1,x2,inBetween,[0.5,0.5,0.5],'FaceAlpha',0.3,'EdgeColor','None','HandleVisibility','off');
      
      plot(ax1,t_data,ICI_med,'o','Color','k','HandleVisibility','off')
      plot(ax1,t_data,exp(log_exp_dist(median(p_decay,1),t_data)),'k--','DisplayName',['\tau_w',sprintf('= %4.2g, [%4.2g,%4.2g]',tau_w,tau_CI(1),tau_CI(2))])%'A_0 exp(-\Delta s_w/\tau_w)')
      legend(ax1,'Location',[0.52,0.85,0.35,0.1])
      
    end
    set(ax1,'YScale','log')
    set(ax1,'FontSize',12)
    xlim(ax1,[0,t_ses(end)])
    xlabel(ax1,'waiting time \Delta s_w','FontSize',14)
    ylim(ax1,[10^(-2),1])
    ylabel(ax1,'p(\Delta s_w)','FontSize',14)    
    
    if sv
      path = pathcat(pathFigures,sprintf('waiting_times%s.%s',sv_suffix,sv_ext));
      print(path,sprintf('-d%s',sv_ext),'-r300')
      disp(sprintf('Figure saved as %s',path))
    end
    
  end
    
    
    
    
  if plot_fig(32)
    %%% plot discont vs cont
    
    cdf_discont = zeros(nT,para.nbin/2+1)*NaN;
    cdf_cont = zeros(nT,para.nbin/2+1)*NaN;
    
    ICI_data = [arrays.ICI_s{:}];
    shift_data = [arrays.shift{:}];
    
    idx_discont = [arrays.r_coding{:}] == 0;
    ICI_discont = ICI_data(idx_discont);
    shift_discont = shift_data(idx_discont);
    
    idx_cont = [arrays.r_coding{:}] == 1;
    ICI_cont = ICI_data(idx_cont);
    shift_cont = shift_data(idx_cont);
    
    for ds = 1:nT
      %% discont
      idx_ds = ICI_discont==t_data(ds);
      if sum(idx_ds) > 20
        [cdf_tmp,x] = ecdf(round(abs(shift_discont(idx_ds))));
        cdf_discont(ds,x(2:end)+1) = cdf_tmp(2:end);
      end
      
      %% cont
      idx_ds = ICI_cont==t_data(ds);
      if sum(idx_ds) > 20
        [cdf_tmp,x] = ecdf(round(abs(shift_cont(idx_ds))));
        cdf_cont(ds,x(2:end)+1) = cdf_tmp(2:end);
      end
    end
    
    x_arr = linspace(0,para.nbin/2,para.nbin/2+1);
    x_arr2 = linspace(-para.nbin/2,para.nbin/2,para.nbin+1);
    D_cont = zeros(nT,1)*NaN;
    D_dc = zeros(nT,1)*NaN;
    
    cdf_uniform = linspace(0,1,para.nbin/2+1);
    D_max = max(cdf_cont(1,:) - cdf_uniform);
    figure('position',[100 100 400 250])
    
    ax(1) = subplot(1,1,1);
%      ax(2) = subplot(1,2,2);
    hold(ax(1),'on')
%      hold(ax(2),'on')
    
    plot([7,7],[0,1],'k:','HandleVisibility','off')
    plot([0,0],[0,0],'k-','LineWidth',2,'DisplayName',['\Delta s=',sprintf('%d',t_data(1))])
    plot(x_arr,cdf_uniform,'k--','HandleVisibility','off')
    for ds = 2:nT
      [~,pos_max] = max(abs(cdf_cont(1,:)-cdf_cont(ds,:)));
      D_cont(ds) = cdf_cont(1,pos_max)-cdf_cont(ds,pos_max);
      
      [~,pos_max] = max(abs(cdf_cont(1,:)-cdf_discont(ds,:)));
      D_dc(ds) = cdf_cont(1,pos_max)-cdf_discont(ds,pos_max);
      
      mask_disc = ~isnan(cdf_discont(ds,:));
      mask_cont = ~isnan(cdf_cont(ds,:));
      if ismember(t_data(ds),[t_data(2),t_data(4),t_data(6)])%mod(ds,2) && ds < 12
        col0 = ds/10;
        stairs(ax(1),x_arr(mask_cont),cdf_cont(ds,mask_cont),'Color',[col0,col0,1],'LineWidth',2,'HandleVisibility','off')
        stairs(ax(1),x_arr(mask_disc),cdf_discont(ds,mask_disc),'Color',[1,col0,col0],'DisplayName',['\Delta s=',sprintf('%d',t_data(ds))],'LineWidth',2)
        
%          if ismember(t_data(ds),[24,72])
%            text(ax(1),x_arr(round(para.nbin/2*0.3)),cdf_discont(ds,round(para.nbin/2*0.3))+0.05,,'FontSize',10)
%          end
%          r_stable = 0.8-0.1*(ds-1)
%          p_theory = F_est([0,2,1-r_stable,r_stable],x_arr2);
%          p_theory(para.nbin/2+2:end) = p_theory(para.nbin/2+2:end) + fliplr(p_theory(1:para.nbin/2));   %% take absolute
%          p_theory(1:para.nbin/2) = [];                                                        %% remove negative side of distribution
%          plot(ax(2),x_arr2(para.nbin/2+1:end),cumsum(p_theory),'r','LineWidth',0.5,'DisplayName','Model')
      end
    end
    mask_cont = ~isnan(cdf_cont(1,:));
    stairs(ax(1),x_arr(mask_cont),cdf_cont(1,mask_cont),'Color','k','LineWidth',2,'HandleVisibility','off')
%      for i=1:2
    xlim(ax(1),[0,40])
    ylim(ax(1),[0,1])
    xlabel(ax(1),'|\Delta x| [cm]','FontSize',14)
    set(ax(1),'XTick',linspace(0,40,3),'XTickLabels',linspace(0,50,3))
%      end
    ylabel(ax(1),'cdf(|\Delta x|)','FontSize',14)
    lgd = legend(ax(1),'Location','South');
    lgd.FontSize = 8;
    
    x_arr = linspace(-para.nbin/2,para.nbin/2,para.nbin+1);
    r_stable = 0.7;
    p_theory = F_est([0,2,1-r_stable,r_stable],x_arr);
    p_theory(para.nbin/2+2:end) = p_theory(para.nbin/2+2:end) + fliplr(p_theory(1:para.nbin/2));   %% take absolute
    p_theory(1:para.nbin/2) = [];                                                        %% remove negative side of distribution
%      plot(ax(1),x_arr(para.nbin/2+1:end),cumsum(p_theory),'g--','LineWidth',2,'DisplayName','Model')          
    
    ax_zoom = axes('position',[0.75,0.45,0.2,0.2]);
    hold on
    mask_disc = ~isnan(D_dc);
    mask_cont = ~isnan(D_cont);
    plot([0,nSes],[0,0],'k:')
    plot([0,nSes],[D_max,D_max],'k--')
    plot(ax_zoom,t_data(mask_cont),D_cont(mask_cont),'b','DisplayName','continuous')
    plot(ax_zoom,t_data(mask_disc),D_dc(mask_disc),'r','DisplayName','discontinuous')
    xlim([0,20])
%      xlim([0,t_ses(end)])
    ylim([-0.3,0.5])
    xlabel('\Delta s','FontSize',10)
    ylabel('D_{KS}','FontSize',10)
    
    if sv
      path = pathcat(pathFigures,sprintf('contvsdiscont2%s.%s',sv_suffix,sv_ext));
      print(path,sprintf('-d%s',sv_ext),'-r300')
      disp(sprintf('Figure saved as %s',path))
    end
    
%      waitforbuttonpress
    plt_stuff = true;
    
    if plt_stuff
      options = statset('MaxIter',1000, 'MaxFunEvals',2000,'Display','off');
      
      %% initial guesses
      mu_0 = 0;
      std_0 = 5;
      offset_0 = 0;
      ampl_0 = 1;
      p_init=[mu_0 std_0 offset_0 ampl_0];
      
      lb = [-40 0 0 0];                     %% lower bounds
      ub = [40, 10, 1, 1];                  %% upper bounds
      
      mu = zeros(nT,3,3)*NaN;
      mu_CI = zeros(nT,3,2,3)*NaN;
      stdx = zeros(nT,3,3)*NaN;
      stdx_CI = zeros(nT,3,2,3)*NaN;
      offset = zeros(nT,3,3)*NaN;
      offset_CI = zeros(nT,3,2,3)*NaN;
      ampl = zeros(nT,3,3)*NaN;
      ampl_CI = zeros(nT,3,2,3)*NaN;
      
      
      x_arr = linspace(-40.5,40.5,para.nbin+2);
      x_arr_p = linspace(-40,40,para.nbin+1);
      
  %      mu = zeros(t_ses(end),4,2,2)
      j=0;
      p=1;
%        figure('position',[100 100 800 200])
      
      
      idx_cont = [arrays.r_coding{:}] == 1;
      idx_disc = [arrays.r_coding{:}] < 1 & [arrays.r_coding{:}] > 0;
      idx_disc_0 = [arrays.r_coding{:}] == 0;
      
      ICI_data = [arrays.ICI_s{:}];
      shift_data = [arrays.shift{:}];
      
      ICI_cont = ICI_data(idx_cont);
      ICI_disc = ICI_data(idx_disc);
      ICI_disc_0 = ICI_data(idx_disc_0);
      
      shift_cont = shift_data(idx_cont);
      shift_disc = shift_data(idx_disc);
      shift_disc_0 = shift_data(idx_disc_0);
      
      N_data_cont = length(ICI_cont);
      N_data_disc = length(ICI_disc);
      N_data_disc_0 = length(ICI_disc_0);
      
      p_final_cont = zeros(nT,N_bs,4)*NaN;
      p_final_disc = zeros(nT,N_bs,4)*NaN;
      p_final_disc_0 = zeros(nT,N_bs,4)*NaN;
      
      for i=2:6
        dt = t_data(i);
        
        bs_idx_cont = randi(N_data_cont,N_bs,N_data_cont);
        bs_idx_disc = randi(N_data_disc,N_bs,N_data_disc);
        bs_idx_disc_0 = randi(N_data_disc_0,N_bs,N_data_disc_0);
        
        for L=1:N_bs
          
          ICI_cont_L = ICI_cont(bs_idx_cont(L,:));
          shift_cont_L = shift_cont(bs_idx_cont(L,:));
          idx_cont_ds = ICI_cont_L==t_data(i);
          N_cont = sum(idx_cont_ds);
          if N_cont > 10
            hist_cont = histcounts(shift_cont_L(idx_cont_ds),x_arr,'Normalization','probability');
            [p_final_cont(i,L,:),resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F_est,p_init,x_arr_p,hist_cont,lb,ub,options);    %% fit data
  %            [p_final,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F_est,p_init,x_arr_p,hist_cont,lb,ub,options);    %% fit data
%              ci = nlparci(p_final,resid,'jacobian',J);
            
%              mu_CI(i,p,:,1) = ci(1,:);
%              stdx_CI(i,p,:,1) = ci(2,:);
%              offset_CI(i,p,:,1) = ci(3,:);
%              ampl_CI(i,p,:,1) = ci(4,:);
%              
%              mu(i,p,1) = p_final(1);
%              stdx(i,p,1) = p_final(2);
%              offset(i,p,1) = p_final(3);
%              ampl(i,p,1) = p_final(4);
          end
          
          ICI_disc_L = ICI_disc(bs_idx_disc(L,:));
          shift_disc_L = shift_disc(bs_idx_disc(L,:));
          idx_disc_ds = ICI_disc_L==t_data(i);
          N_discont = sum(idx_disc_ds);
          if N_discont > 10
            hist_disc = histcounts(shift_disc_L(idx_disc_ds),x_arr,'Normalization','probability');
            [p_final_disc(i,L,:),resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F_est,p_init,x_arr_p,hist_disc,lb,ub,options);    %% fit data
%              [p_final2,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F_est,p_init,x_arr_p,hist_discont,lb,ub,options);    %% fit data
%              ci = nlparci(p_final2,resid,'jacobian',J);
            
%              mu_CI(i,p,:,2) = ci(1,:);
%              stdx_CI(i,p,:,2) = ci(2,:);
%              offset_CI(i,p,:,2) = ci(3,:);
%              ampl_CI(i,p,:,2) = ci(4,:);
%              
%              mu(i,p,2) = p_final2(1);
%              stdx(i,p,2) = p_final2(2);
%              offset(i,p,2) = p_final2(3);
%              ampl(i,p,2) = p_final2(4);
          end
          
          ICI_disc_0_L = ICI_disc_0(bs_idx_disc_0(L,:));
          shift_disc_0_L = shift_disc_0(bs_idx_disc_0(L,:));
          idx_disc_0_ds = ICI_disc_0_L==t_data(i);
          N_discont_0 = sum(idx_disc_0_ds);
          if N_discont_0 > 10
            hist_disc_0 = histcounts(shift_disc_0_L(idx_disc_0_ds),x_arr,'Normalization','probability');
            [p_final_disc_0(i,L,:),resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F_est,p_init,x_arr_p,hist_disc_0,lb,ub,options);    %% fit data
%              [p_final2,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(F_est,p_init,x_arr_p,hist_discont_0,lb,ub,options);    %% fit data
%              ci = nlparci(p_final2,resid,'jacobian',J);
            
%              mu_CI(i,p,:,3) = ci(1,:);
%              stdx_CI(i,p,:,3) = ci(2,:);
%              offset_CI(i,p,:,3) = ci(3,:);
%              ampl_CI(i,p,:,3) = ci(4,:);
%              
%              mu(i,p,3) = p_final2(1);
%              stdx(i,p,3) = p_final2(2);
%              offset(i,p,3) = p_final2(3);
%              ampl(i,p,3) = p_final2(4);
          end
        
        end
        
        
%          if ismember(t_data(i),[t_data(2),t_data(4),t_data(6)])
%            j=j+1;
%            subplot(1,3,j)
%            hold on
%            if N_cont > 20
%              bar(x_arr_p,hist_cont,'FaceColor','b','FaceAlpha',1);%,'DisplayName',sprintf('N_{cont} = %d',N_cont))
%            end
%            if N_discont > 20
%              bar(x_arr_p,hist_discont,'FaceColor','r','FaceAlpha',0.6)%;,'DisplayName',sprintf('N_{discont} = %d',N_discont))
%            end
%            if j == 1
%              ylabel('p(\Delta x)')
%            end
%            
%            if j == 3
%              legend({'cont.','disc.'},'Location','NorthEast')
%            end
%  %            plot(x_arr_p,F_est(p_final,x_arr_p),'r--','LineWidth',2)
%            
%            pos = get(gca,'Position');
%            pos(2) = 0.2;
%            pos(4) = 0.65;
%            set(gca,'Position',pos)
%            set(gca,'XTick',linspace(-40,40,3),'XTickLabel',linspace(-50,50,3))
%            xlabel('\Delta x [cm]')
%            title(['\Delta',sprintf('s=%d',dt)])
%            ylim([0,0.1])
%          end
        
      end
      
%        ampl_cont = p_final_cont(:,:,4);
      ampl_cont = median(p_final_cont(:,:,4),2,'omitnan');
      ampl_cont_CI = squeeze(prctile(p_final_cont(:,:,4),[2.5,97.5],2));
      
      ampl_disc = median(p_final_disc(:,:,4),2,'omitnan');
      ampl_disc_CI = squeeze(prctile(p_final_disc(:,:,4),[2.5,97.5],2));
      
      ampl_disc_0 = median(p_final_disc_0(:,:,4),2,'omitnan');
      ampl_disc_0_CI = squeeze(prctile(p_final_disc_0(:,:,4),[2.5,97.5],2));
      
%        ampl_disc = p_final_disc(:,:,4);
%        ampl_disc_0 = p_final_disc_0(:,:,4);
      
      if sv
        path = pathcat(pathFigures,sprintf('contvsdiscont%s.%s',sv_suffix,sv_ext));
        print(path,sprintf('-d%s',sv_ext),'-r300')
        disp(sprintf('Figure saved as %s',path))
      end
      
      fig_pos = [100 100 400 200];
      figure('position',fig_pos)
      
      x_arr = linspace(2,6,5);
      hold on
      
      
      
%        bar(1,ampl(1,1,1),0.3,'k')
      ampl_cont(x_arr)-ampl_cont_CI(x_arr,:)
      bar(x_arr-0.2,ampl_cont(x_arr),0.2,'FaceColor',[0.4,0.4,1],'DisplayName','continuously coding')
      errorbar(x_arr-0.2,ampl_cont(x_arr),ampl_cont(x_arr)-ampl_cont_CI(x_arr,1),ampl_cont_CI(x_arr,2)-ampl_cont(x_arr),'k','LineStyle','None','HandleVisibility','off')
      
      bar(x_arr,ampl_disc(x_arr),0.2,'FaceColor',[1,0.4,0.4],'DisplayName','episodes of non-coding')
      errorbar(x_arr,ampl_disc(x_arr),ampl_disc(x_arr)-ampl_disc_CI(x_arr,1),ampl_disc_CI(x_arr,2)-ampl_disc(x_arr),'k','LineStyle','None','HandleVisibility','off')
      
      bar(x_arr+0.2,ampl_disc_0(x_arr),0.2,'FaceColor',[0.4,0.4,0.4],'DisplayName','non-coding')
      errorbar(x_arr+0.2,ampl_disc_0(x_arr),ampl_disc_0(x_arr)-ampl_disc_0_CI(x_arr,1),ampl_disc_0_CI(x_arr,2)-ampl_disc_0(x_arr),'k','LineStyle','None','HandleVisibility','off')
      
%        squeeze(ampl(:,1,:))
%        squeeze(ampl_CI(:,1,:,:))
%        ampl_CI(ampl_CI<=0) = 0.001;
%        plot_CI_as_fill(ampl(:,1,1),squeeze(ampl_CI(:,1,:,1)),t_data,gca,{'b',[0.5,0.5,1]},'cont.')
%        plot_CI_as_fill(ampl(:,1,2),squeeze(ampl_CI(:,1,:,2)),t_data,gca,{'r',[1,0.5,0.5]},'discont.')
%        set(gca,'yscale','log')
%        ylim([2*10^(-1),1])
      xlabel('session difference \Delta s','Fontsize',12)
%        ylabel('r_{stable}','Fontsize',14)
      ylabel('fraction stable','Fontsize',12)
      set(gca,'YTick',linspace(0,1,3))
      xlim([1.5,6.5])
      ylim([0,1.3])
      
      lgd = legend('Location',[0.6 0.73 0.3 0.17]);
      title(lgd,'intermittent function')
      if sv
        save_fig(pathFigures,mouse,'disc_vs_cont',sv_suffix,sv_ext,fig_pos)
      end
%        subplot(4,3,9)
%        hold on
%        plot([0,t_data(12)],[0.8,0.8],'k:')
%        errorbar(t_data,ampl(:,1,1),ampl_CI(:,1,1,1)-ampl(:,1,1),ampl_CI(:,1,2,1)-ampl(:,1,1),'b')
%        errorbar(t_data,ampl(:,1,2),ampl_CI(:,1,1,2)-ampl(:,1,2),ampl_CI(:,1,2,2)-ampl(:,1,2),'r')
%        ylabel('r_{stable}','FontSize',14)
%        ylim([0,1])
%        
%        subplot(4,3,12)
%        hold on
%        plot([0,t_data(12)],[2,2],'k:')
%        errorbar(t_data,stdx(:,1,1),stdx_CI(:,1,1,1)-stdx(:,1,1),stdx_CI(:,1,2,1)-stdx(:,1,1),'k')
%        errorbar(t_data,stdx(:,1,2),stdx_CI(:,1,1,2)-stdx(:,1,2),stdx_CI(:,1,2,2)-stdx(:,1,2),'r')
%        xlabel('\Delta t [h]','FontSize',14)
%        ylabel('\sigma','FontSize',14)
%        ylim([0,6])
      
      
%        x_arr = linspace(-para.nbin/2,para.nbin/2,para.nbin+1);
%        
%        figure
%        hold on
%        for i = 1:5
%          col = [(i-1)/8,0,0];
%          
%          r_stable = 1-0.2*i;
%          p_theory = F_est([0,2,1-r_stable,r_stable],x_arr);
%          p_theory(para.nbin/2+2:end) = p_theory(para.nbin/2+2:end) + fliplr(p_theory(1:para.nbin/2));   %% take absolute
%          p_theory(1:para.nbin/2) = [];                                                        %% remove negative side of distribution
%          plot(x_arr(para.nbin/2+1:end),cumsum(p_theory),'Color',col,'LineWidth',2,'DisplayName','Model')     
%        end
%        
%    %      mask
%    %      shift_dist_cont
%    %      size(shift_dist_cont)
%        fig1 = figure('position',[100 100 600 200]);
%        fig2 = figure('position',[100 100 600 200]);
%        
%        for ds = 3:5%size(shift_dist_cont,1)
%          figure(fig1)
%          subplot(1,3,ds-2)
%          hold on
%          plot([6,6],[0,1],'k:')
%          
%          cdf_cont = cumsum(shift_dist_cont(ds,:)/sum(shift_dist_cont(ds,:)));
%          cdf_discont = cumsum(shift_dist_discont(ds,:)/sum(shift_dist_discont(ds,:)));
%          stairs(linspace(0,para.nbin-1,para.nbin),cdf_cont,'k')
%          stairs(linspace(0,para.nbin-1,para.nbin),cdf_discont,'r')
%          
%          D = sum(cdf_cont-cdf_discont);
%          disp(sprintf('D: %5.3g',D))
%    %        disp('sums')
%    %        sum(shift_dist_cont(ds,:))
%    %        sum(shift_dist_discont(ds,:))
%    %        shift_dist_cont(ds,:)
%    %        [F_cont, x_cont] = ecdf(shift_dist_cont(ds,:));
%    %        [F_dc, x_dc] = ecdf(shift_dist_discont(ds,:));
%    %        x_cont
%    %        x_dc
%    %        stairs(x_cont*1.25,F_cont,'k','DisplayName','cont.')
%    %        stairs(x_dc*1.25,F_dc,'r','DisplayName','discont.')
%          
%    %        
%          title(['\Delta',sprintf('t=%d',mask(ds))])
%    %        xlim([0,20])
%          ylim([0,1])
%          xlabel('|\Delta f| [cm]','FontSize',14)
%          if ds == 3
%            ylabel('cdf(|\Delta f|)','FontSize',14)
%          end  
%          
%          if ds == 5
%            legend('Location','SouthEast')
%          end
%          
%          figure(fig2)
%          subplot(1,3,ds-2)
%          hold on
%          bar(shift_dist_cont(ds,:)/sum(shift_dist_cont(ds,:)),'FaceColor','k','FaceAlpha',0.5)
%          bar(-shift_dist_discont(ds,:)/sum(shift_dist_discont(ds,:)),'FaceColor','r','FaceAlpha',0.5)
%          xlim([0,20])
%        end
%        
%        if sv
%          path = pathcat(pathFigures,sprintf('contvsdiscont%s.%s',sv_suffix,sv_ext));
%          print(path,sprintf('-d%s',sv_ext),'-r300')
%          disp(sprintf('Figure saved as %s',path))
%        end
    end
    
  end
  
  if plot_fig(33)
    %% time warping in HC
    
    figure('position',[100 100 600 400])
%      ax(1) = subplot(1,1,1);
%      hold on
    for ds = 1:5
      ax(ds) = subplot(5,1,ds);
      hold on
      for dt = 0:10
        plot(ax(ds),[dt*24+14,dt*24+14],[0,1],'k:')
      end
    end
    lbl_arr = {'other cell (NRNG)','gate cell (GT)','reward cell (RW)'};
    stable_dt = zeros(nT,sum(t_mask_m),N_bs,3)*NaN;
    
    
    
    
    for p = 1:3
      
      N_data = length(arrays.ICI_s{p});
      bs_idx = randi(N_data,N_bs,N_data);
        
      tic
      for L = 1:N_bs
        ICI_s_L = arrays.ICI_s{p}(bs_idx(L,:));
        ICI_t_L = arrays.ICI_t{p}(bs_idx(L,:));
        shift_L = arrays.shift{p}(bs_idx(L,:));
        
        for ds = 1:s_end
          idxes_ds = ICI_s_L==ds;
          for dt = 1:nT
            idxes_dt = ICI_t_L==t_data_m(dt);
            if any(idxes_ds&idxes_dt)
              stable_dt(ds,dt,L,p) = sum(abs(shift_L(idxes_ds&idxes_dt)) < tolerance)/N_norm(ds,dt,p);
            end
          end
        end
      end
      toc
      stable_mean = mean(stable_dt(:,:,:,p),3);
      stable_CI = prctile(stable_dt(:,:,:,p),[5,95],3);
      
      
      for ds=1:5
%          subplot(1,5,ds)
%          plot(stable_mean(ds,:),'o','Color',col{p})
        errorbar(ax(ds),t_data_m,stable_mean(ds,:),abs(stable_CI(ds,:,1)-stable_mean(ds,:)),abs(stable_CI(ds,:,2)-stable_mean(ds,:)),'o-','Color',col{p})
      end
%        for dt=1:10
%          errorbar(ax(dt),t_data,stable_mean(:,dt),abs(stable_CI(:,dt,1)-stable_mean(:,dt)),abs(stable_CI(:,dt,2)-stable_mean(:,dt)),'o-','Color',col{p})
%        end
      
%        mask = ~isnan(stable_mean);
%        x2 = [t_data(mask) fliplr(t_data(mask))];
%        inBetween = [stable_CI(mask,1)', flipud(stable_CI(mask,2))'];
%        fill(ax1,x2,inBetween,col_fill{p},'FaceAlpha',0.5,'EdgeColor','None','HandleVisibility','off');
%        plot(ax1,t_data,stable_mean,'o','Color',col{p},'DisplayName',lbl_arr{p})
      
%        plot(ax2,t_data,N_norm(:,p),'o-','Color',col{p})
    end
    
%      xlim(ax(1),[0,120])
%      ylim(ax(1),[0,0.5])
    for ds = 1:5
      if ds < 5
        set(ax(ds),'XTick',[])
      end
      set(ax(ds),'FontSize',10)
      text(ax(ds),55,0.45,['\Delta',sprintf('s = %d',ds)],'FontSize',14)
      xlim(ax(ds),[0,120])
%        xlim(ax(ds),[0,15])
      ylim(ax(ds),[0,0.5])
    end
    
%      set(ax1,'FontSize',12)
%      xlim(ax1,[0,t_ses(end)])
%      ylim(ax1,[0,0.55])
    xlabel(ax(5),'\Delta t [h]','FontSize',14)
%      ylabel(ax1,'r^*_{stable}','FontSize',14)
%      legend(ax1,'Location','NorthEast')
    
%      xlim(ax2,[0,t_ses(end)])
%      ylim(ax2,[0,3000])
    
    ax_label = axes('position',[0.12,0.1,0.8,0.8],'visible','off');
    ylabel(ax_label,'r_{stable}','visible','on')
%      set(ax_label,  
    if sv
      path = pathcat(pathFigures,sprintf('time_warp_stability%s.%s',sv_suffix,sv_ext));
      print(path,sprintf('-d%s',sv_ext),'-r300')
      disp(sprintf('Figure saved as %s',path))
    end
    
  end
  
  
  if plot_fig(34)
    
    x_arr = linspace(-para.nbin/2,para.nbin/2,para.nbin+1);
%      figure
    ICI_data = [arrays.ICI_s{:}];
    shift_data = [arrays.shift{:}];
    
    idx_noncoding = [arrays.r_silent{:}] == 0;
    idx_silent = [arrays.r_silent{:}] == 1;
    
    ICI_noncoding = ICI_data(idx_noncoding);
    shift_noncoding = shift_data(idx_noncoding);
    
    ICI_silent = ICI_data(idx_silent);
    shift_silent = shift_data(idx_silent);
    
%      set_norm = 'probability';
    set_norm = 'count';
    
    plt_stuff = true;
    plt_ses = true;
    if plt_ses
%        nT = min(nT,12);
      nT = 10;
      D_KS = zeros(1,nT);
      for ds = 2:2
        %% noncoding
        idx_ds = ICI_noncoding==t_data(ds);
%          if sum(idx_ds) > 20
        hist_noncoding = histcounts(shift_noncoding(idx_ds),x_arr,'Normalization',set_norm);
        [cdf_noncoding,x_noncoding] = ecdf(round(abs(shift_noncoding(idx_ds))));
%          x_noncoding = 0:para.nbin/2;
%          cdf_noncoding = zeros(1,para.nbin/2+1)*NaN;
%          for j = x_noncoding
%            if ismember(j,x_noncoding_tmp)
%              cdf_noncoding(j) = x_noncoding(find(x_noncoding_tmp(2:end)==j)+1);
%            else
%              cdf_noncoding(j) = x_noncoding(find(x_noncoding_tmp(2:end)==j-1)+1);
%            end
%          end
%  %            [cdf_tmp,x] = ecdf(round(abs(shift_noncoding(idx_ds))));
%  %            cdf_noncoding(ds,x(2:end)+1) = cdf_tmp(2:end);
%  %          end
        
        %% silent
        idx_ds = ICI_silent==t_data(ds);
        if sum(idx_ds) > 20
          hist_silent = histcounts(shift_silent(idx_ds),x_arr,'Normalization',set_norm);
          [cdf_silent,x_silent] = ecdf(round(abs(shift_silent(idx_ds))));
        end
%          x_silent = 0:para.nbin/2;
%          cdf_silent = zeros(1,para.nbin/2+1)*NaN;
%          for j = x_silent
%            if ismember(j,x_silent_tmp)
%              cdf_silent(j) = x_silent(find(x_silent_tmp(2:end)==j)+1);
%            else
%              cdf_silent(j) = x_silent(find(x_silent_tmp(2:end)==j-1)+1);
%            end
%          end
%          
%  %            [cdf_tmp,x] = ecdf(round(abs(shift_cont(idx_ds))));
%  %            cdf_cont(ds,x(2:end)+1) = cdf_tmp(2:end);
%  %          end
        
%          [~,max_pos] = max(abs(cdf_silent-cdf_noncoding));
%          D_KS(ds) = cdf_silent(max_pos)-cdf_noncoding(max_pos);
        
        if ds == 2 || plt_stuff
          figure('position',[500 500 400 250])
  %          subplot(1,2,1)
          hold on
          bar(x_arr(1:end-1)+0.5,hist_noncoding,1,'FaceColor','b')
          bar(x_arr(1:end-1)+0.5,-hist_silent,1,'FaceColor','r')
          title(sprintf('ds = %d',ds))
          xlabel('\Delta x')
          ylabel('#(\Delta x)')
          
          ax_cdf = axes('position',[0.7,0.6,0.25,0.3]);
  %          ax_cdf = subplot(1,2,2);
  %          hold(ax_cdf,'on')
  %          figure
  %          ax_cdf = subplot(1,1,1);
          hold on
          stairs(x_noncoding,cdf_noncoding,'Color','b')
          stairs(x_silent,cdf_silent,'Color','r')
          xlabel('|\Delta x|')
          ylabel('cdf(|\Delta x|)')
        end
%          hold(ax_cdf,'off')
      end
      
%        ax_D_KS = axes('position',[0.7,0.2,0.3,0.3])
%        hold on
%        plot([0,nT],[0,0],'k:')
%        plot(1:nT,D_KS,'k')
    else
      
      hist_noncoding = histcounts(shift_noncoding(:),x_arr,'Normalization',set_norm);
      hist_silent = histcounts(shift_silent(:),x_arr,'Normalization',set_norm);
      subplot(1,2,1)
      hold on
      bar(x_arr(1:end-1)+0.5,hist_noncoding,'FaceColor','b')
      bar(x_arr(1:end-1)+0.5,-hist_silent,'FaceColor','r')
      xlabel('\Delta x')
      ylabel('#(\Delta x)')
      
      [cdf_noncoding,x_noncoding] = ecdf(abs(shift_noncoding(:)));
      [cdf_silent,x_silent] = ecdf(abs(shift_silent(:)));
      subplot(1,2,2)
      hold on
      stairs(x_noncoding,cdf_noncoding,'Color','b')
      stairs(x_silent,cdf_silent,'Color','r')
      xlabel('|\Delta x|')
      ylabel('cdf(|\Delta x|)')
      
    end
    
    if sv
      path = pathcat(pathFigures,sprintf('silent_vs_noncoding%s.%s',sv_suffix,sv_ext));
      print(path,sprintf('-d%s',sv_ext),'-r300')
      disp(sprintf('Figure saved as %s',path))
    end
    
  end
  
  
end



%  function dist_poisson(lambda,k)
  
  
  
%  end
  

function set_margins(ax,x_margin,y_margin)
  pos = get(ax, 'Position');
  pos(1) = pos(1)-x_margin;
  pos(3) = pos(3)+2*x_margin;
  pos(2) = pos(2)-y_margin;
  pos(4) = pos(4)+2*y_margin;
  set(ax, 'Position',pos)
end


function [p_out] = nchoosek_sum(n,k,p)
  
  p_out = 0;
  for i=k:n
    p_out = p_out + nchoosek(n,i)*p^i*(1-p)^i;
  end
end


function [p_tmp] = expfit_jackknife(X,Y,W)
  %% jackknifing an exponential fit
  if nargin < 3
    W = ones(size(Y));
  end
  mask = Y>0;
  Y = log(Y);
  
  N_data = sum(mask);
  p_tmp = zeros(N_data,2);
  j=1;
  for i=find(mask)'
    mask_ampl = mask;
    mask_ampl(i) = false;
    
    p_tmp(j,:) = lscov(X(mask_ampl,:),Y(mask_ampl),W(mask_ampl));
    j=j+1;
  end
  
end


function data_out = hsm(data)
  %%% adapted from python version of caiman
  %%% Robust estimator of the mode of a data set using the half-sample mode.
  %%% versionadded: 1.0.3
    
  %%% Create the function that we can use for the half-sample mode
  %%% needs input of sorted data
  
  Ndat = length(data);
  if Ndat == 1
      data_out = data(1);
  elseif Ndat == 2
      data_out = mean(data);
  elseif Ndat == 3
      i1 = data(2) - data(1);
      i2 = data(3) - data(2);
      if i1 < i2
          data_out = mean(data(1:2));
      elseif i2 > i1
          data_out = mean(data(2:end));
      else
          data_out = data(2);
      end
  else
      
      wMin = inf;
      N = floor(Ndat/2) + mod(Ndat,2);
      for i = 1:N
          w = data(i+N-1) - data(i);
          if w < wMin
              wMin = w;
              j = i;
          end
      end
      data_out = hsm(data(j:j+N-1));
  end
end


function [spikeNr,md,sd_r] = get_spikeNr(data,time)
  md = hsm(sort(data));       % Find the mode
  
  % only consider values under the mode to determine the noise standard deviation
  ff1 = data - md;
  ff1 = -ff1 .* (ff1 < 0);
  
  % compute 25 percentile
  ff1 = sort(ff1);
  ff1(ff1==0) = NaN;
  Ns = round(sum(ff1>0) * .5);
  
  % approximate standard deviation as iqr/1.349
  iqr_h = ff1(end-Ns);
  sd_r = 2 * iqr_h / 1.349;
  data_thr = md+2*sd_r;
  spikeNr = sum(floor(data/data_thr));
  
end



function plot_CI_as_fill(mean,CI,x_arr,ax,color,plt_label)
  
  hold(ax,'on')
  %% make sure, arrays are properly sized
  if size(mean,1) > size(mean,2)
    mean = mean';
  end
  if size(CI,1) > size(CI,2)
    CI = CI';
  end
  if size(x_arr,1) > 1
    x_arr = x_arr';
  end
  
  if size(CI,1) == 1  %% CI provided as symmetric value (e.g. variance)
    CI = [mean-CI;mean+CI];
  end
  mask = ~isnan(mean);
  x = [x_arr(mask) fliplr(x_arr(mask))];
  inBetween = [CI(1,mask), fliplr(CI(2,mask))];
  
  fill(ax,x,inBetween,color{2},'FaceAlpha',0.5,'EdgeColor','None','HandleVisibility','off');
  if isempty(plt_label)
    plot(ax,x_arr(mask),mean(mask),'-','Color',color{1},'LineWidth',2,'HandleVisibility','off')
  else
    plot(ax,x_arr(mask),mean(mask),'-','Color',color{1},'LineWidth',2,'DisplayName',plt_label)
  end
end


function [bs_median, bs_CI] = bootstrapping(data,N_bs,mode)
  
  %%% data should have entries:
  %% 1st dim: data to be bootstrapped
  %% 2nd dim: # independent points to be bootstrapped
  
  N_samples = size(data,1);
  N_data = size(data,2);
  bs_stats = zeros(N_data,N_bs);  %% mean
  for ds = 1:N_data
    N_s = N_samples-ds;
    bs_samples = randi(N_s,N_s,N_bs);
    for l=1:N_bs
      dat_bs = data(bs_samples(:,l),ds);
      if mode=='mean'
        bs_stats(ds,l) = nanmean(dat_bs);
      end
    end
  end
  
  bs_median = median(bs_stats,2,'omitnan');
  bs_CI = prctile(bs_stats,[2.5,97.5],2);
end


function save_fig(pathFigures,mouse,fig_name,sv_suffix,sv_ext,fig_pos)
  gcf
  set(gcf,'OuterPosition',fig_pos,'InnerPosition',fig_pos)
  path = pathcat(pathFigures,sprintf('m%s_%s_%s.%s',mouse,fig_name,sv_suffix,sv_ext));
  print(path,sprintf('-d%s',sv_ext),'-r300')
  disp(sprintf('Figure saved as %s',path))
end