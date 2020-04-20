%%% requires: 
%%%   data:           struct to detect ROI positions
%%%   ROI_cluster:    information about current clusters
%%%   threshold:      both in (1) ROI number and (2) ROI score, to identify ROIs which should be filled up
%%%   (method):       'extrapolation' -           obtain ROIs as weighted ROI from neighbouring sessions
%%%                   'CNMF_all'      -           rerun CNMF method on all sessions (should be quite lengthy...)
%%%                   'CNMF_patches'  - (default) rerun CNMF method in small patches 

%%% the input cluster to this function can still have several ROIs per session, however only if they were marked for merging

%%% should I rather gather all ROIs to be inferred, first and run CNMF then? (or is it too few "force" on specific regions?)
%%% how is it timewise? is there much overhead when starting single iterations? (but can be parallelized nicely)
%%% how does initial guess of A look like? binary? doubles?



function gap_filling(data,ROI_cluster,threshold,method,plt)


  margin = 5;
  
  tic 
  if nargin < 4 || isempty(method)
    method = 'CNMF_patches';
  end
  
  nCluster = length(ROI_cluster);
  
  %%% find all ROIs that should be filled up
  nSes = size(ROI_cluster(1).list,1);
  scores = [ROI_cluster.score];
  counts = [ROI_cluster.ct];
  cluster_status = counts > threshold(1) & scores > threshold(2);
  good_clusters = find(cluster_status);
  fill_ct = 0;
  
  data(1).imSize = [512,512];
  for s=1:nSes
    data(s).A_fill = zeros(data(1).imSize(1),data(1).imSize(2));
    data(s).fill_ct = 0;
    data(s).fill_idx = false(nCluster,1);
  end
    
  extents_prototype = ones(2,2);
  extents_prototype(:,1) = data(1).imSize;
  
  for c = good_clusters
%      disp(sprintf('filling up cluster %04d, \t count: %02d, \t score: %4.2g',c,counts(c),scores(c)))
    ROI_cluster(c).extents = extents_prototype;
    ROI_cluster(c).A_cluster = sparse(data(1).imSize(1),data(1).imSize(2));
    
    %% identify region of interest (x and y extents +/- some margin) and get merged ROI
    [s_idx w_idx] = find(ROI_cluster(c).list);
    entries = length(s_idx);
    for i = 1:entries;
      s = s_idx(i);
      n = ROI_cluster(c).list(s,w_idx(i));
      
      A_tmp = reshape(data(s).A(:,n),data(1).imSize(1),data(1).imSize(2));
      [y_idx x_idx] = find(A_tmp);
      
      if strcmp(method,'CNMF_patches')
        %% only needed if CNMF is run on patches
        ROI_cluster(c).extents(1,1) = max(1,min(ROI_cluster(c).extents(1,1),min(y_idx)-margin));
        ROI_cluster(c).extents(1,2) = min(data(s).imSize(1),max(ROI_cluster(c).extents(1,2),max(y_idx)+margin));
        ROI_cluster(c).extents(2,1) = max(1,min(ROI_cluster(c).extents(2,1),min(x_idx)-margin));
        ROI_cluster(c).extents(2,2) = min(data(s).imSize(2),max(ROI_cluster(c).extents(2,2),max(x_idx)+margin));
      end
      ROI_cluster(c).A_cluster = ROI_cluster(c).A_cluster + A_tmp;
      
    end
    
    ROI_cluster(c).A_cluster = ROI_cluster(c).A_cluster/sum(ROI_cluster(c).A_cluster(:));   % normalize
    ROI_cluster(c).centroid = [sum((1:data(1).imSize(1))*ROI_cluster(c).A_cluster),sum(ROI_cluster(c).A_cluster*(1:data(1).imSize(2))')];
    
    %%% shouldn't I rather first go through all clusters to identify locations of rerunning CNMF, then identifying small patches, uniting those and rerun CNMF on those. (so a nearby, still missing ROI will not interfere with the method, but can instead be identified in the same run)
    %%% gather all, then treat new and old ROIs same when intializing the CNMF
    
    for s=1:nSes
      if ~nnz(ROI_cluster(c).list(s,:))
        data(s).fill_idx(c) = true;
        if plt
          data(s).fill_ct = data(s).fill_ct + 1;
          data(s).A_fill = data(s).A_fill + ROI_cluster(c).A_cluster;
        end
      end
    end
  
  end
  toc
  
  
  
  %%% display all ROIs inferred ROI positions on which to rerun CNMF method
  %%% construct overall image of first guesses
  %%% remove sub-optimal ROIs here?  
  
  n_good = data(1).cluster_neuron(cluster_status);      %% neurons from this session to be retained
  n_good = n_good(n_good>0);
  n_good_sz = length(n_good);
  
  n_bad = data(1).cluster_neuron(~cluster_status);          %% neurons from this session to be discarded
  n_bad = n_bad(n_bad>0);
  
  c_add = find(data(1).fill_idx)';
  c_add_sz = length(c_add);
  
  disp('constructing new footprint guess')
  
  newData = struct;
  
  newData.A = zeros(prod(size(ROI_cluster(c).A_cluster)),n_good_sz+c_add_sz);
  newData.centroid = zeros(n_good_sz+c_add_sz,2); 
  newData.Astatus = false(n_good_sz+c_add_sz);
  
  n_ct = 0;
  
%  %    load('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/Session01/resultsCNMF_MF1_LK1.mat','C2','S2')
  
  for n=n_good
    n_ct = n_ct + 1;
    newData.A(:,n_ct) = data(1).A(:,n)./sum(data(1).A(:,n));
    newData.centroid(n_ct,:) = data(1).centroid(n,:);
%      newData.C2(n_ct,:) = C2(n,:);
%      newData.S2(n_ct,:) = S2(n,:);
    newData.status(n_ct) = 0;
  end
  
%    for n=n_bad
%      n_ct = n_ct + 1;
%      newData.A(:,n_ct) = data(1).A(:,n)./sum(data(1).A(:,n));
%      newData.centroid(n_ct,:) = data(1).centroid(n,:);
%      newData.C2(n_ct,:) = C2(n,:);
%      newData.S2(n_ct,:) = S2(n,:);
%      newData.status(n_ct) = 0;
%    end
  
  
  for c=c_add
    n_ct = n_ct + 1;
    newData.A(:,n_ct) = sparse(ROI_cluster(c).A_cluster(:));
    newData.centroid(n_ct,:) = ROI_cluster(c).centroid;
%      newData.Astatus(n_ct) = true;
%      newData.C2(n_ct,:) = NaN;
%      newData.S2(n_ct,:) = NaN;
    newData.status(n_ct) = 1;
  end
  
%    pathSave = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/Session01/results_filling.mat';
%    save(pathSave,'-struct','newData','-v7.3')
%    disp(sprintf('ROIs saved @ path %s',pathSave))
  
  
  
  
  
%    newData.A(:,c_add_sz+1:end) = data(1).A(:,n_good)./sum(data(1).A(:,n_good));  % get (normalized) ROI-footprints that will be retained
%    newData.centroid(c_add_sz+1:end,:) = data(1).centroid(n_good,:);
  
%    newData.Astatus = true(c_add_sz);
%    newData.Astatus(c_add_sz+1:end) = false;  
%      %% calculate correlation of added ROIs with sorted out ones - no need, I'm trying to get better ones anyway (and doesn't seem super correlated anyway)


  if ~plt
    disp('now fill!')
    
    for c=c_add
  %      c = c_add(randi(length(c_add)))
  %      c = 889;
      c
  %      figure
  %      imagesc(ROI_cluster(c).A_cluster)
  %      xlim([ROI_cluster(c).extents(2,1),ROI_cluster(c).extents(2,2)])
  %      ylim([ROI_cluster(c).extents(1,1),ROI_cluster(c).extents(1,2)])
      
      %%% here, sort out out-of-extent ROIs already and mark the "new" ones, that are completely enclosed (get idx)
      %%% take care, that indices remain valid (CaImAn just removes non-interesting ROIs)
  %      newData.Astatus
      
      A = struct;
      A.extents = ROI_cluster(c).extents;
      
      A.footprint = zeros(A.extents(1,2)-A.extents(1,1)+1,A.extents(2,2)-A.extents(2,1)+1,1);
      d = size(A.footprint);
      A.total = zeros(A.extents(1,2)-A.extents(1,1)+1,A.extents(2,2)-A.extents(2,1)+1);
      A.ct = 0;
      for n=1:size(newData.A,2)
        A_tmp = reshape(newData.A(:,n),data(1).imSize(1),data(1).imSize(2));
        A_tmp = A_tmp(A.extents(1,1):A.extents(1,2),A.extents(2,1):A.extents(2,2));
        if nnz(A_tmp) > 10
          A.ct = A.ct + 1;
          A.footprint(:,:,A.ct) = A_tmp;
          A.centroid(A.ct,:) = newData.centroid(n,:);
          A.total = A.total + A_tmp;
          
          if nnz(A_tmp) == nnz(newData.A(:,n)) && newData.status(n)==1
            A.status(A.ct) = 2;
          else
            A.status(A.ct) = newData.status(n);
          end
        end
      end
      
%        A.footprint
%        sum(A.footprint)
      [A_out, C_out, status, fitness] = run_CNMF_fill('',A);
      
%        nam = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/Session01/images/ImagingData_MF1_LK1.h5';
%        Y = double(read_file_crop(nam,ROI_cluster(c).extents));
%        disp('reading done')
%        A_tmp = full(ROI_cluster(c).A_cluster(ROI_cluster(c).extents(1,1):ROI_cluster(c).extents(1,2),ROI_cluster(c).extents(2,1):ROI_cluster(c).extents(2,2)));
%        C = zeros(1,size(Y,3));
%        for t=1:size(Y,3)
%          C(t) = sum(sum(A_tmp.*Y(:,:,t)));
%        end
%        size(C)
%        plot(C)
%        
%        options = struct('N_samples_exc',6,'robust_std',0);
%        fitness = compute_event_exceptionality(C,options.N_samples_exc,options.robust_std);
%        fitness
      
%        
%        status
%        A.status
%        fitness
      inside = sum(sum(A.footprint));
  %      inside = inside(status);
      
%        disp('CNMF done')
%        status(find(A.status==2),3)
      if status(find(A.status==2),3)
        
        load('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/884/Session01/reduced_MF1_LK1.mat','max_im')
        
        figure('position',[300 300 1200 900])
        subplot(2,2,1)
        %% plot of input
        hold on
        imagesc(max_im)
        for k=1:A.ct
          A_tmp = zeros(data(1).imSize(1),data(1).imSize(2));
          A_tmp(A.extents(1,1):A.extents(1,2),A.extents(2,1):A.extents(2,2)) = A.footprint(:,:,k);
          
          if A.status(k) == 2
              contour(A_tmp, [0.1 0.1]*max(A_tmp(:)), 'Color','r')
          elseif A.status(k) == 1
            contour(A_tmp, [0.1 0.1]*max(A_tmp(:)), 'Color','y')
          else
            contour(A_tmp, [0.1 0.1]*max(A_tmp(:)), 'Color','g')
          end
        end
        hold off
        xlim([A.extents(2,1),A.extents(2,2)])
        ylim([A.extents(1,1),A.extents(1,2)])
        
        
        subplot(2,2,2)
        %% plot of output
        hold on
        imagesc(max_im)
        for k=1:size(A_out,2)
          A_tmp = zeros(data(1).imSize(1),data(1).imSize(2));
          A_tmp(A.extents(1,1):A.extents(1,2),A.extents(2,1):A.extents(2,2)) = reshape(A_out(:,k),d(1),d(2));
          if status(k,3)
            if A.status(k) == 2
              contour(A_tmp, [0.1 0.1]*max(A_tmp(:)), 'Color','r')
            elseif A.status(k) == 1
              contour(A_tmp, [0.1 0.1]*max(A_tmp(:)), 'Color','y')
            end
          else
            contour(A_tmp, [0.1 0.1]*max(A_tmp(:)), 'Color','g')
          end
        end
        hold off
        xlim([A.extents(2,1),A.extents(2,2)])
        ylim([A.extents(1,1),A.extents(1,2)])
        
        i_ct = 0;
    %      status
    %      find(status)
        for i = 1:A.ct%find(status(:,3))'
          i_ct = i_ct + 1;
          
          subplot(2*A.ct,1,A.ct+i_ct)
%            subplot(2*sum(status(:,3)),1,sum(status(:,3))+i_ct)
    %        size(C_out)
    %        i
          if inside(i)>0.95 && A.status(i)
    %          if size(C_out,1)==1
    %            plot(C_out(:),'k')
    %          else
              plot(C_out(i_ct,:),'k')
    %          end
          else
    %          if size(C_out,1)==1
    %            plot(C_out(:),'b')
    %          else
              plot(C_out(i_ct,:),'b')
    %          end
          end
        end
        
        pause(0.01)
      end
      %%% preprocess:             should be done (yes!): noise estimation, time constants, etc...
      %%% greedy ROI init:        skip
      %%% refinement with HALS:   do not skip
      %%% then, go in CNMF.fit!    
      %%% obtain only the ROIs which are completely represented and get only the ones which had to be newly inferred
    end
  end
    
  
  
  
  
  
  if plt
    
%      disp('plotting all sessions fill up-ROIs')
%      
%      figure('position',[100 100 1600 1200])
%      for s=1:nSes
%        subplot(floor(sqrt(nSes)),ceil(sqrt(nSes)),s)
%        hold on
%        imagesc(data(s).A_fill)
%        for c = good_clusters
%          if ~nnz(ROI_cluster(c).list(s,:))
%            plot(ROI_cluster(c).centroid(2),ROI_cluster(c).centroid(1),'r+')
%          end
%        end
%        hold off
%        set(gca,'YDir','normal')
%        xlim([1,data(1).imSize(2)])
%        ylim([1,data(1).imSize(1)])
%        title(sprintf('session %d, fills: %d',s,data(s).fill_ct))
%      end
%      plotPath = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/session_fills.png';
%      print(plotPath,'-dpng','-r300')
%      disp(sprintf('figure saved under %s',plotPath))
%      
    
    
    
    disp('plotting cluster examples')
    tic
    fig_ROI_examples = figure('position',[100 100 1600 1200]);
    r = 15;
    for i = 1:16
      subplot(4,4,i)
      c = good_clusters(randi(length(good_clusters)));
      [s_idx w_idx] = find(ROI_cluster(c).list);
      entries = length(s_idx);
      
      hold on
      for j = 1:entries;
        s = s_idx(j);
        n = ROI_cluster(c).list(s,w_idx(j));
        
        col = ones(3,1)*4*s/(5.*nSes);
        
        A_tmp = reshape(data(s).A(:,n),512,512);
        
        labelstr = sprintf('s=%02d',s);
        contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color',col,'DisplayName',labelstr)
      end
      hold off
      xlim([ROI_cluster(c).centroid(2)-r ROI_cluster(c).centroid(2)+r])
      ylim([ROI_cluster(c).centroid(1)-r ROI_cluster(c).centroid(1)+r])
      title(sprintf('#%d, count: %d, score: %4.2g',c,ROI_cluster(c).ct,ROI_cluster(c).score))
    end
    plotPath = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/884/ROI_examples.png';
    print(plotPath,'-dpng','-r300')
    disp(sprintf('figure saved under %s',plotPath))
    toc
    
    
%      disp('plotting whole session')
%      load('/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/Session01/reduced_MF1_LK1.mat','max_im')
%      fig = figure('position',[100 100 1600 1200]);
%      hold on
%      imagesc(max_im)
%      colormap(gca,'gray')
%      set(gca,'clim',[0 max(max_im(:))])
%      
%      plot_blobs(gca,reshape(full(data(1).A(:,n_good)),512,512,[]),[],0.25,'-','g');
%      
%  %      for n=n_good
%  %        A_tmp = reshape(data(1).A(:,n),data(1).imSize(1),data(1).imSize(2));
%  %        contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color','g')
%  %      end
%      disp('good plotting done')
%      
%  %      for n=n_bad
%  %        A_tmp = reshape(data(1).A(:,n),data(1).imSize(1),data(1).imSize(2));
%  %        contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color','r')
%  %      end
%      plot_blobs(gca,reshape(full(data(1).A(:,n_bad)),512,512,[]),[],0.25,':','r');
%      disp('bad plotting done')
%  %      
%  %      for c=c_add
%  %        A_tmp = ROI_cluster(c).A_cluster;
%  %    %      A_tmp = reshape(data.A(:,n),data(1).imSize(1),data(1).imSize(2));
%  %        contour(A_tmp,[0.1 0.1]*max(A_tmp(:)),'Color','y')
%  %      end
%  %      datA = cat(3,full(ROI_cluster(c_add).A_cluster))
%  %      plot_blobs(gca,reshape(full(datA),512,512,[]),[],0.25,'--','y');
%      disp('fill plotting done')
%      hold off
%      xlim([1,data(1).imSize(2)])
%      ylim([1,data(1).imSize(1)])
%      
%      title(sprintf('Session 1, retained neurons: %d, discarded neurons: %d, to be inferred neurons: %d, (really found neurons: %d)',length(n_good),length(n_bad),length(c_add),0),'fontSize',16)
%      plotPath = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/fill_gaps_s1.png';
%      print(plotPath,'-dpng','-r300')
%      disp(sprintf('figure saved under %s',plotPath))
  end
end