
function pathData = cell_registration(data,xdata,model,para,p_thr,nSes,mouse,basePath)

  %%% -------------------------------------- Cell registration ---------------------------------------
    
    %% implements also check for one-sided spatial correlation
    %% include obtaining an overall spatial filter and an overall centroid (from almost surely matched cells), and match other cells according to that (or get maximum probability, given single ROIs and overall ROI)
    %% include, that cell shape and position can change over time, therefore give higher significance to cells that can be matched from neighbouring sessions!
    
    disp('preparing data...')
    tic
    for s = 1:nSes
      data(s).registered = false(data(s).nROI,1);
      data(s).matched = false(data(s).nROI,1);
      for sm = 1:nSes;
        xdata(s,sm).p_same_dist = sparse(data(s).nROI,data(sm).nROI);
        xdata(s,sm).p_same_joint = sparse(data(s).nROI,data(sm).nROI);
      end
    end
    toc
    
    %% here, all cells that are not surely different are put together in ROI_clusters
    disp('writing data...')
    tic
    for s = 1:nSes
      for sm = s+1:nSes
        for n = 1:data(s).nROI
          neighbours = xdata(s,sm).neighbours(n,:) > 0;
          idx_nb = find(neighbours);
          idx_dist = max(1,ceil(para.nbins*xdata(s,sm).dist(n,neighbours)/para.dist_max));
          idx_corr = max(1,ceil(para.nbins*xdata(s,sm).fp_corr(n,neighbours)/para.corr_max));
          
          val_tmp = model.p_same_dist(idx_dist);
          xdata(s,sm).p_same_dist(n,neighbours) = val_tmp;
          xdata(sm,s).p_same_dist(neighbours,n) = val_tmp;
          
          for i=1:length(idx_dist)
            val_tmp = model.p_same_joint(idx_dist(i),idx_corr(i));
            xdata(s,sm).p_same_joint(n,idx_nb(i)) = val_tmp;
            xdata(sm,s).p_same_joint(idx_nb(i),n) = val_tmp;
          end
        end
      end
    end
    toc
    
    disp('ROI_clustering...')
    nCluster = 0;
    ROI_cluster = struct;
    tic
    for s = 1:nSes
      for sm = 1:nSes
        if sm == s
          continue
        end
        
        for n = 1:data(s).nROI
          
          if ~data(s).registered(n)   %% add new ROI_cluster if not already belonging to one
            nCluster = nCluster + 1;
            data(s).cluster(n).ID = nCluster;
            data(s).registered(n) = true;
            ROI_cluster(data(s).cluster(n).ID).list = zeros(nSes,1);
            ROI_cluster(data(s).cluster(n).ID).list(s) = n;
            ROI_cluster(data(s).cluster(n).ID).ct = 1;
          end
          
          cluster_candidates = find(xdata(s,sm).p_same_dist(n,:)>(1-p_thr));    %% all ROIs in sm that are candidates to be same as ROI (s,n)
          
          for m = cluster_candidates
            if ~data(sm).registered(m)
              data(sm).cluster(m).ID = data(s).cluster(n).ID;
              data(sm).registered(m) = true;
            end
            
            for c = data(s).cluster(n).ID         %% in case, there are several from one session that fall into one cluster
              if ~ismember(m,ROI_cluster(c).list(sm,:))
                idx = nnz(ROI_cluster(c).list(sm,:))+1;
                ROI_cluster(c).list(sm,idx) = m;
                if idx == 1
                  ROI_cluster(c).ct = ROI_cluster(c).ct + 1;
                end
              end
            end
          end
        end
      end
    end
    
    nMatches = [ROI_cluster.ct];
    disp(sprintf('number of ROI_clusters: %d',nCluster))
    disp(sprintf('number of real ROI_clusters: %d',sum(nMatches > 1)))
    toc
    disp('pre-clustering done')
    
    pathData = sprintf('%s%d/cluster_registered.mat',basePath,mouse);
    save(pathData,'ROI_cluster','data','xdata','model','para','-v7.3')
    disp(sprintf('data saved under %s',pathData))
    
    figure
    histogram(nMatches)
    
end

