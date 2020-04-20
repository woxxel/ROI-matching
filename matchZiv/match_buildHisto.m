

function [histo] = match_buildHisto(paths,xdata,para)
  
  
  if ~exist(paths.histo,'file')
    disp(sprintf('building and saving histo to %s',paths.histo))
    histo = struct;
    
    dbin_dist = para.dist_max/(para.nbins*2);
    histo.dist_x = linspace(dbin_dist,para.dist_max-dbin_dist,para.nbins)';
    
    dbin_corr = 1./(para.nbins*2);
    histo.corr_x = linspace(dbin_corr,1-dbin_corr,para.nbins)';
    
    %%% histograms have 2 entries: 1st for nearest neighbour, 2nd for other neighbours
    histo.dist = zeros(para.nbins,2);
    histo.corr = zeros(para.nbins,2);
    
    if para.bool_shift_corr
      histo.shift_corr = zeros(para.nbins,2);
      histo.shift_1w_corr_min = zeros(para.nbins,2);
      histo.shift_1w_corr_max = zeros(para.nbins,2);
    end
    
    histo.joint = zeros(para.nbins,para.nbins,2);
    
    tic
    for s = 1:para.nSes
      for sm = 1:para.nSes
        if sm-s == 1
  %          if s~=sm
            histo = update_histogram(xdata(s,sm),histo,para);
  %          end
  %        elseif s > sm
  %          histo = update_histogram(xdata(s,sm),histo,para);
        end
      end
    end
    toc
    save(paths.histo,'histo','-v7.3')
  else
    disp(sprintf('loading histo from %s',paths.histo))
    load(paths.histo)
  end
  
end



function [histo] = update_histogram(xdata,histo,para)
    
    nROI = size(xdata.dist,1);
    %%% update distance-, fp_correlation- and joint-histogram
    idx_dist = max(1,ceil(para.nbins*xdata.dist(xdata.neighbours>0)/para.dist_max));
    idx_corr = max(1,ceil(para.nbins*xdata.corr(xdata.neighbours>0)/para.corr_max));
    
    if para.bool_shift_corr
      idx_shift_corr = max(1,ceil(para.nbins*xdata.shift_corr(xdata.neighbours>0)/para.corr_max));
      idx_shift_1w_corr_min = max(1,ceil(para.nbins*xdata.shift_1w_corr_min(xdata.neighbours>0)/para.corr_max));
      idx_shift_1w_corr_max = max(1,ceil(para.nbins*xdata.shift_1w_corr_max(xdata.neighbours>0)/para.corr_max));
    end
    
    for i=1:length(idx_dist)
      histo.dist(idx_dist(i),2) = histo.dist(idx_dist(i),2) + 1;
      histo.corr(idx_corr(i),2) = histo.corr(idx_corr(i),2) + 1;
      histo.joint(idx_dist(i),idx_corr(i),2) = histo.joint(idx_dist(i),idx_corr(i),2) + 1;
      
      if para.bool_shift_corr
        histo.shift_corr(idx_shift_corr(i),2) = histo.shift_corr(idx_shift_corr(i),2) + 1;
        histo.shift_1w_corr_min(idx_shift_1w_corr_min(i),2) = histo.shift_1w_corr_min(idx_shift_1w_corr_min(i),2) + 1;
        histo.shift_1w_corr_max(idx_shift_1w_corr_max(i),2) = histo.shift_1w_corr_max(idx_shift_1w_corr_max(i),2) + 1;
      end
    end
    
    %%% find nearest neighbours
    for n = 1:nROI
      if nnz(xdata.neighbours(n,:))
        min_dist = min(xdata.dist(n,xdata.neighbours(n,:)>0));
        m = find((xdata.dist(n,:)==min_dist) & (xdata.neighbours(n,:)>0));
        
        if xdata.corr(n,m) > 0.0    %% only those with some overlap are considered to be candidates for same cells = nearest neighbours
%            xdata.neighbours(n,m) = 2;
          
          NN_idx_dist = max(1,ceil(para.nbins*min_dist/para.dist_max));
          histo.dist(NN_idx_dist,1) = histo.dist(NN_idx_dist,1) + 1;
          histo.dist(NN_idx_dist,2) = histo.dist(NN_idx_dist,2) - 1;
          
          NN_idx_corr = max(1,ceil(para.nbins*xdata.corr(n,m)/para.corr_max));
          histo.corr(NN_idx_corr,1) = histo.corr(NN_idx_corr,1) + 1;
          histo.corr(NN_idx_corr,2) = histo.corr(NN_idx_corr,2) - 1;
          
          histo.joint(NN_idx_dist,NN_idx_corr,1) = histo.joint(NN_idx_dist,NN_idx_corr,1) + 1;
          histo.joint(NN_idx_dist,NN_idx_corr,2) = histo.joint(NN_idx_dist,NN_idx_corr,2) - 1;
          
          if para.bool_shift_corr
            %% additional stuff
            NN_idx_corr = max(1,ceil(para.nbins*xdata.shift_corr(n,m)/para.corr_max));
            histo.shift_corr(NN_idx_corr,1) = histo.shift_corr(NN_idx_corr,1) + 1;
            histo.shift_corr(NN_idx_corr,2) = histo.shift_corr(NN_idx_corr,2) - 1;
            
            NN_idx_corr = max(1,ceil(para.nbins*xdata.shift_1w_corr_min(n,m)/para.corr_max));
            histo.shift_1w_corr_min(NN_idx_corr,1) = histo.shift_1w_corr_min(NN_idx_corr,1) + 1;
            histo.shift_1w_corr_min(NN_idx_corr,2) = histo.shift_1w_corr_min(NN_idx_corr,2) - 1;
            
            NN_idx_corr = max(1,ceil(para.nbins*xdata.shift_1w_corr_max(n,m)/para.corr_max));
            histo.shift_1w_corr_max(NN_idx_corr,1) = histo.shift_1w_corr_max(NN_idx_corr,1) + 1;
            histo.shift_1w_corr_max(NN_idx_corr,2) = histo.shift_1w_corr_max(NN_idx_corr,2) - 1;
          end
        end
      end
    end
    
    %%% calculate mean here as well?!
end