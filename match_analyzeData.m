
function [xdata, histo, para] = match_analyzeData(data,nSes,dist_max)
    
    para.microns_per_pixel = 530.684/512;
    para.corr_max = 1;
    para.dist_max = dist_max;
    para.nbins = 100;
    
    dbin_dist = para.dist_max/(para.nbins*2);
    histo = struct;
    histo.dist_x = linspace(dbin_dist,para.dist_max-dbin_dist,para.nbins)';
    
    dbin_corr = 1./(para.nbins*2);
    histo.corr_x = linspace(dbin_corr,1-dbin_corr,para.nbins)';
    
    tic
    %%% histograms have 2 entries: 1st for nearest neighbour, 2nd for other neighbours
    histo.dist = zeros(para.nbins,2);
    histo.corr = zeros(para.nbins,2);
    histo.joint = zeros(para.nbins,para.nbins,2);
    
    xdata = struct;
    %% now, process them
    for s = 1:nSes
      
      for sm = 1:nSes
        if sm > s
  %            xdata(s,sm).neighbour = sparse(data(s).nROI,data(sm).nROI);
          xdata(s,sm).fp_corr = sparse(data(s).nROI,data(sm).nROI);
          xdata(s,sm).dist = sparse(data(s).nROI,data(sm).nROI);
          xdata(s,sm).p_same = sparse(data(s).nROI,data(sm).nROI);
          
          dist_tmp = zeros(data(s).nROI,data(sm).nROI);
          for n = 1:data(s).nROI
            dist_tmp(n,:) = para.microns_per_pixel*sqrt((data(s).centroid(n,1) - data(sm).centroid(:,1)).^2 + (data(s).centroid(n,2) - data(sm).centroid(:,2)).^2);
          end
          xdata(s,sm).neighbours = 1*(dist_tmp < para.dist_max);
          xdata(s,sm).dist(xdata(s,sm).neighbours>0) = dist_tmp(xdata(s,sm).neighbours>0);
          
          for n = 1:data(s).nROI
            for m = find(xdata(s,sm).neighbours(n,:)>0)
              xdata(s,sm).fp_corr(n,m) = dot(data(s).A(:,n),data(sm).A(:,m))/(data(s).norm(n)*data(sm).norm(m));
            end
          end
          
          [histo,xdata_tmp] = update_histogram(xdata(s,sm),histo,para);
          xdata(s,sm) = xdata_tmp;
          
        elseif s > sm
          %% mostly copying, just nearest neighbours have to be recalculated
          xdata(s,sm).dist = xdata(sm,s).dist';
          xdata(s,sm).fp_corr = xdata(sm,s).fp_corr';
          xdata(s,sm).neighbours = 1*(xdata(sm,s).neighbours>0)';  %% possibly different nearest neighbours
          
          [histo,xdata_tmp] = update_histogram(xdata(s,sm),histo,para);
          xdata(s,sm) = xdata_tmp;
        end
      end
    end    
    toc
    
end


function [histo,xdata] = update_histogram(xdata,histo,para)
    
    nROI = size(xdata.dist,1);
    %%% update distance-, fp_correlation- and joint-histogram
    idx_dist = max(1,ceil(para.nbins*xdata.dist(xdata.neighbours>0)/para.dist_max));
    idx_corr = max(1,ceil(para.nbins*xdata.fp_corr(xdata.neighbours>0)/para.corr_max));
    
    for i=1:length(idx_dist)
      histo.dist(idx_dist(i),2) = histo.dist(idx_dist(i),2) + 1;
      histo.corr(idx_corr(i),2) = histo.corr(idx_corr(i),2) + 1;
      histo.joint(idx_dist(i),idx_corr(i),2) = histo.joint(idx_dist(i),idx_corr(i),2) + 1;
    end
    
    %%% find nearest neighbours
    for n = 1:nROI
      if nnz(xdata.neighbours(n,:))
        min_dist = min(xdata.dist(n,xdata.neighbours(n,:)>0));
        m = find((xdata.dist(n,:)==min_dist) & (xdata.neighbours(n,:)>0));
        
        if xdata.fp_corr(n,m) > 0.0    %% only those with some overlap are considered to be candidates for same cells = nearest neighbours
          xdata.neighbours(n,m) = 2;
          
          NN_idx_dist = ceil(para.nbins*min_dist/para.dist_max);
          histo.dist(NN_idx_dist,1) = histo.dist(NN_idx_dist,1) + 1;
          histo.dist(NN_idx_dist,2) = histo.dist(NN_idx_dist,2) - 1;
          
          NN_idx_corr = max(1,ceil(para.nbins*xdata.fp_corr(n,m)/para.corr_max));
          histo.corr(NN_idx_corr,1) = histo.corr(NN_idx_corr,1) + 1;
          histo.corr(NN_idx_corr,2) = histo.corr(NN_idx_corr,2) - 1;
          
          histo.joint(NN_idx_dist,NN_idx_corr,1) = histo.joint(NN_idx_dist,NN_idx_corr,1) + 1;
          histo.joint(NN_idx_dist,NN_idx_corr,2) = histo.joint(NN_idx_dist,NN_idx_corr,2) - 1;
        end
      end
    end
    
    %%% calculate mean here as well?!
end