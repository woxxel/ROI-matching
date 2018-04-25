

function [merge_status] = check_merge(data,xdata,s_ref,n_ref,s,n)
  
  %% get ROI idxes and throw them all together to see merging result
  merge_status = true;
  microns_per_pixel = 530.684/512;
  
  A_thr = 400;
  
    
  imSize = [512,512];
  
  %% merge ROIs and compute statistics to compare to
  A_merge = sparse(imSize(1)*imSize(2),1);
  fp_corr_oneway_1_ref = zeros(length(n),1);
  fp_corr_oneway_2_ref = zeros(length(n),1);
  
  for i = 1:length(n)
    A_merge = A_merge + data(s).A(:,n(i));
    
    fp_corr_ref(i) = xdata(s,s_ref).fp_corr(n(i),n_ref);
    
    %% compute one-sided correlation
    idx = find(data(s).A(:,n(i)));
    fp_corr_oneway_1_ref(i) = full(dot(data(s).A(idx,n(i)),data(s_ref).A(idx,n_ref))/(data(s).norm(n(i))*norm(data(s_ref).A(idx,n_ref))));
    
    idx = find(data(s_ref).A(:,n_ref));
    fp_corr_oneway_2_ref(i) = full(dot(data(s_ref).A(idx,n_ref),data(s).A(idx,n(i)))/(data(s_ref).norm(n_ref)*norm(data(s).A(idx,n(i)))));
  end
  
  A_area = nnz(A_merge);
  A_ref_area = nnz(data(s_ref).A(:,n_ref));
  
  %% check whether they are connected
  A_test = bwconncomp(full(reshape(A_merge,imSize(1),imSize(2)))>0,8);
  %% check whether 1. ROIs are connected, 2. size of ROI becomes too large, 3. size difference of ROIs is too large
  if A_test.NumObjects > 1
%      disp('dont merge, they are not connected / merged ROI is too large!')
%      A_test.NumObjects
    merge_status = false;
    return
  end
  
  if nnz(A_merge) > A_thr
%      disp('ROI became too large')
%      nnz(A_merge)
    merge_status = false;
    return
  end
  
  if abs(A_area - A_ref_area)/A_ref_area > 0.5 % get area threshold from 95% distribution or so
%      disp('ROI become too much larger than ref')
%      A_area
%      A_ref_area
%      abs(A_area - A_ref_area)/A_ref_area
    merge_status = false;
    return
  end
  
  
  A_norm = norm(A_merge);
  A_tmp = A_merge/sum(A_merge(:));
  A_centroid = [sum((1:imSize(1))*reshape(A_tmp,imSize(1),imSize(2))),sum(reshape(A_tmp,imSize(1),imSize(2))*(1:imSize(2))')];
%                      A_merge = reshape(A_merge,imSize(1)*imSize(2),1);

  dist = microns_per_pixel*sqrt((A_centroid(1) - data(s_ref).centroid(n_ref,1)).^2 + (A_centroid(2) - data(s_ref).centroid(n_ref,2)).^2);
  fp_corr = full(dot(A_merge,data(s_ref).A(:,n_ref))/(A_norm*data(s_ref).norm(n_ref)));
  
  %% check new one-way correlation
  idx = find(A_merge);
  fp_corr_oneway_1 = full(dot(A_merge(idx),data(s_ref).A(idx,n_ref))/(A_norm*norm(data(s_ref).A(idx,n_ref))));
  
  idx = find(data(s_ref).A(:,n_ref));
  fp_corr_oneway_2 = full(dot(data(s_ref).A(idx,n_ref),A_merge(idx))/(data(s_ref).norm(n_ref)*norm(A_merge(idx))));
  
  %                      if 3/4 of 1-way correlations increase && overall correlation increases %% distance stays about same or better, mark as merging candidate and add to this cluster (the final question of merging comes later, when the whole cluster was processed)
  
  if (sum(fp_corr_oneway_1>fp_corr_oneway_1_ref)+sum(fp_corr_oneway_2 > fp_corr_oneway_2_ref)) > (3./2*length(n)) && all(fp_corr > fp_corr_ref)
%      disp('seems alright - compare Ca traces!')
  else
    merge_status = false;
%      disp('nope - dont merge!')
    return
  end
  
  
  plt=false;
  if plt
    %%% prepare figure
    %% display new vs old stats
    string_stat1 = '';
    for i=1:length(n)
      string_stat1 = sprintf('%s \n n=%d \t old distance: %4.2g, old correlation %4.2g',string_stat1,n(i),full(xdata(s,s_ref).dist(n(i),n_ref)),full(xdata(s,s_ref).fp_corr(n(i),n_ref)));
    end
    string_stat1 = sprintf('%s \n new distance: %4.2g, new correlation %4.2g',string_stat1,dist,fp_corr);
    
    string_stat2 = '';
    for i=1:length(n)
      string_stat2 = sprintf('%s \n n=%d \t old 1-way corr(1): %4.2g, old 1-way corr(2): %4.2g',string_stat2,n(i),fp_corr_oneway_1_ref(i),fp_corr_oneway_2_ref(i));
    end
    string_stat2 = sprintf('%s \n new 1-way corr(1): %4.2g, new 1-way corr(2) %4.2g',string_stat2,fp_corr_oneway_1,fp_corr_oneway_2);
    
  %    disp('---------------- merging? ---------------------')
    %% if this gives comparable or better results, check temporal correlation
    
  
    figure
    hold on
    for n_merge = n
      contour(reshape(data(s).A(:,n_merge),imSize(1),imSize(2)),[0.01,0.01]*max(data(s).A(:,n_merge)),'b')
    end
    contour(reshape(A_merge,imSize(1),imSize(2)),[0.01,0.01]*max(A_merge),'r')
    contour(reshape(data(s_ref).A(:,n_ref),imSize(1),imSize(2)),[0.01,0.01]*max(data(s_ref).A(:,n_ref)),'k')
    hold off
    text(A_centroid(2)-18,A_centroid(1)+15,string_stat1,'fontsize',10,'fontweight','bold');
    text(A_centroid(2)-18,A_centroid(1)-15,string_stat2,'fontsize',10,'fontweight','bold');
    
    xlim([A_centroid(2)-20 A_centroid(2)+20])
    ylim([A_centroid(1)-20 A_centroid(1)+20])
  end
  %% ask user if this should merge? (maybe at end of matching, then display all ROIs (in 3D?) and ask if it looks alright
end