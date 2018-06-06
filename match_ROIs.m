
function [ROI_cluster,ROI_data] = match_ROIs(nSes,basePath,mouse)
%  function [histo,model,para] = match_ROIs(nSes,basePath,mouse)

  p_thr = 0.8;
  max_dist = 12;
  
  [data] = match_loadSessions(nSes,mouse,basePath);
  
  [xdata, histo, para] = match_analyzeData(data,nSes,max_dist);
  
  [model,histo] = match_buildModel(xdata,histo,para,nSes,mouse,basePath);
  
%    estimate_model_accuracy(histo,model,para,mouse,basePath);
  
%  %    savePath = sprintf('%s%d/matching_data.m',basePath,mouse);
%  %    save(savePath,'data','xdata','histo','model','para','-v7.3')
  
  pathData = cell_registration(data,xdata,model,para,p_thr,nSes,mouse,basePath);
  
  [ROI_cluster, ROI_data] = cell_matching(basePath,mouse,0.9);
  
end