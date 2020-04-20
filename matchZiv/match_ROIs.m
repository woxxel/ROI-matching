
function match_ROIs(nSes,basePath,mouse)
%  function [histo,model,para] = match_ROIs(nSes,basePath,mouse)

  p_thr = 0.8;
  max_dist = 12;
  
%    [data] = match_loadSessions(nSes,mouse,basePath);
  
  paths = struct;
  paths.mouse = pathcat(basePath,mouse);
  paths.footprints_all = pathcat(paths.mouse,'footprints.mat');
  paths.xdata = pathcat(paths.mouse,'xdata.mat');
  paths.histo = pathcat(paths.mouse,'histo.mat');
  paths.results = pathcat(paths.mouse,'matching_results.mat');
  
  [footprints] = match_loadSessions(paths,nSes,[]);
  [xdata, para] = match_buildXData(paths,footprints,footprints.data.nSes,12);      %% calculate distances and footprint correlation
  [histo] = match_buildHisto(paths,xdata,para);
  [model,histo] = match_buildModel(histo,para,footprints.data.nSes,paths.mouse);
  
  [xdata] = match_assign_prob(xdata,footprints.data,model,para,paths.xdata);
  setappdata(0,'xdata',xdata)
  
  match_matching(footprints.data,0.5);
  
%    estimate_model_accuracy(histo,model,para,mouse,basePath);
  
%  %    savePath = sprintf('%s%d/matching_data.m',basePath,mouse);
%  %    save(savePath,'data','xdata','histo','model','para','-v7.3')
  
%    pathData = cell_registration(data,xdata,model,para,p_thr,nSes,mouse,basePath);
  
%    [ROI_cluster, ROI_data] = cell_matching(basePath,mouse,0.9);
  
end