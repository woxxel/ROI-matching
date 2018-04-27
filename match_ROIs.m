
function [data,xdata,ROI_cluster,model,para] = match_ROIs(nSes,mouse)
  
  basePath = '/media/wollex/AS2/';
  
  p_thr = 0.95;
  max_dist = 12;
  
  [data] = match_loadSessions(nSes,mouse,basePath);
  
  [xdata, histo, para] = match_analyzeData(data,nSes,max_dist);
  
  [model,histo] = match_buildModel(xdata,histo,para,nSes,mouse,basePath);
  
%    savePath = sprintf('%s%d/matching_data.m',basePath,mouse);
%    save(savePath,'data','xdata','histo','model','para','-v7.3')
  
  [data,xdata,ROI_cluster] = cell_registration(data,xdata,model,para,p_thr,nSes,mouse,basePath);
  
end