

function switch_sessions(pathMouse,s1,s2)
    
    nSes = 15;
    
    %%% --- footprints ---
    pathFP = pathcat(pathMouse,'footprints.mat');
    load(pathFP);
    footprints_old = footprints;
    
    footprints.session(s1) = footprints_old.session(s2);
    footprints.session(s2) = footprints_old.session(s1);
    
    footprints.data.session(s1) = footprints_old.data.session(s2);
    footprints.data.session(s2) = footprints_old.data.session(s1);
    
%      pathFP_try = pathcat(pathMouse,'footprints_try.mat');
    save(pathFP,'footprints','-v7.3')
    
    
    %%% --- xdata ---
    pathXData = pathcat(pathMouse,'xdata.mat');
    load(pathXData);
    xdata_old = xdata;
    
    for s = 1:nSes
        if s~=s1 && s~=s2
            xdata(s1,s) = xdata_old(s2,s);
            xdata(s,s1) = xdata_old(s,s2);
            
            xdata(s2,s) = xdata_old(s1,s);
            xdata(s,s2) = xdata_old(s,s1);
        end
    end
    
    xdata(s1,s1) = xdata_old(s2,s2);
    xdata(s2,s2) = xdata_old(s1,s1);
    
    xdata(s1,s2) = xdata_old(s2,s1);
    xdata(s2,s1) = xdata_old(s1,s2);
    
%      pathXData_try = pathcat(pathMouse,'xdata_try.mat');
    save(pathXData,'xdata','-v7.3')
    
    
    %%% --- clusters ---
    pathC = pathcat(pathMouse,'matching_results.mat');
    load(pathC);
    
    % data
    data_old = data;
    data.session(s1) = data_old.session(s2);
    data.session(s2) = data_old.session(s1);
    
    % status
    status_old = status;
    status.session(s1) = status_old.session(s2);
    status.session(s2) = status_old.session(s1);
    
    for i = 1:length(status.session(s1).manipulation)
        for j = 1:length(status.session(s1).manipulation(i).pre)
            status.session(s1).manipulation(i).pre(j).ID(1) = s1;
        end
        for j = 1:length(status.session(s1).manipulation(i).post)
            status.session(s1).manipulation(i).post(j).ID(1) = s1;
        end
    end
    for i = 1:length(status.session(s2).manipulation)
        for j = 1:length(status.session(s2).manipulation(i).pre)
            status.session(s2).manipulation(i).pre(j).ID(1) = s2;
        end
        for j = 1:length(status.session(s2).manipulation(i).post)
            status.session(s2).manipulation(i).post(j).ID(1) = s2;
        end
    end
    
    for i = 1:length(status.manipulation)
        if status.manipulation(i).idx(1) == s1
            status.manipulation(i).idx(1) = s2;
        elseif status.manipulation(i).idx(1) == s2;
            status.manipulation(i).idx(1) = s1;
        end
    end
    
    % clusters
    clusters_old = clusters_sv;
    nC = length(clusters_sv);
    
    for c = 1:nC
        clusters_sv(c).session(s1) = clusters_old(c).session(s2);
        clusters_sv(c).session(s2) = clusters_old(c).session(s1);
    end
    
%      pathC_try = pathcat(pathMouse,'matching_results_try.mat');
    save(pathC,'clusters_sv','status','data','-v7.3')
    %% session folders (manually)
    
    
end