

function [data] = match_loadSessions(nSes,mouse)

    rot_max = 1;
    rot = linspace(-rot_max,rot_max,10*rot_max+1);
    
    data(nSes) = struct('shift',[],'rotation',[],'nROI',[],'A',[],'centroid',[],'norm',[]);
    
    basePath = '/media/mizuta/AS2/';
%      basePath = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data';
    
    tic
    
    path_bg = sprintf('%s%d/Session01/reduced_MF1_LK1.mat',basePath,mouse);
    loadDat = load(path_bg,'max_im');
    bg_ref = loadDat.max_im;
    imSize = size(bg_ref);
    
    gcp;
    %% register all sessions
    parfor s = 1:nSes
      disp(sprintf('loading session %02d',s))
      
%        path_ROI = sprintf('/media/wollex/Analyze_AS1/245/Session%02d/resultsCNMF_MF1_LK1.mat',s);
      path_ROI = sprintf('%s%d/Session%02d/resultsCNMF_MF1_LK1.mat',basePath,mouse,s);
      loadDat = load(path_ROI,'A2');
      data(s).nROI = size(loadDat.A2,2);
      
      A_tmp = reshape(full(loadDat.A2),imSize(1),imSize(2),data(s).nROI);
      
      if s == 1
        
        data(s).shift = [0,0,0];
        data(s).rotation = 0;
        
      else
        
%          path_bg = sprintf('/media/wollex/Analyze_AS1/245/Session%02d/imageStack/reduced_MF1_LK1.mat',s);
        path_bg = sprintf('%s%d/Session%02d/reduced_MF1_LK1.mat',basePath,mouse,s);
        loadDat = load(path_bg,'max_im');
        bg_tmp = loadDat.max_im;
        
        max_C = 0;
        rot_tmp = -rot_max;
        for r = rot
          bg_rot = imrotate(bg_tmp,r,'crop');
          C = fftshift(real(ifft2(fft2(bg_ref).*fft2(rot90(bg_rot,2)))));
          if max(C(:)) > max_C;
            max_C = max(C(:));
            rot_tmp = r;
            [ind_y,ind_x] = find(C == max_C);
          elseif max(C(:)) == max_C;
            rot_tmp = [rot_tmp r];
            disp('same')
          end
          [ind_y_tmp,ind_x_tmp] = find(C == max(C(:)));
%            disp(sprintf('max: %5.3g, rot: %5.3g, x/y: %d/%d',max(C(:)),r,(floor(imSize(1)/2) - ind_y_tmp),(floor(imSize(2)/2) - ind_x_tmp)))
        end
        %% need to adjust for possibly different shifts as well!
        disp(sprintf('final shift: rot: %5.3g, x/y: %d/%d',rot_tmp,(floor(imSize(1)/2) - ind_y),(floor(imSize(2)/2) - ind_x)))
        rot_tmp = mean(rot_tmp);
        
        %% imtranslate takes [x,y,z] vector
        data(s).shift(3) = 0;
        data(s).shift(2) = floor(imSize(1)/2) - ind_y;
        data(s).shift(1) = floor(imSize(2)/2) - ind_x;
        data(s).rotation = rot_tmp;
        
        A_tmp = imtranslate(A_tmp,-data(s).shift(:));
        if data(s).rotation ~= 0
          A_tmp = imrotate(A_tmp,data(s).rotation,'crop');
        end
      end
      
      data(s).A = sparse(reshape(A_tmp,imSize(1)*imSize(2),data(s).nROI));
      data(s).centroid = zeros(data(s).nROI,2);
      data(s).norm = zeros(data(s).nROI,1);
      
      for n = 1:data(s).nROI
        %% calculate centroid position
        A_tmp_norm = sparse(A_tmp(:,:,n)/sum(data(s).A(:,n)));
        
        data(s).centroid(n,:) = [sum((1:imSize(1))*A_tmp_norm),sum(A_tmp_norm*(1:imSize(2))')];
        data(s).norm(n) = norm(data(s).A(:,n));
      end
      
      %% normalize
%        data(s).A = data(s).A ./ sum(data(s).A);
      
%        for n = 1:data(s).nROI
%          data(s).norm(n) = norm(data(s).A(:,n));
%        end
    end
    disp('loading done')
    toc
    
end
