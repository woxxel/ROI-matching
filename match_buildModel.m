
function [model,histo] = match_buildModel(xdata,histo,para,nSes)
    
    tic
    
    %% predefine functions to obtain parameters and fits
    logn_mu = @(m,s) log((m^2)/sqrt(s^2+m^2));
    logn_sigma = @(m,s) sqrt(log(s^2/(m^2)+1));
    
    beta_a = @(m,s) m^2*(1-m)/s^2 - m;
    beta_b = @(m,alpha) alpha*(1-m)/m;
    
    sigm_func = @(xdata,lin_slope,sig_slope,sig_centre) lin_slope*xdata./(1+exp(-sig_slope.*(xdata-sig_centre)));
    rev_sigm = @(n,a,b,c,x) n*(1-1./(1+exp(-a.*(x-c))))+b;  %% fit to average correlation / distance
    gauss_offs = @(n,s,mu,b,x) n*1/sqrt(2*pi*s.^2)*exp(-(x-mu).^2./(2*s.^2)) + b; %% fit to std(correlation) / distance
    
    model = struct;
    
    %%% ------------------------------- get joint model -------------------------------------
    
%      dist_array = [];
%      corr_array = [];
%      for s = 1:nSes
%        for sm = s+1:nSes
%          dist_array = [dist_array; xdata(s,sm).dist(xdata(s,sm).neighbours>0)];
%          corr_array = [corr_array; xdata(s,sm).fp_corr(xdata(s,sm).neighbours>0)];
%        end
%      end
    
%%% ------------------- estimate fitting functions for distances ---------------------------
    
    %% initialize fitting function
    F_dist = @(x,xdata)...
        x(1)*lognpdf(xdata,x(2),x(3))...                                %% lognormal distribution
        + x(4)*sigm_func(xdata,x(5),x(6),x(7));                                   %% linear * sigmoid function
    
    %% construct joint and normalized histogram
    histo.dist_tot = sum(histo.dist,2)/sum(histo.dist(:))*para.nbins/para.dist_max;
    
    %% get first estimate of parameters
    norm_0 = 1;                                                                         %% weight of both functions
    w_1 = 1;               %% weight of logn-function
    w_2 = 1;               %% weight of lin-sigmoid-function
    
    mu_guess = sum(histo.dist(:,1).*histo.dist_x) / sum(histo.dist(:,1));
    s_guess = sqrt( sum(histo.dist(:,1) .* (histo.dist_x-mu_guess).^2) / (sum(histo.dist(:,1))-1));
    mu_0 = logn_mu(mu_guess,s_guess);
    s_0 = logn_sigma(mu_guess,s_guess);
    
    m_0 = 1;                                                                            %% slope of sigmoid function
    c_0 = 0.75*para.dist_max;                                                           %% x value of half height of sigmoid function
    slope_0 = 1;
        
    p_init=[w_1 mu_0 s_0 w_2 slope_0 m_0 c_0];
    lb = [0 -Inf 0 0 0 0 0];                                                            %% lower bounds
    ub = [Inf Inf Inf Inf Inf Inf para.dist_max];                                                 %% upper bounds
    options = statset('MaxIter',1000, 'MaxFunEvals',2000,'Display','off');
    
    p_dist_final = lsqcurvefit(F_dist,p_init,histo.dist_x,histo.dist_tot,lb,ub,options);                %% fit data
    
    % calculating the distributions
    model.dist_same = p_dist_final(1)*lognpdf(histo.dist_x,p_dist_final(2),p_dist_final(3));
    model.dist_others = p_dist_final(4)*sigm_func(histo.dist_x,p_dist_final(5),p_dist_final(6),p_dist_final(7));
    model.dist_joint = model.dist_same + model.dist_others;
    
    model.w_dist_same = sum(model.dist_same)/sum(model.dist_same + model.dist_others);
    model.w_dist_others = sum(model.dist_others)/sum(model.dist_same + model.dist_others);
    
    model.dist_same = model.dist_same/sum(model.dist_joint)*para.nbins/para.dist_max;
    model.dist_others = model.dist_others/sum(model.dist_joint)*para.nbins/para.dist_max;
    model.dist_joint = model.dist_joint/sum(model.dist_joint)*para.nbins/para.dist_max;
    
    % finding the intersection between same cells and different cells:
    idx_search=find(histo.dist_x>1 & histo.dist_x<10);
    [~,idx_intersect]=min(abs(model.dist_same(idx_search)-model.dist_others(idx_search)));
    model.dist_intersect=round(100*histo.dist_x(idx_intersect+idx_search(1)-1))/100;      %% why round?
    
    % calculating the P_same of the model:
    model.p_same_dist = model.dist_same./model.dist_joint;
    model.p_same_dist(1)=model.p_same_dist(2); % avoid p_same going to 0 at 0 distance because of the lognormal distribution
    
    
%%% ------------------- estimate fitting functions for correlation ---------------------------
    
    %% construct joint histogram
    histo.corr_tot = sum(histo.corr,2)/sum(histo.corr(2:end))*para.nbins/para.corr_max;
%      histo.corr_tot(1) = 0;    %% remove peak at 0 correlation
    
    %% should implement restriction, that f_joint is normalized
    F_corr = @(x,xdata) ...
        x(1)*lognpdf(xdata,x(2),x(3)) ...          %% lognormal distribution (same cell)
        + x(4)*betapdf(xdata,x(5),x(6));           %% beta-distribution (other cells)    %% rather change to reverse sigmoid function??
%          + (1-x(4))*sigm_func(xdata,x(5),x(6)));                          %% sigmoid function
    
    w_1 = 1;
    w_2 = 1;
    
    mu_guess = sum(histo.corr(:,1).*histo.corr_x) / sum(histo.corr(:,1));
    s_guess = sqrt( sum(histo.corr(:,1) .* (histo.corr_x-mu_guess).^2) / (sum(histo.corr(:,1))-1));
    mu_0 = logn_mu(1-mu_guess,s_guess);
    s_0 = logn_sigma(1-mu_guess,s_guess);
    
    mu_guess = sum(histo.corr_x.*histo.corr(:,2))/sum(histo.corr(:,2));
    s_guess = sqrt( sum(histo.corr(:,2) .* (histo.corr_x-mu_guess).^2) ./ (sum(histo.corr(:,2))-1));
    alpha = beta_a(1-mu_guess,s_guess);
    beta = beta_b(1-mu_guess,alpha);
    
    lb = [0 -Inf 0 0 -Inf -Inf];
    ub = [Inf Inf Inf Inf Inf Inf];
    p_init = [w_1 mu_0 s_0 w_2 alpha beta];
    
    options = statset('MaxIter',1000, 'MaxFunEvals',2000,'Display','off');
    
    % finding the parameters that best fit the data:
    p_corr_final = lsqcurvefit(F_corr,p_init,1-histo.corr_x,histo.corr_tot,lb,ub,options);                %% fit data
    
    % calculating the distributions
    model.corr_same = p_corr_final(1)*lognpdf(1-histo.corr_x,p_corr_final(2),p_corr_final(3));
    model.corr_others = p_corr_final(4)*betapdf(1-histo.corr_x,p_corr_final(5),p_corr_final(6));
    model.corr_joint = model.corr_same + model.corr_others;
    
    model.w_corr_same = sum(model.corr_same)/sum(model.corr_joint);
    model.w_corr_others = sum(model.corr_others)/sum(model.corr_joint);
    
    model.corr_same = model.corr_same/sum(model.corr_joint)*para.nbins/para.corr_max;
    model.corr_others = model.corr_others/sum(model.corr_joint)*para.nbins/para.corr_max;
    model.corr_joint = model.corr_joint/sum(model.corr_joint)*para.nbins/para.corr_max;
    
    % finding the intersection between same cells and different cells:
    idx_search=find(histo.corr_x>0.3 & histo.corr_x<0.95);
    [~,idx_intersect]=min(abs(model.corr_same(idx_search)-model.corr_others(idx_search)));
    model.corr_intersect = round(100*histo.corr_x(idx_intersect+idx_search(1)-1))/100;
    
    % calculating the P_same of the model:
    model.p_same_corr = model.corr_same./model.corr_joint;
    model.p_same_corr(end) = 1;
    
    
    toc
    
%%% -------------------------------------------- compute the joint model -----------------------------------------------
    
    NN_mean = zeros(para.nbins,1);
    nNN_mean = zeros(para.nbins,1);
    
    NN_std = zeros(para.nbins,1);
    nNN_std = zeros(para.nbins,1);
    
    
    model.sub_same = zeros(para.nbins,para.nbins);
    model.sub_others = zeros(para.nbins,para.nbins);
    model.sub_total = zeros(para.nbins,para.nbins);
    
    model.p_same_joint = zeros(para.nbins,para.nbins);
    model.P_dcsame = zeros(para.nbins,para.nbins);
    model.P_dcothers = zeros(para.nbins,para.nbins);
    
    tStart = tic;
    
    hist_arr = [linspace(0,1-1./para.nbins,para.nbins) inf];
    
    %% calculate mean and variance of both populations
    NN_mean_0 = histo.joint(:,:,1)*histo.corr_x ./ sum(histo.joint(:,:,1),2);
    nNN_mean_0 = histo.joint(:,:,2)*histo.corr_x ./ sum(histo.joint(:,:,2),2);
    NN_mask = ~isnan(NN_mean_0);
    nNN_mask = ~isnan(nNN_mean_0);
    %% should be some kind of sigmoid function (same ROIs can only be so far, different ROIs can only be so close)
    rev_sigm_fit = fittype( rev_sigm );
    NN_mean_fit = fit(histo.dist_x(NN_mask),NN_mean_0(NN_mask),rev_sigm_fit,'Weights',histo.dist(NN_mask,1),'Lower',[0.5,0,0,0],'Upper',[1,Inf,0.5,para.dist_max],'StartPoint',[1,0.5,0,para.dist_max/2]);
    nNN_mean_fit = fit(histo.dist_x(nNN_mask),nNN_mean_0(nNN_mask),rev_sigm_fit,'Weights',histo.dist(nNN_mask,2),'Lower',[0.5,0,0,0],'Upper',[1,Inf,0.5,para.dist_max],'StartPoint',[1,0.5,0,para.dist_max/2]);
    model.NN_mean = rev_sigm(NN_mean_fit.n,NN_mean_fit.a,NN_mean_fit.b,NN_mean_fit.c,histo.dist_x);
    model.nNN_mean = rev_sigm(nNN_mean_fit.n,nNN_mean_fit.a,nNN_mean_fit.b,nNN_mean_fit.c,histo.dist_x);
    
    NN_std_0 = zeros(para.nbins,1);
    nNN_std_0 = zeros(para.nbins,1);
    for i = 1:para.nbins
      NN_std_0(i) = sqrt( histo.joint(i,:,1) * (histo.corr_x-NN_mean_0(i)).^2 / (sum(histo.joint(i,:,1))-1));
      nNN_std_0(i) = sqrt( histo.joint(i,:,2) * (histo.corr_x-nNN_mean_0(i)).^2 / (sum(histo.joint(i,:,2))-1));
    end
    NN_mask = ~isnan(NN_std_0);
    nNN_mask = ~isnan(nNN_std_0);
    %% should be some kind of offset normal distribution (low std at high and low correlation, as it's bounded)
    gauss_offs_fit = fittype( gauss_offs );
    NN_std_fit = fit(histo.dist_x(NN_mask),NN_std_0(NN_mask),gauss_offs_fit,'Weights',histo.dist(NN_mask,1),'Lower',[0 0 0 0],'Upper',[Inf para.dist_max/2 para.dist_max 1],'StartPoint',[nanmax(NN_std_0) para.dist_max/4 para.dist_max/2 0]);
    nNN_std_fit = fit(histo.dist_x(nNN_mask),nNN_std_0(nNN_mask),gauss_offs_fit,'Weights',histo.dist(nNN_mask,2),'Lower',[0 0 0 0],'Upper',[Inf para.dist_max/2 para.dist_max 1],'StartPoint',[nanmax(nNN_std_0) para.dist_max/4 para.dist_max/2 0]);
    model.NN_std = gauss_offs(NN_std_fit.n,NN_std_fit.s,NN_std_fit.mu,NN_std_fit.b,histo.dist_x);
    model.nNN_std = gauss_offs(nNN_std_fit.n,nNN_std_fit.s,nNN_std_fit.mu,nNN_std_fit.b,histo.dist_x);
    
    NN_frac = histo.joint(:,:,1)./sum(histo.joint,3);
    
    %% for each distance, obtain correlation model    %% (maybe better distinction (first guess) between NN and nNN?)
    for i = 1:para.nbins
      
      %% obtain (reverse) lognpdf for NN distribution 
      mu = logn_mu(1-model.NN_mean(i),model.NN_std(i));
      sigma = logn_sigma(1-model.NN_mean(i),model.NN_std(i));
      model.sub_same(i,:) = lognpdf(1-histo.corr_x,mu,sigma);
      model.sub_same(i,:) = model.sub_same(i,:)/sum(model.sub_same(i,:))*para.nbins/para.corr_max;  %% normalize
      model.P_dcsame(i,:) = model.sub_same(i,:)*model.dist_same(i)/model.w_dist_same;           %% P(dist,corr|same)
      
      %% obtain betapdf for nNN distribution
      alpha = beta_a(model.nNN_mean(i),model.nNN_std(i));
      beta = beta_b(model.nNN_mean(i),alpha);
      model.sub_others(i,:) = betapdf(histo.corr_x,alpha,beta);
      model.sub_others(i,:) = model.sub_others(i,:)/sum(model.sub_others(i,:))*para.nbins/para.corr_max;  %% normalize
      model.P_dcothers(i,:) = model.sub_others(i,:)*model.dist_others(i)/model.w_dist_others;     %% P(dist,corr|others)
      
    end
    
    %% and use as distributions to obtain the joint model
    model.p_same_joint = model.P_dcsame*model.w_dist_same ./ (model.P_dcsame*model.w_dist_same + model.P_dcothers*model.w_dist_others);
    
%      overall = model.P_dcsame*model.w_dist_same + model.P_dcothers*model.w_dist_others;
%      model.p_same_joint(model.P_dcsame < 10^(-3)) = 0;
    
    toc(tStart)
    
    
%      for i = 1:para.nbins
%        model.sub_same(i,:) = lognpdf(1-histo.corr_x,NN_mean(i),NN_var(i))
      
%        subplot(ceil(sqrt(para.nbins)),ceil(sqrt(para.nbins)),i)
      
%        hold on
%        b_corr = bar(histo.corr_x,[NN_histo,nNN_histo],1);
%        b_corr(1).FaceColor = 'green';
%        b_corr(2).FaceColor = 'red';
      
%        plot(histo.corr_x,model.sub_total(i,:),'k')
%        plot(histo.corr_x,model.sub_same(i,:),'green')
%        plot(histo.corr_x,model.sub_others(i,:),'r')
%          ylim([0,para.nbins/2])
%        hold off
%      end
    
    
    cmrg = [linspace(1,0,50)',linspace(0,1,50)',linspace(0,0,50)'];
    cmg = [linspace(0,0,50)',linspace(1,1,50)',linspace(0,0,50)'];
    cmr = [linspace(1,1,50)',linspace(0,0,50)',linspace(0,0,50)'];
    cmbw = [linspace(1,0,50)',linspace(1,0,50)',linspace(1,0,50)'];
    
    
    fig_histo = figure('position',[1000 100 1800 1500]);
    
%      subplot(3,1,1)
%      hold on
%      b_dist = bar(histo.dist_x,histo.dist,1);
%      b_dist(1).FaceColor = 'green';
%      b_dist(2).FaceColor = 'red';
%      
%      plot(histo.dist_x,model.dist_same*sum(histo.dist(:))*para.dist_max/para.nbins,'g')
%      plot(histo.dist_x,model.dist_others*sum(histo.dist(:))*para.dist_max/para.nbins,'r')
%      
%      hold off
%      xlim([0 para.dist_max])
    
    subplot(2,3,1)
    surf(histo.corr_x,histo.dist_x,histo.joint(:,:,1))
    caxis([0 max(max(histo.joint(:,:,1)))*0.05])
    colormap(cmrg)
    xlim([0 para.corr_max])
    ylim([0 para.dist_max])
    view(-151,44)
    title('data NN')
    
    subplot(2,3,2)
    surf(histo.corr_x,histo.dist_x,histo.joint(:,:,2))
    caxis([0 max(max(histo.joint(:,:,2)))*0.05])
    colormap(cmrg)
    xlim([0 para.corr_max])
    ylim([0 para.dist_max])
    view(-151,44)
    title('data others')
    
    ax_model = subplot(4,3,3);
    hold on
    plot(histo.dist_x,model.NN_mean,'g-')
    plot(histo.dist_x,NN_mean_0,'gx')
    plot(histo.dist_x,model.nNN_mean,'r-')
    plot(histo.dist_x,nNN_mean_0,'rx')
    xlim([0 para.dist_max])
    ylabel('mean')
    title('correlation mean model')
    hold off
    
    subplot(4,3,6)
    hold on
    plot(histo.dist_x,model.NN_std,'g-')
    plot(histo.dist_x,NN_std_0,'gx')
    plot(histo.dist_x,model.nNN_std,'r-')
    plot(histo.dist_x,nNN_std_0,'rx')
    xlim([0 para.dist_max])
    xlabel('centroid distance')
    ylabel('STD')
    title('correlation STD model')
    hold off
    
    
    ax_same = subplot(2,3,4);
    surf(histo.corr_x,histo.dist_x,model.sub_same)
    caxis([0 max(model.sub_same(:))*0.05])
    colormap(ax_same,cmrg)
    xlim([0 para.corr_max])
    ylim([0 para.dist_max])
    view(-151,44)
    title('model NN (normalized)')
    
    ax_others = subplot(2,3,5);
    surf(histo.corr_x,histo.dist_x,model.sub_others)
    caxis([0 max(model.sub_others(:))*0.05])
    colormap(ax_others,cmrg)
    xlim([0 para.corr_max])
    ylim([0 para.dist_max])
    view(-151,44)
    title('model others (normalized)')
    
    ax4 = subplot(2,3,6);
    surf(histo.corr_x,histo.dist_x,model.p_same_joint)
    colormap(ax4,cmrg)
    xlim([0,1])
    ylim([0,para.dist_max])
    zlim([0,1])
    view([-151,44])
    xlabel('footprint correlation','rotation',-10)
    ylabel('centroid distance','rotation',60)
    zlabel('P_{same}')
    
    path = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/joint_model.png';
    print(path,'-dpng','-r300')
%      saveas(fig_model,path,'jpg')
    
    disp(sprintf('saved under %s',path))
    
    
    
    
    fig_model = figure('position',[100 100 1200 1000]);
    
    %% set up alpha value to display "amount of data per pixel"
    I = log(sum(histo.joint,3));
    I = I/max(I(:));
    
%      ax1 = subplot(2,2,1);
%      pos = get(gca,'Position');
%      h_NNfrac = imagesc(ax1,histo.corr_x,histo.dist_x,NN_frac);
%      colormap(ax1,cmrg)
%      set(ax1,'Ydir','normal')
%      cb1 = colorbar(ax1,'horiz','Position',[pos(1) pos(2)+pos(4)+.01 pos(3) 0.02]);
%      set(h_NNfrac,'AlphaData',I);
%      cb1.Label.String = 'Fraction of nearest neighbours';
%      ylabel('centroid distance')
%      
%      ax_model = subplot(4,2,2);
%      hold on
%      plot(histo.dist_x,model.NN_mean,'g-')
%      plot(histo.dist_x,NN_mean_0,'gx')
%      plot(histo.dist_x,model.nNN_mean,'r-')
%      plot(histo.dist_x,nNN_mean_0,'rx')
%      xlim([0,para.dist_max])
%      ylabel('mean')
%      title('correlation mean model')
%      hold off
%      
%      subplot(4,2,4)
%      hold on
%      plot(histo.dist_x,model.NN_std,'g-')
%      plot(histo.dist_x,NN_std_0,'gx')
%      plot(histo.dist_x,model.nNN_std,'r-')
%      plot(histo.dist_x,nNN_std_0,'rx')
%      xlim([0,para.dist_max])
%      xlabel('centroid distance')
%      ylabel('STD')
%      title('correlation STD model')
%      hold off
%      
%      ax_pjoint = subplot(2,2,3);
    ax_pjoint = gca;
    pos = get(ax_pjoint,'Position');    
    h_pjoint = imagesc(ax_pjoint,histo.corr_x,histo.dist_x,model.p_same_joint);
    set(h_pjoint,'AlphaData',I);
    colormap(ax_pjoint,cmrg)
    cb21 = colorbar(ax_pjoint,'horiz','Position',[pos(1) pos(2)+pos(4)+.01 pos(3) 0.02]);
    cb21.Label.String = 'P_{same}';
    xlabel('footprint correlation')
    ylabel('centroid distance')
    
    ax_contour = axes('position',pos);
    contour(ax_contour,histo.corr_x,histo.dist_x,model.p_same_joint,[0.05,0.1,0.3,0.5,0.7,0.9,0.95])
    ax_contour.Visible = 'off';
    ax_contour.XTick = [];
    ax_contour.YTick = [];
    linkaxes([ax_pjoint,ax_contour])
    colormap(ax_contour,cmbw)
%      cb22 = colorbar(ax_contour,'horiz','Position',[pos(1) pos(2)+pos(4)+.02 pos(3) 0.02]);
    cb22.Label.String = 'P_{same}';
    set([ax_pjoint,ax_contour],'YDir','normal')
    
%      set(ax_contour,'YDir','normal')
    title('P_{same}(dist,corr)')
    
    
%      ax4 = axes('Position',[0.575 0.1 0.3 0.35],'Parent',fig_model);
%      surf(histo.corr_x,histo.dist_x,model.p_same_joint)
%      colormap(ax4,cmrg)
%      xlim([0,1])
%      ylim([0,para.dist_max])
%      zlim([0,1])
%      view([-160,35])
%      xlabel('footprint correlation','rotation',-10)
%      ylabel('centroid distance','rotation',60)
%      zlabel('P_{same}')
    
    
%      I = log(histo.joint(:,:,1));
%      I = I/max(I(:));
%      ax3 = subplot(2,2,3);
%      pos = get(gca,'Position');
%      h_Psame = imagesc(ax3,histo.corr_x,histo.dist_x,log(model.P_dcsame),[-5,0]);% max(model.P_dcsame(:))]);
%      colormap(ax3,cmg)
%      cb3 = colorbar(ax3,'horiz','Position',[pos(1) pos(2)-0.07 pos(3) 0.02]);
%      set(h_Psame,'AlphaData',I);
%  %      title('P(dist,corr|same)')
%      
%      I = log(histo.joint(:,:,2));
%      I = I/max(I(:));
%      
%      ax4 = axes('Position',pos);%subplot(2,2,4)
%      h_Pothers = imagesc(ax4,histo.corr_x,histo.dist_x,log(model.P_dcothers),[-5,0]);% max(model.P_dcothers(:))]);
%      ax4.Visible = 'off';
%      ax4.XTick = [];
%      ax4.YTick = [];
%      linkaxes([ax3,ax4])
%      set(h_Pothers,'AlphaData',I);
%      colormap(ax4,cmr)
%      set([ax3,ax4],'YDir','normal')
%  %      cb4 = colorbar(ax4,'horiz','Position',[pos(1) pos(2)-0.07 pos(3) 0.02]);
%  %      title('P(dist,corr|others)')
    
    path = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/joint_model_contour.png';
    print(path,'-dpng','-r300')
%      saveas(fig_model,path,'jpg')
    
    disp(sprintf('saved under %s',path))
    
end
