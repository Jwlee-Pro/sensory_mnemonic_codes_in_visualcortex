clear; close all;
set(0, 'defaultfigurewindowstyle', 'docked')

addpath('./Package_JS')
load('./Tab_all.mat')

ID = zeros(size(Tab,1),1);
for i = 1:size(Tab,1) 
   ID(i) = str2double(Tab.Subj(i,:));
end

Tab.ID = ID;

IDs = unique(ID);
nSbj = length(IDs);

Index = table(IDs);


%% make bias corrected estimation error

Tab.error_corrected = nan(size(Tab,1),1);

card_filter = repelem([1,-1,1,-1], 6);

for sbj = 1:nSbj
    idx = ID == IDs(sbj);
    T = Tab(idx,:);
    %plot(T.stim, T.error, 'o')
    x = T.stim;
    y = T.error;
    
    yy1 = smooth(x,y, 0.2, 'rloess');
    yy2 = smooth(x,y, 0.01, 'rloess'); 
    
    % using robust (outlier weight 0) lowess with almost no span (not like smoothing but more like conditional mean)
    yy  = yy2;
    Tab.error_corrected(idx) = y-yy;
    
    % quantifying individual's oblique bias
    ux = unique(x);
    uyy = ux*nan;
    for i = 1:length(ux)
       uyy(i) = nanmean(yy(x == ux(i)));
    end
    
    uyy = uyy-mean(uyy);
    
    Index.oblique_bias(sbj) = mean(uyy.*card_filter');
    
    % quantifying individual's estimation density categorization (by cardinal)
    stim = T.stim;
    est = T.esti;
    dist_from_card = min([abs(est-0), abs(est-90), abs(est-180)], [], 2);
    d1 = nansum(dist_from_card<22.5);
    d2 = nansum(dist_from_card>=22.5);
    
    Index.density_cate(sbj) = (d2-d1)/(d1+d2);
    
    
    
    if false
        figure(sbj);
        clf
        subplot(2,1,1);
        plot(x,y, 'o'); hold on;
        [sx, ord] = sort(x);
        sy1 = yy1(ord);
        sy2 = yy2(ord);
        plot(sx,sy1, 'color', 'red');
        plot(sx,sy2, 'color', 'green');
        title(IDs(sbj));
        ylim([-35, 35]);
        hold off;
        
        X = [T.esti; T.esti+180; T.esti-180];
        [f,xi] = ksdensity(X, 'function', 'pdf', 'bandwidth', 5);
        subplot(2,2,3);
        cla;
        histogram(T.esti, 24, 'Normalization', 'probability'); hold on;
        plot(xi, f*30)
        xlim([0, 180]);
        
        subplot(2,2,4);
        histogram(dist_from_card, 20)
        
    end
    
end

% erase stimulus driven bias




%% Index for precision and estimation precision

for sbj = 1:nSbj
    idx = ID == IDs(sbj);
    T = Tab(idx,:);
    
    q = 180/pi;
    %% decision precision indexing
    ref = T.ref;
    choice = T.choice == 2;
    timing = T.Timing;
    
    fxn1 = @(par) -LL_cdf_choice(par, ref, choice);
    par0 = [0,5,0.01, 0.01];
    par_fitted = fminsearchbnd(fxn1, par0, [-90, 0.01, 0, 0], [90, 100, 0.2, 0.2]);
    
    Index.decision_prec(sbj) = q/par_fitted(2);
    
    fxn2 = @(par) -LL_cdf_choice(par([1,2,4,5]), ref(timing == 1), choice(timing == 1))...
        -LL_cdf_choice(par([1,3,4,5]), ref(timing == 2), choice(timing == 2));
    par0 = [0, 5, 5, 0.01, 0.01];
    par_fitted2 = fminsearchbnd(fxn2, par0, [-90, 0.01,0.01, 0, 0], [90, 100, 100, 0.2, 0.2]);
    Index.decision_prec_decay(sbj) = par_fitted2(3)/par_fitted2(2);
    
    
    if false
        
        tmp = {choice(ref == -21), choice(ref == -4),...
            choice(ref == 0), choice(ref == 4),choice(ref == 21)};
        pm = cellfun(@nanmean, tmp);
        xs = linspace(-25, 25, 100);
        ys = normcdf(xs, par_fitted(1), par_fitted(2))*(1-par_fitted(3)-par_fitted(4)) + par_fitted(3);
        figure(sbj); clf;
        plot([-21, -4, 0, 4, 21], pm, 'o'); hold on;
        plot(xs, ys);
        title(IDs(sbj));

    end
    
    % estimation precision quantified by IQR of estimation error
    % distribution
    
    Index.est_prec(sbj) = 90/iqr(T.error);
    Index.est_prec_cor(sbj) = 90/iqr(T.error_corrected);
    
    %% DIV and TIV (reference distnace induced variance and timing induced variance)
    
    err = T.error;
    err_cor = T.error_corrected;
    
    ref = T.ref;
    timing = T.Timing;
    
    idx_r1 = abs(ref)<10;
    idx_r2 = abs(ref)>10;
    
    clf;
    histogram(err(idx_r1)); hold on;
    histogram(err(idx_r2))
    clf;
    histogram(err_cor(idx_r1)); hold on;
    histogram(err_cor(idx_r2))
    
    idx_t1 = timing == 1;
    idx_t2 = timing == 2;
    
    Index.DIV(sbj) = iqr(err(idx_r1))/iqr(err(idx_r2));
    Index.DIV_cor(sbj) = iqr(err_cor(idx_r1))/iqr(err_cor(idx_r2));
    Index.TIV(sbj) = iqr(err(idx_t2))/iqr(err(idx_t1));
    Index.TIV_cor(sbj) = iqr(err_cor(idx_t2))/iqr(err_cor(idx_t1));
    
    %% Decision-Estimation consistency index
    err = T.error;
    ref = T.ref;
    choice = T.choice == 1;
    err_from_ref = err-ref;
    
    idx_r1 = abs(ref)<10;
    x = err_from_ref(idx_r1);
    y = choice(idx_r1);
    
    tmp = ~isnan(x+y);
    x = x(tmp);
    y = y(tmp);
    
    fxn1 = @(par) -LL_cdf_choice(par, x, y);
    par0 = [0,5,0.01, 0.01];
    par_fitted = fminsearchbnd(fxn1, par0, [-90, 0.01, 0, 0], [90, 100, 0.2, 0.2]);
    sig1 = par_fitted(2);
    
    fxn2 = @(par) -LL_pdf_est(par,T.error_corrected);
    par0 = [0,5, 0];
    par_fitted2 = fminsearchbnd(fxn2, par0, [-90, 0.01, 0], [90, Inf, 0.3]);
    sig2 = par_fitted2(2);
    
    plot(x, y, 'o'); hold on;
    xs = linspace(-25, 25, 100);
    ys = normcdf(xs, par_fitted(1), par_fitted(2))*(1-par_fitted(3)-par_fitted(4)) + par_fitted(3);
    plot(xs, ys);
    
    Index.DE_con(sbj) = q/sig1;
    Index.DE_con_rel(sbj) = 1-sig1/sig2;
    
    
    
end

save('./Index.mat', 'Tab', 'Index');

% clf
% plot(Index.DIV, Index.DIV_cor, 'o')
% plot(Index.TIV, Index.TIV_cor, 'o')
% plot(T.error, T.error_corrected, 'o');











