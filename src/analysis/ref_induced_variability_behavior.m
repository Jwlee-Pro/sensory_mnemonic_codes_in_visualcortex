clear all; close all; clc;

% Load sub/roi-list 
load('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/sub_list.mat')
addpath('/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/src/packages/CircStat2012a')
addpath('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/B1_IEM_revisit/customcolormap')

nTR = 14; 
TRinterest = [3 4; 6 7; 9 10; 12 13]; 
stimcond = 0:7.5:179; 


% Colormaps 
st_color = jet(length(stimcond))*0.8 ; 
a = linspace(0,1,8); 
ROI_color = customcolormap(a(1:8), [96 182 110; 146 195 68; 245 233 42; 235 143 49; 217 40 44; 121 77 142; 49 73 140; 13 132 134]/255); 


%% Main parameters
runSub = 1:50; 
nTR = 14; 
refs = [-21, -4, 0, 4, 21]; 
stimcond = 0:7.5:172.5; 

binCenter = -90:5:85 ; 
nBin = length(binCenter); 

matSub = {}; 
matSub.mean = nan(length(stimcond), 2, length(sub_list)); 
matSub.medi = nan(length(stimcond), 2, length(sub_list)); 
matSub.std  = nan(length(stimcond), 2, length(sub_list)); 
matSub.iqr  = nan(length(stimcond), 2, length(sub_list)); 

stim_m = []; 
timing_m = []; 
ref_m = []; 
choice_m = []; 
err_m = []; 
nObs = nan(length(stimcond), 2, length(sub_list)); 
% dec_tr = []; 
% dec_sen = []; 
% dec_mne = []; 
lapse_crit = 30 ;

for isub = runSub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_estimated/VC_sub-' sub_list(isub,:) '_dec.mat'], 'response','choice','stimulus','ref','timing')
    
    % Behavioral errors
    errme = response - stimulus; 
    errme(errme>90) = errme(errme>90) -180; 
    errme(errme<-90) = errme(errme<-90) +180; 
    
    for istim = 1:length(stimcond)
        for ic = 1:2
            ind = find(stimulus==stimcond(istim) & choice==ic & ~isnan(response) & abs(errme)<lapse_crit); 
            
            % Calculate circular mean (or median) and std (or iqr) 
            matSub.mean(istim,ic,isub) = circ_mean(errme(ind)'*2*pi/180)*180/pi/2; 
            matSub.medi(istim,ic,isub) = circ_median(errme(ind)'*2*pi/180)*180/pi/2; 
            matSub.std(istim,ic,isub)  = circ_std(errme(ind)'*2*pi/180)*180/pi/2; 
            matSub.iqr(istim,ic,isub)  = iqr(errme(ind)); 
            
            % Number of samples to calculate 'mean' and 'iqr'
            nObs(istim,ic,isub) = length(ind); 
        end
    end
    
    % All merged
    stim_m = [stim_m stimulus];
    timing_m = [timing_m timing]; 
    ref_m = [ref_m ref]; 
    choice_m = [choice_m choice]; 
    err_m = [err_m errme]; 
end

% [merged] Calculate circular mean (or median) and std (or iqr) 
matMerge = {}; 
matMerge.mean = nan(length(stimcond), 2); 
matMerge.medi = nan(length(stimcond), 2); 
matMerge.std  = nan(length(stimcond), 2); 
matMerge.iqr  = nan(length(stimcond), 2);

for istim = 1:length(stimcond)
    for ic = 1:2
        ind = find(stim_m==stimcond(istim) & choice_m==ic & ~isnan(err_m) & abs(err_m)<lapse_crit); 

        matMerge.mean(istim,ic) = circ_mean(err_m(ind)'*2*pi/180)*180/pi/2; 
        matMerge.medi(istim,ic) = circ_median(err_m(ind)'*2*pi/180)*180/pi/2; 
        matMerge.std(istim,ic)  = circ_std(err_m(ind)'*2*pi/180)*180/pi/2; 
        matMerge.iqr(istim,ic)  = iqr(err_m(ind)); 
    end
end

figure(99); clf; 
SP = subplot(1,1,1); clf; 
hist(err_m(abs(err_m)<lapse_crit),40)
xlim([-45 45]); 
xlabel('error (deg)'); 


%% Figure drawing
set(figure(100),'position',[1 436 264 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matMerge.mean(:,1),'r.-','linewidth',1.3,'markersize',12); 
plot(stimcond, matMerge.mean(:,2),'b.-','linewidth',1.3,'markersize',12); 
xlim([0 180]); ylim([-12 12]); 
xticks(linspace(0,180,5)); 
yticks(linspace(-10, 10,3)); 
xlabel('Orientation'); ylabel('Error bias (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matMerge.mean(:,1)-matMerge.mean(:,2),'k.-','linewidth',1.3,'markersize',12); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('\Delta Error bias (deg)'); 

set(figure(101),'position',[1 436 264 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matMerge.medi(:,1),'r.-','linewidth',1.3,'markersize',12); 
plot(stimcond, matMerge.medi(:,2),'b.-','linewidth',1.3,'markersize',12); 
xlim([0 180]); ylim([-12 12]); 
xticks(linspace(0,180,5)); 
yticks(linspace(-10, 10,3)); 
xlabel('Orientation'); ylabel('Error bias (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matMerge.medi(:,1)-matMerge.medi(:,2),'k.-','linewidth',1.3,'markersize',12); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('\Delta Error bias (deg)'); 


set(figure(102),'position',[1 1 264 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matMerge.std(:,1),'r.-','linewidth',1.3,'markersize',12); 
plot(stimcond, matMerge.std(:,2),'b.-','linewidth',1.3,'markersize',12); 
xlim([0 180]); %ylim([-12 12]); 
xticks(linspace(0,180,5)); 
% yticks(linspace(-10, 10,3)); 
xlabel('Orientation'); ylabel('Error std (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matMerge.std(:,1)-matMerge.std(:,2),'k.-','linewidth',1.3,'markersize',12); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('\Delta Error std (deg)'); 


set(figure(103),'position',[1 1 264 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matMerge.iqr(:,1),'r.-','linewidth',1.3,'markersize',12); 
plot(stimcond, matMerge.iqr(:,2),'b.-','linewidth',1.3,'markersize',12); 
xlim([0 180]); %ylim([-12 12]); 
xticks(linspace(0,180,5)); 
% yticks(linspace(-10, 10,3)); 
xlabel('Orientation'); ylabel('Error iqr (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matMerge.iqr(:,1)-matMerge.iqr(:,2),'k.-','linewidth',1.3,'markersize',12); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('\Delta Error iqr (deg)'); 



set(figure(104),'position',[1 565 263 240]); clf;
plot(matMerge.mean(:), matMerge.medi(:),'ko')
xlabel('circ mean'); ylabel('circ median'); 


set(figure(105),'position',[1 565 263 240]); clf;
plot(matMerge.std(:), matMerge.iqr(:),'ko')
xlabel('circ std'); ylabel('iqr'); 


%% Average across near/far
% [merged] Calculate circular mean (or median) and std (or iqr) 
matNearFar = {}; 
matNearFar.mean = nan(length(stimcond), 2); 
matNearFar.medi = nan(length(stimcond), 2); 
matNearFar.std  = nan(length(stimcond), 2); 
matNearFar.iqr  = nan(length(stimcond), 2);
err_m_cor = err_m; 
for istim = 1:length(stimcond)
    % Calculate std 'after' subtracting mean term
    for ir = 1:length(refs)
        ind1 = find(stim_m==stimcond(istim) & ref_m==refs(ir) & ~isnan(err_m) & abs(err_m)<lapse_crit); 
        temp = err_m(ind1) - circ_mean(err_m(ind1)'*2*pi/180)*180/pi/2; 
        temp(temp>90) = temp(temp>90) -180; 
        temp(temp<-90) = temp(temp<-90) +180;  
        err_m_cor(ind1) = temp; 
    end
    for ic = 1:2
        if ic == 1  % Near
            ref_ind = abs(ref_m)<5; 
        else        % Far
            ref_ind = abs(ref_m)>5; 
        end
        ind = find(stim_m==stimcond(istim) & ref_ind & ~isnan(err_m) & abs(err_m)<lapse_crit); 

        matNearFar.mean(istim,ic) = circ_mean(err_m(ind)'*2*pi/180)*180/pi/2; 
        matNearFar.medi(istim,ic) = circ_median(err_m(ind)'*2*pi/180)*180/pi/2; 
        
        
        matNearFar.std(istim,ic)  = circ_std(err_m_cor(ind)'*2*pi/180)*180/pi/2; 
        matNearFar.iqr(istim,ic)  = iqr(err_m_cor(ind)); 
    end
end


%% Figure drawing (near vs. far)

nearfar_color = [170 121 66; 248 186 0; ]/255; 

set(figure(106),'position',[1 436 264 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matNearFar.mean(:,1),'r.-','linewidth',1.3,'markersize',12,'color',nearfar_color(1,:)); 
plot(stimcond, matNearFar.mean(:,2),'b.-','linewidth',1.3,'markersize',12,'color',nearfar_color(2,:)); 
xlim([0 180]); ylim([-12 12]); 
xticks(linspace(0,180,5)); 
yticks(linspace(-10, 10,3)); 
xlabel('Orientation'); ylabel('Error bias (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matNearFar.mean(:,1)-matNearFar.mean(:,2),'k.-','linewidth',1.3,'markersize',12); 
plot([0 180],[0 0],'k--'); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('\Delta Error bias (deg)'); 

set(figure(107),'position',[1 436 264 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matNearFar.medi(:,1),'r.-','linewidth',1.3,'markersize',12,'color',nearfar_color(1,:)); 
plot(stimcond, matNearFar.medi(:,2),'b.-','linewidth',1.3,'markersize',12,'color',nearfar_color(2,:)); 
xlim([0 180]); ylim([-12 12]); 
xticks(linspace(0,180,5)); 
yticks(linspace(-10, 10,3)); 
xlabel('Orientation'); ylabel('Error bias (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matNearFar.medi(:,1)-matNearFar.medi(:,2),'k.-','linewidth',1.3,'markersize',12); 
plot([0 180],[0 0],'k--'); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('\Delta Error bias (deg)'); 


set(figure(108),'position',[1 1 264 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matNearFar.std(:,1),'r.-','linewidth',1.3,'markersize',12,'color',nearfar_color(1,:)); 
plot(stimcond, matNearFar.std(:,2),'b.-','linewidth',1.3,'markersize',12,'color',nearfar_color(2,:)); 
xlim([0 180]); %ylim([-12 12]); 
xticks(linspace(0,180,5)); 
% yticks(linspace(-10, 10,3)); 
xlabel('Orientation'); ylabel('Error std (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matNearFar.std(:,1)-matNearFar.std(:,2),'k.-','linewidth',1.3,'markersize',12); 
plot([0 180],[0 0],'k--'); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('\Delta Error std (deg)'); 


set(figure(109),'position',[1 1 264 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matNearFar.iqr(:,1),'r.-','linewidth',1.3,'markersize',12,'color',nearfar_color(1,:)); 
plot(stimcond, matNearFar.iqr(:,2),'b.-','linewidth',1.3,'markersize',12,'color',nearfar_color(2,:)); 
xlim([0 180]); %ylim([-12 12]); 
xticks(linspace(0,180,5)); 
% yticks(linspace(-10, 10,3)); 
xlabel('Orientation'); ylabel('Error iqr (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matNearFar.iqr(:,1)-matNearFar.iqr(:,2),'k.-','linewidth',1.3,'markersize',12); 
plot([0 180],[0 0],'k--'); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('\Delta Error iqr (deg)'); 



set(figure(110),'position',[1 565 263 240]); clf;
plot(matNearFar.mean(:), matNearFar.medi(:),'ko')
xlabel('circ mean'); ylabel('circ median'); 


set(figure(111),'position',[1 565 263 240]); clf;
plot(matNearFar.std(:), matNearFar.iqr(:),'ko')
xlabel('circ std'); ylabel('iqr'); 


set(figure(112),'position',[1 1 172 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(matNearFar.medi(:,1),matNearFar.medi(:,2),'ko')
xlabel('medi (near)'); ylabel('medi (far)'); 
plot([-10 10],[-10 10],'k--'); 
xlim([-10 10]);
ylim([-10 10]);
[h,p,ci,stats] = ttest(matNearFar.mean(:,1)-matNearFar.mean(:,2)); 
title(['p = ' num2str(p)]); 

SP = subplot(2,1,2); cla; hold on; 
plot(matNearFar.iqr(:,1),matNearFar.iqr(:,2),'ko')
xlabel('iqr (near)'); ylabel('iqr (far)'); 
plot([5 18],[5 18],'k--'); 
xlim([5 18]);
ylim([5 18]);
[h,p,ci,stats] = ttest(matNearFar.iqr(:,1)-matNearFar.iqr(:,2)); 
title(['p = ' num2str(p)]); 


set(figure(113),'position',[1 1 172 369]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(matMerge.medi(:,1),matMerge.medi(:,2),'ko')
xlabel('medi (ccw)'); ylabel('medi (cw)'); 
plot([-10 10],[-10 10],'k--'); 
xlim([-10 10]);
ylim([-10 10]);
[h,p,ci,stats] = ttest(matMerge.mean(:,1)-matMerge.mean(:,2)); 
title(['p = ' num2str(p)]); 

SP = subplot(2,1,2); cla; hold on; 
plot(matMerge.iqr(:,1),matMerge.iqr(:,2),'ko')
xlabel('iqr (ccw)'); ylabel('iqr (cw)'); 
plot([5 18],[5 18],'k--'); 
xlim([5 18]);
ylim([5 18]);
[h,p,ci,stats] = ttest(matMerge.iqr(:,1)-matMerge.iqr(:,2)); 
title(['p = ' num2str(p)]); 












