clear all; close all; clc;

%% Main parameters
% Load sub/roi-list 
load('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/sub_list.mat')
addpath('/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/src/packages/CircStat2012a')
addpath('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/B1_IEM_revisit/customcolormap')
addpath('/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/src/library'); 

nTR = 14; 
TRinterest = [3 4; 6 7; 9 10; 12 13]; 
stimcond = 0:7.5:179; 


% Colormaps 
st_color = jet(length(stimcond))*0.8 ; 
a = linspace(0,1,8); 
ROI_color = customcolormap(a(1:8), [96 182 110; 146 195 68; 245 233 42; 235 143 49; 217 40 44; 121 77 142; 49 73 140; 13 132 134]/255); 
J_el = [102, 126, 182; 182, 126, 102]/255; 
nearfar_color = [170 121 66; 248 186 0; ]/255; 

runSub = 1:50; 
nTR = 14; 
refs = [-21, -4, 0, 4, 21]; 
stimcond = 0:7.5:172.5; 

binCenter = -90:5:85 ; 
nBin = length(binCenter); 


%% Main 
stim_m = []; 
timing_m = []; 
ref_m = []; 
choice_m = []; 
err_m = []; 
dec_mne = []; 

nObs = nan(length(stimcond), 2, length(sub_list)); 

lapse_crit = 20 ;

errsize_e = nan(nTR, length(sub_list)); 
errsize_l = nan(nTR, length(sub_list)); 

for isub = runSub
    % Memory-code (decoded)
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_earlydelay/VC_sub-' sub_list(isub,:) '_dec.mat'])
    [nTR, nTrials] = size(Decoded_result.est);
    
    % Error (Behavior)
    errme = response - stimulus; 
    errme(errme>90) = errme(errme>90) -180; 
    errme(errme<-90) = errme(errme<-90) +180; 
    
    % Orientation-specific features
    for istim = 1:length(stimcond)
        ind = find(stimulus==stimcond(istim) & ~isnan(response) & abs(errme)<lapse_crit); 

        % Calculate circular mean and std 
        matBeh.ori.mean(istim,isub) = circ_m(errme(ind)'); 
        matBeh.ori.std(istim,isub)  = circ_s(errme(ind)'); 

        % Number of samples to calculate mean and std 
        nObs_stim(istim,isub) = length(ind); 
    end
    % Reference-specific features
    for iref = 1:length(refs)
        ind = find(ref==refs(iref) & ~isnan(response) & abs(errme)<lapse_crit); 

        % Calculate circular mean and std 
        matBeh.ref.mean(iref,isub) = circ_m(errme(ind)'); 
        matBeh.ref.std(iref,isub)  = circ_s(errme(ind)'); 

        % Number of samples to calculate mean and std 
        nObs_ref(iref,isub) = length(ind); 
    end
    
    % Error (fmri decoded)
    errme_fmri = nan(size(Decoded_result.est)); 
    for iTR = 1:nTR
        temp = Decoded_result.est(iTR,:) - stimulus; 
        temp(temp>90) = temp(temp>90) -180; 
        temp(temp<-90) = temp(temp<-90) +180; 
        errme_fmri(iTR,:) = temp ;
    end
    
    % Calculate mean & std for each subject 
    for istim = 1:length(stimcond)
        % All trials
        ind = find(stimulus==stimcond(istim) & ~isnan(response) & abs(errme)<lapse_crit); 
        for iTR = 1:nTR
            matfMRI.ori.mean(iTR,istim,isub) = circ_m(errme_fmri(iTR,ind)'); 
            matfMRI.ori.std(iTR,istim,isub)  = circ_s(errme_fmri(iTR,ind)'); 
        end
        
        % Early Dm trials
        ind = find(stimulus==stimcond(istim) & ~isnan(response) & abs(errme)<lapse_crit & timing==1); 
        for iTR = 1:nTR
            matfMRI.ori.mean_e(iTR,istim,isub) = circ_m(errme_fmri(iTR,ind)'); 
            matfMRI.ori.std_e(iTR,istim,isub)  = circ_s(errme_fmri(iTR,ind)'); 
        end
        
        % Late Dm trials
        ind = find(stimulus==stimcond(istim) & ~isnan(response) & abs(errme)<lapse_crit & timing==2); 
        for iTR = 1:nTR
            matfMRI.ori.mean_l(iTR,istim,isub) = circ_m(errme_fmri(iTR,ind)'); 
            matfMRI.ori.std_l(iTR,istim,isub)  = circ_s(errme_fmri(iTR,ind)'); 
        end
    end
    % Reference-specific features
    for iref = 1:length(refs)
        % All trials
        ind = find(ref==refs(iref) & ~isnan(response) & abs(errme)<lapse_crit); 
        for iTR = 1:nTR
            matfMRI.ref.mean(iTR,iref,isub) = circ_m(errme_fmri(ind)'); 
            matfMRI.ref.std(iTR,iref,isub)  = circ_s(errme_fmri(ind)'); 
        end
        
        % Early Dm trials
        ind = find(ref==refs(iref) & ~isnan(response) & abs(errme)<lapse_crit); 
        for iTR = 1:nTR
            matfMRI.ref.mean_e(iTR,iref,isub) = circ_m(errme_fmri(ind)'); 
            matfMRI.ref.std_e(iTR,iref,isub)  = circ_s(errme_fmri(ind)'); 
        end
        
        % Late Dm trials
        ind = find(ref==refs(iref) & ~isnan(response) & abs(errme)<lapse_crit); 
        for iTR = 1:nTR
            matfMRI.ref.mean_l(iTR,iref,isub) = circ_m(errme_fmri(ind)'); 
            matfMRI.ref.std_l(iTR,iref,isub)  = circ_s(errme_fmri(ind)'); 
        end
    end
    
    
    % Trial merging
    stim_m = [stim_m stimulus];
    timing_m = [timing_m timing]; 
    ref_m = [ref_m ref]; 
    choice_m = [choice_m choice]; 
    err_m = [err_m errme]; 
    dec_mne = [dec_mne errme_fmri];
    
    errsize_e(:,isub) = 45-circ_mean(abs(errme_fmri(:,timing==1))'*2*pi/180)*180/pi/2; 
    errsize_l(:,isub) = 45-circ_mean(abs(errme_fmri(:,timing==2))'*2*pi/180)*180/pi/2; 
    
end

%% Figure: orientation-specific bias/variance
set(figure(1),'position',[1 1169 398 176]); clf; 
SP = subplot(1,2,1); cla; hold on; 
errorbar(stimcond,circ_m(matBeh.ori.mean'),circ_s(matBeh.ori.mean')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
plot([0 180],[0 0]+circ_m(circ_m(matBeh.ori.mean')'),'k--'); 
xlabel('orientation (deg)'); ylabel('bias (deg)');
xlim([0 180]); 

SP = subplot(1,2,2); cla; hold on; 
errorbar(stimcond,circ_m(matBeh.ori.std'),circ_s(matBeh.ori.std')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
xlabel('orientation (deg)'); ylabel('variability (iqr)');
xlim([0 180]); 


set(figure(2),'position',[1 1169 398 176]); clf; 
SP = subplot(1,2,1); cla; hold on; 
errorbar(refs,circ_m(matBeh.ref.mean'),circ_s(matBeh.ref.mean')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
plot([-22, 22],[0 0]+circ_m(circ_m(matBeh.ref.mean')'),'k--'); 
xlabel('reference'); ylabel('bias (deg)');
xlim([-25, 25]); 

SP = subplot(1,2,2); cla; hold on; 
errorbar(refs,circ_m(matBeh.ref.std'),circ_s(matBeh.ref.std')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
xlabel('reference'); ylabel('variability (iqr)');
xlim([-25, 25]); 

%% Figure: orientation-specific bias/variance (fMRI)
set(figure(11),'position',[1 1068 1432 277]); clf; 
for iTR = 1:nTR
    SP = subplot(2,9,iTR); cla; hold on; 
    errorbar(stimcond,circ_m(squeeze(matfMRI.ori.mean(iTR,:,:))'),circ_s(squeeze(matfMRI.ori.mean(iTR,:,:))')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1,'markersize',3); 
%     plot([0 180],[0 0]+circ_m(circ_m(matBeh.ori.mean')'),'k--'); 
    xlabel('orientation (deg)'); ylabel('bias (deg)'); 
    title(['TR=' num2str(iTR)]); 
    xlim([0 180]); 
    xticks(linspace(0,180,3))
end
% half-view
set(figure(12),'position',[1 1068 1432 277]); clf; 
for iTR = 1:nTR
    SP = subplot(2,9,iTR); cla; hold on; 
    errorbar(stimcond,circ_m(squeeze(matfMRI.ori.mean(iTR,:,:))'),circ_s(squeeze(matfMRI.ori.mean(iTR,:,:))')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1,'markersize',3); 
%     plot([0 180],[0 0]+circ_m(circ_m(matBeh.ori.mean')'),'k--'); 
    xlabel('orientation (deg)'); ylabel('bias (deg)'); 
    title(['TR=' num2str(iTR)]); 
    xlim([0 180]); 
    xticks(linspace(0,180,3))
end

matfMRI.ori.mean(:,1:7,:)




set(figure(12),'position',[1 1169 398 176]); clf; 
SP = subplot(1,2,1); cla; hold on; 
errorbar(refs,circ_m(matBeh.ref.mean'),circ_s(matBeh.ref.mean')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
plot([-22, 22],[0 0]+circ_m(circ_m(matBeh.ref.mean')'),'k--'); 
xlabel('reference'); ylabel('bias (deg)');
xlim([-25, 25]); 

SP = subplot(1,2,2); cla; hold on; 
errorbar(refs,circ_m(matBeh.ref.std'),circ_s(matBeh.ref.std')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
xlabel('reference'); ylabel('variability (iqr)');
xlim([-25, 25]); 


%% fMRI results
dec_mne_c1 = nan(size(dec_mne)); 
dec_mne_corr = nan(size(dec_mne)); 
for iTR = 1:nTR
    dec_mne_c1(iTR,:) = cond_correct(err_m, dec_mne(iTR,:), ref_m, stim_m, 'n', 'c'); 
    dec_mne_corr(iTR,:) = cond_correct(err_m, dec_mne(iTR,:), ref_m, stim_m, 'c', 'c'); 
end


for iTR = 1:nTR
    for ir = 1:length(refs)
        ind = find(ref_m==refs(ir) & ~isnan(err_m) & abs(err_m)<lapse_crit); 
        matfMRI.mean(ir,iTR) = circ_m(dec_mne_c1(iTR,ind)'); 
        matfMRI.iqr(ir,iTR)  = circ_s(dec_mne_c1(iTR,ind)');
%         matfMRI.iqr(ir,iTR)  = iqr(dec_mne_c1(iTR,ind)); 
        
        ind = find(ref_m==refs(ir) & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==1); 
        matfMRI.mean_e(ir,iTR) = circ_m(dec_mne_c1(iTR,ind)'); 
        matfMRI.iqr_e(ir,iTR)  = circ_s(dec_mne_c1(iTR,ind)'); 
%         matfMRI.iqr_e(ir,iTR)  = iqr(dec_mne_c1(iTR,ind)); 
        
        ind = find(ref_m==refs(ir) & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==2); 
        matfMRI.mean_l(ir,iTR) = circ_m(dec_mne_c1(iTR,ind)'); 
        matfMRI.iqr_l(ir,iTR)  = circ_s(dec_mne_c1(iTR,ind)'); 
%         matfMRI.iqr_l(ir,iTR)  = iqr(dec_mne_c1(iTR,ind)); 
        
        ind = find(ref_m==refs(ir) & ~isnan(err_m) & abs(err_m)<lapse_crit); 
        matfMRI.mean_corr(ir,iTR) = circ_m(dec_mne_corr(iTR,ind)'); 
        matfMRI.iqr_corr(ir,iTR)  = circ_s(dec_mne_corr(iTR,ind)'); 
%         matfMRI.iqr_corr(ir,iTR)  = iqr(dec_mne_corr(iTR,ind)); 
        
        ind = find(ref_m==refs(ir) & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==1); 
        matfMRI.mean_corr_e(ir,iTR) = circ_m(dec_mne_corr(iTR,ind)'); 
        matfMRI.iqr_corr_e(ir,iTR)  = circ_s(dec_mne_corr(iTR,ind)'); 
%         matfMRI.iqr_corr_e(ir,iTR)  = iqr(dec_mne_corr(iTR,ind)); 
        
        ind = find(ref_m==refs(ir) & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==2); 
        matfMRI.mean_corr_l(ir,iTR) = circ_m(dec_mne_corr(iTR,ind)'); 
        matfMRI.iqr_corr_l(ir,iTR)  = circ_s(dec_mne_corr(iTR,ind)'); 
%         matfMRI.iqr_corr_l(ir,iTR)  = iqr(dec_mne_corr(iTR,ind)); 
    end
end




%% Behavior summary figures 
set(figure(1),'position',[1 1200 398 145]); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(refs,matMerge.mean ,'ko-','markerfacecolor','w'); 
plot([-22, 22],[0 0]+nanmean(matMerge.mean),'k--'); 
xlabel('reference'); ylabel('bias (deg)');
xlim([-25, 25]); 

SP = subplot(1,2,2); cla; hold on; 
plot(refs,matMerge.iqr,'ko-','markerfacecolor','w'); 
xlabel('reference'); ylabel('variability (iqr)');
xlim([-25, 25]); 

set(figure(2),'position',[1 1200 398 145]); clf; 
SP = subplot(1,2,1); cla; hold on; 
errorbar(refs,circ_m(matMerge_sub.mean'),circ_s(matMerge_sub.mean') ,'ko-','markerfacecolor','w','capsize',1); 
plot([-22, 22],[0 0]+nanmean(matMerge.mean),'k--'); 
xlabel('reference'); ylabel('bias (deg)');
xlim([-25, 25]); 

SP = subplot(1,2,2); cla; hold on; 
errorbar(refs,circ_m(matMerge_sub.iqr'),circ_s(matMerge_sub.iqr') ,'ko-','markerfacecolor','w','capsize',1); 
xlabel('reference'); ylabel('variability (iqr)');
xlim([-25, 25]); 


% Correct for stimulus/reference bias 
set(figure(3),'position',[1 1145 190 200]); clf; 
SP = subplot(1,1,1); cla; hold on; 
plot(nanmean(matMerge_sub.iqr(2:4,:),1), nanmean(matMerge_sub.iqr([1 5],:),1),'ko'); 
plot([5 20], [5 20],'k--'); 
xlim([5 20]); 
ylim([5 20]); 
[h,p,ci,stats] = ttest(nanmean(matMerge_sub.iqr(2:4,:),1) - nanmean(matMerge_sub.iqr([1 5],:),1)); 
title(['p = ' num2str(p)]); 
xlabel('IQR (near)'); ylabel('IQR (far)'); 


set(figure(4),'position',[192 1145 190 200]); clf; 
SP = subplot(1,1,1); cla; hold on; 
plot(nanmean(matMerge_sub_refcorr.iqr(2:4,:),1), nanmean(matMerge_sub_refcorr.iqr([1 5],:),1),'ro'); 
plot([4 16], [4 16],'k--'); 
xlim([4 16]); 
ylim([4 16]); 
[h,p,ci,stats] = ttest(nanmean(matMerge_sub_refcorr.iqr(2:4,:),1) - nanmean(matMerge_sub_refcorr.iqr([1 5],:),1)); 
title(['p = ' num2str(p)]); 
xlabel('IQR (near)'); ylabel('IQR (far)'); 



set(figure(5),'position',[192 1145 190 200]); clf; 
SP = subplot(1,1,1); cla; hold on; 
errorbar(refs,circ_m(matMerge_sub.iqr'),circ_s(matMerge_sub.iqr') ,'ko-','markerfacecolor','w','capsize',1); 
xlabel('reference'); ylabel('variability (iqr)');
xlim([-25, 25]); 

set(figure(6),'position',[192 1145 190 200]); clf; 
SP = subplot(1,1,1); cla; hold on; 
errorbar(refs,circ_m(matMerge_sub_refcorr.iqr'),circ_s(matMerge_sub_refcorr.iqr') ,'ro-','markerfacecolor','w','capsize',1); 
xlabel('reference'); ylabel('variability (iqr)');
xlim([-25, 25]); 



% Check fidelity
set(figure(7),'position',[1 1177 206 168]); clf; 
SP = subplot(1,1,1); cla; hold on; 

xval = [1:nTR flip(1:nTR)]; 
temp = errsize_e; 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, J_el(1,:),'facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color',J_el(1,:),'linewidth',1.5)

temp = errsize_l; 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, J_el(2,:),'facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color',J_el(2,:),'linewidth',1.5)
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('45-|error| (deg)'); 


set(figure(8),'position',[1 972 599 373]); clf; 
SP = subplot(2,2,1); cla; hold on; 
plot(1:nTR, nanmean(matfMRI.iqr_e(2:4,:),1),'b--','color',J_el(1,:),'linewidth',1.4); 
plot(1:nTR, nanmean(matfMRI.iqr_e([1 5],:),1),'b-','color',J_el(1,:),'linewidth',1.4); 
xlabel('Time (TR)'); ylabel('IQR (deg)'); 
legend('near','far'); %ylim([60 100]); 

SP = subplot(2,2,3); cla; hold on; 
plot(1:nTR, nanmean(matfMRI.iqr_e(2:4,:),1) - nanmean(matfMRI.iqr_e([1 5],:),1),'ko-','markerfacecolor','w'); 
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('\Delta IQR (deg)'); 

SP = subplot(2,2,2); cla; hold on; 
plot(1:nTR, nanmean(matfMRI.iqr_l(2:4,:),1),'b--','color',J_el(2,:),'linewidth',1.4); 
plot(1:nTR, nanmean(matfMRI.iqr_l([1 5],:),1),'b-','color',J_el(2,:),'linewidth',1.4); 
xlabel('Time (TR)'); ylabel('IQR (deg)'); 
legend('near','far'); %ylim([60 100]); 

SP = subplot(2,2,4); cla; hold on; 
plot(1:nTR, nanmean(matfMRI.iqr_l(2:4,:),1) - nanmean(matfMRI.iqr_l([1 5],:),1),'ko-','markerfacecolor','w'); 
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('\Delta IQR (deg)'); 




set(figure(9),'position',[1 972 599 373]); clf; 
SP = subplot(2,2,1); cla; hold on; 
plot(1:nTR, nanmean(matfMRI.iqr_corr_e(2:4,:),1),'b--','color',J_el(1,:),'linewidth',1.4); 
plot(1:nTR, nanmean(matfMRI.iqr_corr_e([1 5],:),1),'b-','color',J_el(1,:),'linewidth',1.4); 
xlabel('Time (TR)'); ylabel('IQR (deg)'); 
legend('near','far'); %ylim([60 90]); 

SP = subplot(2,2,3); cla; hold on; 
plot(1:nTR, nanmean(matfMRI.iqr_corr_e(2:4,:),1) - nanmean(matfMRI.iqr_corr_e([1 5],:),1),'ko-','markerfacecolor','w'); 
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('\Delta IQR (deg)'); 

SP = subplot(2,2,2); cla; hold on; 
plot(1:nTR, nanmean(matfMRI.iqr_corr_l(2:4,:),1),'b--','color',J_el(2,:),'linewidth',1.4); 
plot(1:nTR, nanmean(matfMRI.iqr_corr_l([1 5],:),1),'b-','color',J_el(2,:),'linewidth',1.4); 
xlabel('Time (TR)'); ylabel('IQR (deg)'); 
legend('near','far'); %ylim([60 90]); 

SP = subplot(2,2,4); cla; hold on; 
plot(1:nTR, nanmean(matfMRI.iqr_corr_l(2:4,:),1) - nanmean(matfMRI.iqr_corr_l([1 5],:),1),'ko-','markerfacecolor','w'); 
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('\Delta IQR (deg)'); 

%

set(figure(10),'position',[1 972 599 373]); clf; 
SP = subplot(2,2,1); cla; hold on; 
plot(1:nTR, nanmean(nanmean(matfMRI_sub.iqr_e(:,2:4,:),2),3),'b--','color',J_el(1,:),'linewidth',1.4); 
plot(1:nTR, nanmean(nanmean(matfMRI_sub.iqr_e(:,[1 5],:),2),3),'b-','color',J_el(1,:),'linewidth',1.4); 
xlabel('Time (TR)'); ylabel('IQR (deg)'); 
legend('near','far'); %ylim([60 90]); 

SP = subplot(2,2,3); cla; hold on; 
plot(1:nTR, nanmean(nanmean(matfMRI_sub.iqr_e(:,2:4,:),2),3) - nanmean(nanmean(matfMRI_sub.iqr_e(:,[1 5],:),2),3),'ko-','markerfacecolor','w'); 
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('\Delta IQR (deg)'); 

SP = subplot(2,2,2); cla; hold on; 
plot(1:nTR, nanmean(nanmean(matfMRI_sub.iqr_l(:,2:4,:),2),3),'b--','color',J_el(2,:),'linewidth',1.4); 
plot(1:nTR, nanmean(nanmean(matfMRI_sub.iqr_l(:,[1 5],:),2),3),'b-','color',J_el(2,:),'linewidth',1.4); 
xlabel('Time (TR)'); ylabel('IQR (deg)'); 
legend('near','far'); %ylim([60 90]); 

SP = subplot(2,2,4); cla; hold on; 
plot(1:nTR, nanmean(nanmean(matfMRI_sub.iqr_l(:,2:4,:),2),3) - nanmean(nanmean(matfMRI_sub.iqr_l(:,[1 5],:),2),3),'ko-','markerfacecolor','w'); 
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('\Delta IQR (deg)'); 



set(figure(11),'position',[1 972 599 373]); clf; 
SP = subplot(2,2,1); cla; hold on; 
plot(1:nTR, nanmean(nanmean(matfMRI_sub_refcorr.iqr_e(:,2:4,:),2),3),'b--','color',J_el(1,:),'linewidth',1.4); 
plot(1:nTR, nanmean(nanmean(matfMRI_sub_refcorr.iqr_e(:,[1 5],:),2),3),'b-','color',J_el(1,:),'linewidth',1.4); 
xlabel('Time (TR)'); ylabel('IQR (deg)'); 
legend('near','far'); %ylim([60 90]); 

SP = subplot(2,2,3); cla; hold on; 
plot(1:nTR, nanmean(nanmean(matfMRI_sub_refcorr.iqr_e(:,2:4,:),2),3) - nanmean(nanmean(matfMRI_sub_refcorr.iqr_e(:,[1 5],:),2),3),'ko-','markerfacecolor','w'); 
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('\Delta IQR (deg)'); 

SP = subplot(2,2,2); cla; hold on; 
plot(1:nTR, nanmean(nanmean(matfMRI_sub_refcorr.iqr_l(:,2:4,:),2),3),'b--','color',J_el(2,:),'linewidth',1.4); 
plot(1:nTR, nanmean(nanmean(matfMRI_sub_refcorr.iqr_l(:,[1 5],:),2),3),'b-','color',J_el(2,:),'linewidth',1.4); 
xlabel('Time (TR)'); ylabel('IQR (deg)'); 
legend('near','far'); %ylim([60 90]); 

SP = subplot(2,2,4); cla; hold on; 
plot(1:nTR, nanmean(nanmean(matfMRI_sub_refcorr.iqr_l(:,2:4,:),2),3) - nanmean(nanmean(matfMRI_sub_refcorr.iqr_l(:,[1 5],:),2),3),'ko-','markerfacecolor','w'); 
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('\Delta IQR (deg)'); 
asdf


%% Figure drawing
% Avg 
figure(12); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(stimcond, matMerge_fmri.mean_avg(:,1,1),'r.-','linewidth',1.3,'markersize',12); 
plot(stimcond, matMerge_fmri.mean_avg(:,2,1),'b.-','linewidth',1.3,'markersize',12); 
plot([0 180],[0 0],'k--'); 

SP = subplot(1,2,2); cla; hold on; 
plot(stimcond, matMerge_fmri.mean_avg(:,1,2),'r.-','linewidth',1.3,'markersize',12); 
plot(stimcond, matMerge_fmri.mean_avg(:,2,2),'b.-','linewidth',1.3,'markersize',12); 
plot([0 180],[0 0],'k--'); 


% Early Dm
set(figure(100),'position',[1 976 1652 369]); clf; 
for iTR = 1:7
    SP = subplot(2,7,iTR); cla; hold on; 
    plot(stimcond, matMerge_fmri.mean(:,1,1,iTR),'r.-','linewidth',1.3,'markersize',12); 
    plot(stimcond, matMerge_fmri.mean(:,2,1,iTR),'b.-','linewidth',1.3,'markersize',12); 
    plot([0 180],[0 0],'k--'); 
    xlim([0 180]); %ylim([-12 12]); 
    xticks(linspace(0,180,5)); 
    %yticks(linspace(-10, 10,3)); 
    xlabel('Orientation'); ylabel('Error bias (deg)'); 
    title(['TR=' num2str(iTR)]); 
    
    SP = subplot(2,7,iTR+7); cla; hold on; 
    temp = matMerge_fmri.mean(:,1,1,iTR)-matMerge_fmri.mean(:,2,1,iTR); 
    temp(temp>90) = temp(temp>90) -180; 
    temp(temp<-90) = temp(temp<-90) +180; 
    plot(stimcond, temp,'k.-','linewidth',1.3,'markersize',12); 
    plot([0 180],[0 0],'k--'); 
    xlim([0 180]); 
    xticks(linspace(0,180,5)); 
    xlabel('Orientation'); ylabel('\Delta Error bias (deg)'); 
end

set(figure(101),'position',[1 533 1652 369]); clf; 
for iTR1 = 8:14
    iTR = iTR1 - 7; 
    SP = subplot(2,7,iTR); cla; hold on; 
    plot(stimcond, matMerge_fmri.mean(:,1,1,iTR1),'r.-','linewidth',1.3,'markersize',12); 
    plot(stimcond, matMerge_fmri.mean(:,2,1,iTR1),'b.-','linewidth',1.3,'markersize',12); 
    plot([0 180],[0 0],'k--'); 
    xlim([0 180]); %ylim([-12 12]); 
    xticks(linspace(0,180,5)); 
    %yticks(linspace(-10, 10,3)); 
    xlabel('Orientation'); ylabel('Error bias (deg)'); 
    title(['TR=' num2str(iTR1)]); 
    
    SP = subplot(2,7,iTR+7); cla; hold on; 
    temp = matMerge_fmri.mean(:,1,1,iTR1)-matMerge_fmri.mean(:,2,1,iTR1); 
    temp(temp>90) = temp(temp>90) -180; 
    temp(temp<-90) = temp(temp<-90) +180; 
    plot(stimcond, temp,'k.-','linewidth',1.3,'markersize',12); 
    plot([0 180],[0 0],'k--'); 
    xlim([0 180]); 
    xticks(linspace(0,180,5)); 
    xlabel('Orientation'); ylabel('\Delta Error bias (deg)'); 
end

% Late Dm
set(figure(102),'position',[1 976 1652 369]); clf; 
for iTR = 1:7
    SP = subplot(2,7,iTR); cla; hold on; 
    plot(stimcond, matMerge_fmri.mean(:,1,2,iTR),'r.-','linewidth',1.3,'markersize',12); 
    plot(stimcond, matMerge_fmri.mean(:,2,2,iTR),'b.-','linewidth',1.3,'markersize',12); 
    plot([0 180],[0 0],'k--'); 
    xlim([0 180]); %ylim([-12 12]); 
    xticks(linspace(0,180,5)); 
    %yticks(linspace(-10, 10,3)); 
    xlabel('Orientation'); ylabel('Error bias (deg)'); 
    title(['TR=' num2str(iTR)]); 
    
    SP = subplot(2,7,iTR+7); cla; hold on; 
    temp = matMerge_fmri.mean(:,1,2,iTR)-matMerge_fmri.mean(:,2,2,iTR); 
    temp(temp>90) = temp(temp>90) -180; 
    temp(temp<-90) = temp(temp<-90) +180; 
    plot(stimcond, temp,'k.-','linewidth',1.3,'markersize',12); 
    plot([0 180],[0 0],'k--'); 
    xlim([0 180]); 
    xticks(linspace(0,180,5)); 
    xlabel('Orientation'); ylabel('\Delta Error bias (deg)'); 
end

set(figure(103),'position',[1 533 1652 369]); clf; 
for iTR1 = 8:14
    iTR = iTR1 - 7; 
    SP = subplot(2,7,iTR); cla; hold on; 
    plot(stimcond, matMerge_fmri.mean(:,1,2,iTR1),'r.-','linewidth',1.3,'markersize',12); 
    plot(stimcond, matMerge_fmri.mean(:,2,2,iTR1),'b.-','linewidth',1.3,'markersize',12); 
    plot([0 180],[0 0],'k--'); 
    xlim([0 180]); %ylim([-12 12]); 
    xticks(linspace(0,180,5)); 
    %yticks(linspace(-10, 10,3)); 
    xlabel('Orientation'); ylabel('Error bias (deg)'); 
    title(['TR=' num2str(iTR1)]); 
    
    SP = subplot(2,7,iTR+7); cla; hold on; 
    temp = matMerge_fmri.mean(:,1,2,iTR1)-matMerge_fmri.mean(:,2,2,iTR1); 
    temp(temp>90) = temp(temp>90) -180; 
    temp(temp<-90) = temp(temp<-90) +180; 
    plot(stimcond, temp,'k.-','linewidth',1.3,'markersize',12); 
    plot([0 180],[0 0],'k--'); 
    xlim([0 180]); 
    xticks(linspace(0,180,5)); 
    xlabel('Orientation'); ylabel('\Delta Error bias (deg)'); 
end




% % Average across near/far
% % [merged] Calculate circular mean (or median) and std (or iqr) 
% matNearFar_fmri = {}; 
% matNearFar_fmri.mean = nan(length(stimcond), 2); 
% matNearFar_fmri.iqr  = nan(length(stimcond), 2);
% dec_mne_cor = dec_mne; 
% for istim = 1:length(stimcond)
%     % Calculate std 'after' subtracting mean term
%     for ir = 1:length(refs)
%         ind1 = find(stim_m==stimcond(istim) & ref_m==refs(ir) & ~isnan(dec_mne(iTR,:)) & abs(err_m)<lapse_crit); 
%         temp = err_m(ind1) - circ_mean(err_m(ind1)'*2*pi/180)*180/pi/2; 
%         temp(temp>90) = temp(temp>90) -180; 
%         temp(temp<-90) = temp(temp<-90) +180;  
%         err_m_cor(ind1) = temp; 
%     end
%     for ic = 1:2
%         if ic == 1  % Near
%             ref_ind = abs(ref_m)<5; 
%         else        % Far
%             ref_ind = abs(ref_m)>5; 
%         end
%         ind = find(stim_m==stimcond(istim) & ref_ind & ~isnan(err_m) & abs(err_m)<lapse_crit); 
% 
%         matNearFar_fmri.mean(istim,ic) = circ_mean(err_m(ind)'*2*pi/180)*180/pi/2; 
%         matNearFar_fmri.iqr(istim,ic)  = iqr(err_m_cor(ind)); 
%     end
% end
% 
% # 아래 몇 줄을 토대로 위를 고쳐야 함 
% for istim = 1:length(stimcond)
%     for ic = 1:2
%         ind = find(stim_m==stimcond(istim) & choice_m==ic & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==1 & abs(ref_m)<5); 
%         matMerge_fmri.mean(istim,ic,1,:) = circ_mean(dec_mne(:,ind)'*2*pi/180)*180/pi/2; 
%         for iTR = 1:14
%             matMerge_fmri.iqr(istim,ic,1,iTR)  = iqr(dec_mne(iTR,ind)); 
%         end
%         ind = find(stim_m==stimcond(istim) & choice_m==ic & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==2 & abs(ref_m)<5); 
%         matMerge_fmri.mean(istim,ic,2,:) = circ_mean(dec_mne(:,ind)'*2*pi/180)*180/pi/2; 
%         for iTR = 1:14
%             matMerge_fmri.iqr(istim,ic,2,iTR)  = iqr(dec_mne(iTR,ind)); 
%         end
%     end
% end




%% Average time-series (for insanity check)
stim_m = []; 
timing_m = []; 
ref_m = []; 
choice_m = []; 
err_m = []; 
dec_mne = []; 
lapse_crit = 30 ;
dec_val = []; 
dec_val_h1 = []; 
dec_val_h2 = []; 
for isub = runSub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_earlydelay/VC_sub-' sub_list(isub,:) '_meannorm_dec.mat'])
    
    % Behavioral errors
    errme = response - stimulus; 
    errme(errme>90) = errme(errme>90) -180; 
    errme(errme<-90) = errme(errme<-90) +180; 
    
    % fmri decoded error 
    temp = Decoded_result.est' - stimulus; 
    temp(temp>90) = temp(temp>90) -180; 
    temp(temp<-90) = temp(temp<-90) +180; 
    errme_fmri = temp ;
    
    % All merged
    stim_m = [stim_m stimulus];
    timing_m = [timing_m timing]; 
    ref_m = [ref_m ref]; 
    choice_m = [choice_m choice]; 
    err_m = [err_m errme]; 
    dec_mne = [dec_mne errme_fmri];
    dec_val = [dec_val Decoded_result.est']; 
    
    temp = Decoded_result.est_firstHalf' - stimulus; 
    temp(temp>90) = temp(temp>90) -180; 
    temp(temp<-90) = temp(temp<-90) +180; 
    dec_val_h1 = [dec_val_h1 temp]; 
    temp = Decoded_result.est_secondHalf' - stimulus; 
    temp(temp>90) = temp(temp>90) -180; 
    temp(temp<-90) = temp(temp<-90) +180; 
    dec_val_h2 = [dec_val_h2 temp]; 
end


for istim = 1:length(stimcond)
    
    ind = find(stim_m==stimcond(istim) & ~isnan(err_m) & abs(err_m)<lapse_crit); 
    matMerge_behav.mean(istim) = circ_mean(err_m(ind)'*2*pi/180)*180/pi/2; 
    matMerge_behav.iqr(istim) = iqr(err_m(ind)); 
    
    matMerge.mean(istim) = circ_mean(dec_mne(ind)'*2*pi/180)*180/pi/2; 
    matMerge.iqr(istim) = iqr(dec_mne(ind)); 
    
    matMerge.mean_1h(istim) = circ_mean(dec_val_h1(ind)'*2*pi/180)*180/pi/2; 
    matMerge.iqr_1h(istim) = iqr(dec_val_h1(ind)); 
    matMerge.mean_2h(istim) = circ_mean(dec_val_h2(ind)'*2*pi/180)*180/pi/2; 
    matMerge.iqr_2h(istim) = iqr(dec_val_h2(ind)); 
    for ic = 1:2
        ind = find(stim_m==stimcond(istim) & choice_m==ic & ~isnan(err_m) & abs(err_m)<lapse_crit); 

        matMerge.mean_c(istim,ic) = circ_mean(dec_mne(ind)'*2*pi/180)*180/pi/2; 
        matMerge.iqr_c(istim,ic)  = iqr(dec_mne(ind)); 
    end
end

% Figure #1: cardinal bias and variance in BEHAVIOR
set(figure(1000),'position',[1 1029 323 316]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matMerge_behav.mean - circ_mean(matMerge_behav.mean'*2*pi/180)*180/pi/2,'k-','linewidth',1.5); 
plot([0 180],[0 0],'k--'); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('error bias (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matMerge_behav.iqr,'k-','linewidth',1.5); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('error iqr (deg)'); 


% Figure #2: decoded vs true orientation (mnemonic code) 
set(figure(1001),'position',[1 636 340 333]); clf; 
SP = subplot(1,1,1); cla; hold on; 
scatter1 = scatter(stim_m, dec_val,'ko','markerfacecolor','k'); 
scatter1.MarkerFaceAlpha = 0.05;
scatter1.MarkerEdgeAlpha = 0;
xlim([0 180]); 
ylim([0 180]); 
xticks(linspace(0,180,5)); 
yticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('Decoded orientation'); 


% Figure #3: check existence of cardinal bias and variance in decoded orientation
set(figure(1002),'position',[325 1029 323 316]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matMerge.mean - circ_mean(matMerge.mean'*2*pi/180)*180/pi/2,'r-','linewidth',1.5); 
plot([0 180],[0 0],'k--'); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('error bias (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matMerge.iqr,'r-','linewidth',1.5); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('error iqr (deg)'); 

% Figure #4: compare bias 
set(figure(1003),'position',[325 1029 323 316]); clf; 
SP = subplot(2,1,1); cla; hold on; 
plot(stimcond, matMerge.mean_1h - circ_mean(matMerge.mean_1h'*2*pi/180)*180/pi/2,'b-','linewidth',1.5); 
plot(stimcond, matMerge.mean_2h - circ_mean(matMerge.mean_2h'*2*pi/180)*180/pi/2,'r-','linewidth',1.5); 
% legend('5-9TR','10-13TR');
plot([0 180],[0 0],'k--'); 
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('error bias (deg)'); 

SP = subplot(2,1,2); cla; hold on; 
plot(stimcond, matMerge.iqr_1h,'b-','linewidth',1.5); 
plot(stimcond, matMerge.iqr_2h,'r-','linewidth',1.5); 
legend('5-9TR','10-13TR');
xlim([0 180]); 
xticks(linspace(0,180,5)); 
xlabel('Orientation'); ylabel('error iqr (deg)'); 
