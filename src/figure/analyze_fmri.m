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
            matfMRI.ref.mean(iTR,iref,isub) = circ_m(errme_fmri(iTR,ind)'); 
            matfMRI.ref.std(iTR,iref,isub)  = circ_s(errme_fmri(iTR,ind)'); 
        end
        
        % Early Dm trials
        ind = find(ref==refs(iref) & ~isnan(response) & abs(errme)<lapse_crit & timing==1); 
        for iTR = 1:nTR
            matfMRI.ref.mean_e(iTR,iref,isub) = circ_m(errme_fmri(iTR,ind)'); 
            matfMRI.ref.std_e(iTR,iref,isub)  = circ_s(errme_fmri(iTR,ind)'); 
        end
        
        % Late Dm trials
        ind = find(ref==refs(iref) & ~isnan(response) & abs(errme)<lapse_crit & timing==2); 
        for iTR = 1:nTR
            matfMRI.ref.mean_l(iTR,iref,isub) = circ_m(errme_fmri(iTR,ind)'); 
            matfMRI.ref.std_l(iTR,iref,isub)  = circ_s(errme_fmri(iTR,ind)'); 
        end
    end
    
    % Choice-dependent features
    for ic = 1:2
        for iTR = 1:nTR
            % All trials
            ind = find(choice==ic & ~isnan(response) & abs(errme)<lapse_crit & abs(ref)<5); 
            matfMRI.choice.mean(iTR,ic,isub) = circ_m(errme_fmri(iTR,ind)'); 
            
            % Early Dm trials
            ind = find(choice==ic & ~isnan(response) & abs(errme)<lapse_crit & abs(ref)<5 & timing==1); 
            matfMRI.choice.mean_e(iTR,ic,isub) = circ_m(errme_fmri(iTR,ind)'); 
            
            % Late Dm trials
            ind = find(choice==ic & ~isnan(response) & abs(errme)<lapse_crit & abs(ref)<5 & timing==2); 
            matfMRI.choice.mean_l(iTR,ic,isub) = circ_m(errme_fmri(iTR,ind)'); 
        end
    end
    
    % overall variance level
    ind = find(~isnan(response) & abs(errme)<lapse_crit); 
    for iTR = 1:nTR
        std_trend(iTR,isub)  = circ_s(errme_fmri(iTR,ind)'); 
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
xlabel('orientation (deg)'); ylabel('variability (deg)');
xlim([0 180]); 


set(figure(2),'position',[1 1169 398 176]); clf; 
SP = subplot(1,2,1); cla; hold on; 
errorbar(refs,circ_m(matBeh.ref.mean'),circ_s(matBeh.ref.mean')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
plot([-22, 22],[0 0]+circ_m(circ_m(matBeh.ref.mean')'),'k--'); 
xlabel('reference'); ylabel('bias (deg)');
xlim([-25, 25]); 

SP = subplot(1,2,2); cla; hold on; 
errorbar(refs,circ_m(matBeh.ref.std'),circ_s(matBeh.ref.std')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
xlabel('reference'); ylabel('variability (deg)');
xlim([-25, 25]); 

%% Figure: orientation-specific bias/variance (fMRI)
set(figure(11),'position',[1 1068 1432 277]); clf; 
for iTR = 1:nTR
    SP = subplot(2,9,iTR); cla; hold on; 
    errorbar(stimcond,circ_m(squeeze(matfMRI.ori.mean(iTR,:,:))'),circ_s(squeeze(matfMRI.ori.mean(iTR,:,:))')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1,'markersize',3); 
    plot([0 180],[0 0]+circ_m(circ_m(squeeze(matfMRI.ori.mean(iTR,:,:))')'),'k--'); 
    xlabel('orientation (deg)'); ylabel('bias (deg)'); 
    title(['TR=' num2str(iTR)]); 
    xlim([0 180]); 
    xticks(linspace(0,180,3))
end
set(figure(12),'position',[1 1068 1432 277]); clf; 
for iTR = 1:nTR
    SP = subplot(2,9,iTR); cla; hold on; 
    errorbar(stimcond,circ_m(squeeze(matfMRI.ori.std(iTR,:,:))'),circ_s(squeeze(matfMRI.ori.std(iTR,:,:))')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1,'markersize',3); 
    xlabel('orientation (deg)'); ylabel('variability (deg)'); 
    title(['TR=' num2str(iTR)]); 
    xlim([0 180]); 
    xticks(linspace(0,180,3))
end

% Align to cardinal orientation (0 & 90)
stimHalf = [stimcond(19:24)-180 stimcond(1:7)]; 
set(figure(3),'position',[1 1169 398 176]); clf; 
SP = subplot(1,2,1); cla; hold on; 
errorbar(stimHalf,circ_m(ori_cardinal(matBeh.ori.mean)'),circ_s(ori_cardinal(matBeh.ori.mean)')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
plot([-45 45],[0 0]+circ_m(circ_m(ori_cardinal(matBeh.ori.mean)')'),'k--'); 
xlabel('orientation from cardinal (deg)'); ylabel('bias (deg)');
xlim([-45 45]); xticks(linspace(-45,45,3));

SP = subplot(1,2,2); cla; hold on; 
errorbar(stimHalf,circ_m(ori_cardinal(matBeh.ori.std)'),circ_s(ori_cardinal(matBeh.ori.std)')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1); 
xlabel('orientation from cardinal (deg)'); ylabel('variability (deg)');
xlim([-45 45]); xticks(linspace(-45,45,3));

set(figure(12),'position',[1 1068 1432 277]); clf; 
for iTR = 1:nTR
    SP = subplot(2,9,iTR); cla; hold on; 
    temp = ori_cardinal(squeeze(matfMRI.ori.mean(iTR,:,:)));
    errorbar(stimHalf,circ_m(temp'),circ_s(temp')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1,'markersize',3); 
    plot([-45 45],[0 0]+circ_m(circ_m(temp')'),'k--'); 
    xlabel('orientation (deg)'); ylabel('bias (deg)'); 
    title(['TR=' num2str(iTR)]); 
    xlim([-45 45]); 
    xticks(linspace(-45,45,3));
end

set(figure(13),'position',[1 1068 1432 277]); clf; 
for iTR = 1:nTR
    SP = subplot(2,9,iTR); cla; hold on; 
    temp = ori_cardinal(squeeze(matfMRI.ori.std(iTR,:,:)));
    errorbar(stimHalf,circ_m(temp'),circ_s(temp')/sqrt(length(sub_list)-1) ,'ko-','markerfacecolor','w','capsize',1,'markersize',3); 
    xlabel('orientation (deg)'); ylabel('variability (deg)'); 
    title(['TR=' num2str(iTR)]); 
    xlim([-45 45]); 
    xticks(linspace(-45,45,3));
end



set(figure(14),'position',[1 379 398 176]); clf; 
% early stage (TR5~9) vs. late stage (TR10~13)
temp1 = matfMRI.ori.mean(5:13,:,:);
temp2 = matfMRI.ori.std(5:13,:,:);
for isub = 1:length(sub_list)
    e1(:,isub) = circ_m(temp1(:,:,isub)); 
    k1(:,isub) = circ_m(temp2(:,:,isub)); 
end

SP = subplot(1,2,1); cla; hold on; 
errorbar(stimHalf,circ_m(ori_cardinal(e1)'),circ_s(ori_cardinal(e1)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1); 
% plot([-45 45],[0 0]+circ_m(circ_m(ori_cardinal(e1)')'),'k--'); 
xlabel('orientation from cardinal (deg)'); ylabel('bias (deg)');
xlim([-45 45]); xticks(linspace(-45,45,3));
plot([-45 45],[0 0]+circ_m(circ_m(ori_cardinal(e1)')'),'k--'); 
SP = subplot(1,2,2); cla; hold on; 
errorbar(stimHalf,circ_m(ori_cardinal(k1)'),circ_s(ori_cardinal(k1)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1); 
% plot([-45 45],[0 0]+circ_m(circ_m(ori_cardinal(e1)')'),'k--'); 
xlabel('orientation from cardinal (deg)'); ylabel('variability (deg)');
xlim([-45 45]); xticks(linspace(-45,45,3));





set(figure(14),'position',[1 379 398 176]); clf; 
% early stage (TR5~9) vs. late stage (TR10~13)
temp1 = matfMRI.ori.mean(5:9,:,:);
temp2 = matfMRI.ori.mean(10:13,:,:);
for isub = 1:length(sub_list)
    e1(:,isub) = circ_m(temp1(:,:,isub)); 
    e2(:,isub) = circ_m(temp2(:,:,isub)); 
end

SP = subplot(1,2,1); cla; hold on; 
errorbar(stimHalf,circ_m(ori_cardinal(e1)'),circ_s(ori_cardinal(e1)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1); 
errorbar(stimHalf,circ_m(ori_cardinal(e2)'),circ_s(ori_cardinal(e2)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1,'color',[0 0 0]+0.7); 
% plot([-45 45],[0 0]+circ_m(circ_m(ori_cardinal(e1)')'),'k--'); 
xlabel('orientation from cardinal (deg)'); ylabel('bias (deg)');
xlim([-45 45]); xticks(linspace(-45,45,3));
legend('5-9TR','10-13TR');

% early stage (TR5~9) vs. late stage (TR10~13)
temp1 = matfMRI.ori.std(5:9,:,:);
temp2 = matfMRI.ori.std(10:13,:,:);
for isub = 1:length(sub_list)
    k1(:,isub) = circ_m(temp1(:,:,isub)); 
    k2(:,isub) = circ_m(temp2(:,:,isub)); 
end

SP = subplot(1,2,2); cla; hold on; 
errorbar(stimHalf,circ_m(ori_cardinal(k1)'),circ_s(ori_cardinal(k1)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1); 
errorbar(stimHalf,circ_m(ori_cardinal(k2)'),circ_s(ori_cardinal(k2)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1,'color',[0 0 0]+0.7); 
% plot([-45 45],[0 0]+circ_m(circ_m(ori_cardinal(e1)')'),'k--'); 
xlabel('orientation from cardinal (deg)'); ylabel('variability (deg)');
xlim([-45 45]); xticks(linspace(-45,45,3));




set(figure(15),'position',[1 129 398 176]); clf; 
SP = subplot(1,2,1); cla; hold on; 
pre = ori_cardinal(e1); 
post= ori_cardinal(e2); 
% plot(circ_m([pre(8:9,:); -pre(5:6,:)]), circ_m([post(8:9,:); -post(5:6,:)]),'ko');
plot(nanmean([pre(8:end,:); -pre(1:6,:)],1), nanmean([post(8:end,:); -post(1:6,:)],1),'ko');
[h,p,ci,stats] = ttest(nanmean([pre(8:end,:); -pre(1:6,:)],1)- nanmean([post(8:end,:); -post(1:6,:)],1))
title(['p=' num2str(p)]);
plot([-30 30], [-30 30],'k--'); 
xlim([-30 30]); ylim([-30 30]); 
xlabel('cardinal bias 5~9TR'); ylabel('cardinal bias 10~13TR');

SP = subplot(1,2,2); cla; hold on; 
plot(nanmean(std_trend(5:9,:),1), nanmean(std_trend(10:13,:),1),'ko');
[h,p,ci,stats] = ttest(nanmean(std_trend(5:9,:),1) - nanmean(std_trend(10:13,:),1))
title(['p=' num2str(p)]);
plot([20 40], [20 40],'k--'); 
xlim([20 40]); ylim([20 40]); 
xlabel('variability 5~9TR'); ylabel('variability 10~13TR');




% Choice-dependent bias
pte = nan(1,nTR); 
ptl = nan(1,nTR); 
for iTR = 1:nTR
    [h,pte(iTR),ci,stats] = ttest(squeeze(matfMRI.choice.mean_e(iTR,1,:)) - squeeze(matfMRI.choice.mean_e(iTR,2,:)));
    [h,ptl(iTR),ci,stats] = ttest(squeeze(matfMRI.choice.mean_l(iTR,1,:)) - squeeze(matfMRI.choice.mean_l(iTR,2,:)));
end

set(figure(16),'position',[1 632 457 173]); clf; 
SP = subplot(1,2,1); cla; hold on; 
xval = [1:nTR flip(1:nTR)]; 
temp = squeeze(matfMRI.choice.mean_e(:,1,:)); 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, 'r','facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color','r','linewidth',1.5)

temp = squeeze(matfMRI.choice.mean_e(:,2,:)); 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, 'b','facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color','b','linewidth',1.5)
plot([1 14],[0 0],'k--'); 
ylim([-20 20]); 
xlabel('Time (TR)'); ylabel('bias (deg)');
plot(find(pte<0.05), 18*ones(1,sum(pte<0.05)),'k*'); 

SP = subplot(1,2,2); cla; hold on; 
xval = [1:nTR flip(1:nTR)]; 
temp = squeeze(matfMRI.choice.mean_l(:,1,:)); 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, 'r','facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color','r','linewidth',1.5)

temp = squeeze(matfMRI.choice.mean_l(:,2,:)); 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, 'b','facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color','b','linewidth',1.5)
plot([1 14],[0 0],'k--'); 
ylim([-20 20]); 
xlabel('Time (TR)'); ylabel('bias (deg)');
plot(find(pte<0.05), 18*ones(1,sum(ptl<0.05)),'k*'); 



% Near vs. Far condition comparisons 
for isub = 1:length(sub_list)
    temp = matfMRI.ref.std_e(:,:,isub); 
    near_e(:,isub) = circ_m(temp(:,2:4)'); 
    far_e(:,isub) = circ_m(temp(:,[1 5])'); 
    
    temp = matfMRI.ref.std_l(:,:,isub); 
    near_l(:,isub) = circ_m(temp(:,2:4)'); 
    far_l(:,isub) = circ_m(temp(:,[1 5])'); 
end
pte = nan(1,nTR); 
ptl = nan(1,nTR); 
for iTR = 1:nTR
    [h,pte(iTR),ci,stats] = ttest(near_e(iTR,:) - far_e(iTR,:));
    [h,ptl(iTR),ci,stats] = ttest(near_l(iTR,:) - far_l(iTR,:));
end

set(figure(17),'position',[459 632 457 173]); clf; 
SP = subplot(1,2,1); cla; hold on; 
xval = [1:nTR flip(1:nTR)]; 
temp = near_e; 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, nearfar_color(1,:),'facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color',nearfar_color(1,:),'linewidth',1.5)

temp = far_e; 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, nearfar_color(2,:),'facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color',nearfar_color(2,:),'linewidth',1.5)

xlabel('Time (TR)'); ylabel('variablity (deg)'); 
plot(find(pte<0.05), 36*ones(1,sum(pte<0.05)),'k*'); 

SP = subplot(1,2,2); cla; hold on; 
temp = near_l; 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, nearfar_color(1,:),'facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color',nearfar_color(1,:),'linewidth',1.5)

temp = far_l; 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval+meanval +flip(stdval+meanval)] ; 
patch(xval, yval, nearfar_color(2,:),'facealpha',0.2,'edgecolor','none'); 
plot(1:nTR,meanval,'color',nearfar_color(2,:),'linewidth',1.5)

xlabel('Time (TR)'); ylabel('variablity (deg)');
plot(find(ptl<0.05), 36*ones(1,sum(ptl<0.05)),'k*'); 


% Merging TRs (6~11TR) in early Dm to see overall effect
set(figure(18),'position',[459 351 457 204]); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(nanmean(near_e(6:11,:),1), nanmean(far_e(6:11,:),1),'ko') 
plot([24 40],[24 40],'k--');
[h,p,ci,stats] = ttest(nanmean(near_e(6:11,:),1) - nanmean(far_e(6:11,:),1))
xlabel('Near variability'); ylabel('Far variability');
xlim([24 40]); ylim([24 40]); 
title(['p=' num2str(p) ]); 

SP = subplot(1,2,2); cla; hold on; 
plot(nanmean(near_l(9:11,:),1), nanmean(far_l(9:11,:),1),'ko') 
[h,p,ci,stats] = ttest(nanmean(near_l(9:11,:),1) - nanmean(far_l(9:11,:),1))
plot([24 40],[24 40],'k--');
xlabel('Near variability'); ylabel('Far variability');
xlim([24 40]); ylim([24 40]); 
title(['p=' num2str(p) ]); 

