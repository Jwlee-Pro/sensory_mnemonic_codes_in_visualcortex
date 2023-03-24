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
dec_mne = []; 
lapse_crit = 30 ;

for isub = runSub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_earlydelay/VC_sub-' sub_list(isub,:) '_dec.mat'])
    
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
    
    % fmri decoded error 
    errme_fmri = nan(size(Decoded_result.est)); 
    for iTR = 1:14
        temp = Decoded_result.est(iTR,:) - stimulus; 
        temp(temp>90) = temp(temp>90) -180; 
        temp(temp<-90) = temp(temp<-90) +180; 
        errme_fmri(iTR,:) = temp ;
    end
    
    % All merged
    stim_m = [stim_m stimulus];
    timing_m = [timing_m timing]; 
    ref_m = [ref_m ref]; 
    choice_m = [choice_m choice]; 
    err_m = [err_m errme]; 
    dec_mne = [dec_mne errme_fmri];
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



% fmri analysis
matMerge_fmri = {}; 
matMerge_fmri.mean = nan(length(stimcond), 2, 2, 14); 
matMerge_fmri.iqr  = nan(length(stimcond), 2, 2, 14);

% for istim = 1:length(stimcond)
%     for ic = 1:2
%         ind = find(stim_m==stimcond(istim) & choice_m==ic & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==1); 
%         matMerge_fmri.mean(istim,ic,1,:) = circ_mean(dec_mne(:,ind)'*2*pi/180)*180/pi/2; 
%         for iTR = 1:14
%             matMerge_fmri.iqr(istim,ic,1,iTR)  = iqr(dec_mne(iTR,ind)); 
%         end
%         ind = find(stim_m==stimcond(istim) & choice_m==ic & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==2); 
%         matMerge_fmri.mean(istim,ic,2,:) = circ_mean(dec_mne(:,ind)'*2*pi/180)*180/pi/2; 
%         for iTR = 1:14
%             matMerge_fmri.iqr(istim,ic,2,iTR)  = iqr(dec_mne(iTR,ind)); 
%         end
%     end
% end

% consider only near condition 
for istim = 1:length(stimcond)
    for ic = 1:2
        ind = find(stim_m==stimcond(istim) & choice_m==ic & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==1 & abs(ref_m)<5); 
        matMerge_fmri.mean(istim,ic,1,:) = circ_mean(dec_mne(:,ind)'*2*pi/180)*180/pi/2; 
        for iTR = 1:14
            matMerge_fmri.iqr(istim,ic,1,iTR)  = iqr(dec_mne(iTR,ind)); 
        end
        ind = find(stim_m==stimcond(istim) & choice_m==ic & ~isnan(err_m) & abs(err_m)<lapse_crit & timing_m==2 & abs(ref_m)<5); 
        matMerge_fmri.mean(istim,ic,2,:) = circ_mean(dec_mne(:,ind)'*2*pi/180)*180/pi/2; 
        for iTR = 1:14
            matMerge_fmri.iqr(istim,ic,2,iTR)  = iqr(dec_mne(iTR,ind)); 
        end
    end
end


%% Figure drawing
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


% Check fidelity
figure(); 
hold on; 
plot(45-circ_mean(abs(dec_mne(:,timing_m==1))'*2*pi/180)*180/pi/2, 'bo-','markerfacecolor','w'); 
plot(45-circ_mean(abs(dec_mne(:,timing_m==2))'*2*pi/180)*180/pi/2, 'ro-','markerfacecolor','w'); 
legend('Early','Late');
plot([0 14],[0 0],'k--'); 
xlabel('Time (TR)'); ylabel('45-|err|'); 



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














