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

%% P2-(1) 
% Sensory information during encoding, Dm, reproduction periods decoded with sensory-trained coding scheme.
% Hypothesis (Decoded_#_current_previous)
clear ModelName
ModelName{1} = 'Decoded_1_None_Stim'; 
ModelName{2} = 'Decoded_2_None_Esti'; 
ModelName{3} = 'Decoded_3_Stim_None'; 
ModelName{4} = 'Decoded_4_Stim_Stim'; 
ModelName{5} = 'Decoded_5_Stim_Esti'; 
ModelName{6} = 'Decoded_6_Esti_None'; 
ModelName{7} = 'Decoded_7_Esti_Stim'; 
ModelName{8} = 'Decoded_8_Esti_Esti'; 

runSub = 1:50; 

iModels = 3; 

nTR = 14; 
refs = [-21, -4, 0, 4, 21]; 
stimcond = 0:7.5:172.5; 

binCenter = -90:5:85 ;
nBin = length(binCenter); 

ROIs = {'V1','V2','V3','V4','V3ab','LO1','IPS0','IPS1'}; 


iROI = 3;



shiftCH_e = nan(nTR,120,length(sub_list)); 
shiftCH_l = nan(nTR,120,length(sub_list)); 
shiftCHest_e = nan(nTR,120,length(sub_list)); 
shiftCHest_l = nan(nTR,120,length(sub_list)); 
shiftCH_ref_e = nan(nTR,5,120,length(sub_list)); 
shiftCH_ref_l = nan(nTR,5,120,length(sub_list)); 


decoded_merge = []; 
decErr_merge = []; 
stimulus_merge = []; 
timing_merge = []; 

corrval = nan(nTR,length(sub_list)); 
corrval_e = nan(nTR,length(sub_list)); 
corrval_l = nan(nTR,length(sub_list)); 
corrStartDeg = nan(nTR,length(sub_list)); 
corrStartDeg_e = nan(nTR,length(sub_list)); 
corrStartDeg_l = nan(nTR,length(sub_list)); 
corrEstiDeg = nan(nTR,length(sub_list)); 
corrEstiDeg_e = nan(nTR,length(sub_list)); 
corrEstiDeg_l = nan(nTR,length(sub_list)); 

fidelity_decx = cell(1,length(sub_list)); 
fidelityest_decx = cell(1,length(sub_list)); 

aDecAngle_sen_e = nan(nTR, length(stimcond), length(sub_list)); 
aDecAngle_sen_l = nan(nTR, length(stimcond), length(sub_list)); 
aDecAngle_mne_e = nan(nTR, length(stimcond), length(sub_list)); 
aDecAngle_mne_l = nan(nTR, length(stimcond), length(sub_list)); 
aDecAngle_TR_e = nan(nTR, length(stimcond), length(sub_list)); 
aDecAngle_TR_l = nan(nTR, length(stimcond), length(sub_list)); 

aDecErr_sen_e = nan(nTR, 2, length(sub_list)); 
aDecErr_sen_l = nan(nTR, 2, length(sub_list)); 
aDecErr_mne_e = nan(nTR, 2, length(sub_list)); 
aDecErr_mne_l = nan(nTR, 2, length(sub_list)); 
aDecErr_TR_e  = nan(nTR, 2, length(sub_list)); 
aDecErr_TR_l  = nan(nTR, 2, length(sub_list)); 

temp_choice_sen_e = nan(nTR,2,length(sub_list));
temp_choice_sen_l = nan(nTR,2,length(sub_list));
temp_choice_mne_e = nan(nTR,2,length(sub_list));
temp_choice_mne_l = nan(nTR,2,length(sub_list));
tempx_choice_sen_e = nan(nTR,2,length(sub_list));
tempx_choice_sen_l = nan(nTR,2,length(sub_list));
tempx_choice_mne_e = nan(nTR,2,length(sub_list));
tempx_choice_mne_l = nan(nTR,2,length(sub_list));

temp_ref_sen_e = nan(nTR,5,length(sub_list));
temp_ref_sen_l = nan(nTR,5,length(sub_list));
temp_ref_mne_e = nan(nTR,5,length(sub_list));
temp_ref_mne_l = nan(nTR,5,length(sub_list));
std_ref_sen_e = nan(nTR,5,length(sub_list));
std_ref_sen_l = nan(nTR,5,length(sub_list));
std_ref_mne_e = nan(nTR,5,length(sub_list));
std_ref_mne_l = nan(nTR,5,length(sub_list));

Behav_std_ref = nan(5,length(sub_list));
Behav_mean_ref = nan(5,length(sub_list));


timing_m = []; 
ref_m = []; 
choice_m = []; 
dec_tr = []; 
dec_sen = []; 
dec_mne = []; 
for isub = runSub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_estimated/VC_sub-' sub_list(isub,:) '_dec.mat'])
    [nTestTR, nTrials] = size(Decoded_result{1}.est); 
    
    % Behavioral errors
    errme = response - stimulus; 
    errme(errme>90) = errme(errme>90) -180; 
    errme(errme<-90) = errme(errme<-90) +180; 
    for ir = 1:length(refs)
        Behav_mean_ref(ir,isub) = circ_mean((errme(ref==refs(ir) & ref==refs(ir) & ~isnan(errme)))'*2*pi/180)*180/pi/2;
        Behav_std_ref(ir,isub) = circ_std((errme(ref==refs(ir) & ref==refs(ir) & ~isnan(errme)))'*2*pi/180)*180/pi/2;
    end
    
    % Sensory decoded
    errme_sensory = nan(length(Decoded_result{3}.est(1,:)), nTR); 
    for iTRx = 1:nTR
        tempx = circ_mean([Decoded_result{3}.est(iTRx,:); Decoded_result{4}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        errme = tempx - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        errme_sensory(:,iTRx) = errme; 
    end
    temp_sen_mean_e(:,isub) = 45-circ_mean(abs(errme_sensory(timing==1,:))*2*pi/180)*180/pi/2;
    temp_sen_mean_l(:,isub) = 45-circ_mean(abs(errme_sensory(timing==2,:))*2*pi/180)*180/pi/2;
    temp_choice_sen_e(:,1,isub) = circ_mean((errme_sensory(timing==1 & choice==1 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    temp_choice_sen_e(:,2,isub) = circ_mean((errme_sensory(timing==1 & choice==2 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    temp_choice_sen_l(:,1,isub) = circ_mean((errme_sensory(timing==2 & choice==1 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    temp_choice_sen_l(:,2,isub) = circ_mean((errme_sensory(timing==2 & choice==2 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    
    for ir = 1:5
        temp_ref_sen_e(:,ir,isub) = circ_mean((errme_sensory(timing==1 & ref==refs(ir),:))*2*pi/180)*180/pi/2;
        std_ref_sen_e(:,ir,isub)  = circ_std((errme_sensory(timing==1 & ref==refs(ir),:))*2*pi/180)*180/pi/2;
        temp_ref_sen_l(:,ir,isub) = circ_mean((errme_sensory(timing==2 & ref==refs(ir),:))*2*pi/180)*180/pi/2;
        std_ref_sen_l(:,ir,isub)  = circ_std((errme_sensory(timing==2 & ref==refs(ir),:))*2*pi/180)*180/pi/2;
    end
    
    % Memory decoded
    errme_memory = nan(length(Decoded_result{3}.est(1,:)), nTR); 
    for iTRx = 1:nTR
        tempx = circ_mean([Decoded_result{9}.est(iTRx,:); Decoded_result{10}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        errme = tempx - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        errme_memory(timing==1,iTRx) = errme(timing==1); 
        
        tempx = circ_mean([Decoded_result{6}.est(iTRx,:); Decoded_result{7}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        errme = tempx - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        errme_memory(timing==2,iTRx) = errme(timing==2); 
    end
    
    for ir = 1:5
        temp_ref_mne_e(:,ir,isub) = circ_mean((errme_memory(timing==1 & ref==refs(ir),:))*2*pi/180)*180/pi/2;
        std_ref_mne_e(:,ir,isub)  = circ_std((errme_memory(timing==1 & ref==refs(ir),:))*2*pi/180)*180/pi/2;
        temp_ref_mne_l(:,ir,isub) = circ_mean((errme_memory(timing==2 & ref==refs(ir),:))*2*pi/180)*180/pi/2;
        std_ref_mne_l(:,ir,isub)  = circ_std((errme_memory(timing==2 & ref==refs(ir),:))*2*pi/180)*180/pi/2;
    end
    
    % choice-dependent bias
%     temp_mean_e(:,isub) = 45-circ_mean(errme_memory(timing==1,:)*2*pi/180)*180/pi/2;
    temp_mne_mean_e(:,isub) = 45-circ_mean(abs(errme_memory(timing==1,:))*2*pi/180)*180/pi/2;
    temp_mne_mean_l(:,isub) = 45-circ_mean(abs(errme_memory(timing==2,:))*2*pi/180)*180/pi/2;
    
    temp_choice_mne_e(:,1,isub) = circ_mean((errme_memory(timing==1 & choice==1 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    temp_choice_mne_e(:,2,isub) = circ_mean((errme_memory(timing==1 & choice==2 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    temp_choice_mne_l(:,1,isub) = circ_mean((errme_memory(timing==2 & choice==1 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    temp_choice_mne_l(:,2,isub) = circ_mean((errme_memory(timing==2 & choice==2 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    
    
    % TR-TR
    errme_TR = nan(length(Decoded_result{3}.est(1,:)), nTR); 
    for iTRx = 1:nTR
        errme = Decoded_result{iTRx}.est(iTRx,:) - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        errme_TR(:,iTRx) = errme; 
    end
    
    % All merged
    timing_m = [timing_m timing]; 
    ref_m = [ref_m ref]; 
    choice_m = [choice_m choice]; 
    dec_tr = [dec_tr; errme_TR]; 
    dec_sen = [dec_sen; errme_sensory]; 
    dec_mne = [dec_mne; errme_memory]; 
    
    
    
    % choice-dependent bias
    temp_TR_mean_e(:,isub) = 45-circ_mean(abs(errme_TR(timing==1,:))*2*pi/180)*180/pi/2;
    temp_TR_mean_l(:,isub) = 45-circ_mean(abs(errme_TR(timing==2,:))*2*pi/180)*180/pi/2;
    
    temp_choice_TR_e(:,1,isub) = circ_mean((errme_TR(timing==1 & choice==1 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    temp_choice_TR_e(:,2,isub) = circ_mean((errme_TR(timing==1 & choice==2 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    temp_choice_TR_l(:,1,isub) = circ_mean((errme_TR(timing==2 & choice==1 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    temp_choice_TR_l(:,2,isub) = circ_mean((errme_TR(timing==2 & choice==2 & abs(ref)<5,:))*2*pi/180)*180/pi/2;
    
    errme_TR = nan(length(Decoded_result{3}.est(1,:)), nTR); 
    for iTRx = 1:nTR
        errme = Decoded_result{iTRx}.est(iTRx,:) - double(stimulus) - double(ref); 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        errme_TR(:,iTRx) = errme; 
        
        avgPred_e(iTRx,isub) = (sum(errme(timing==1)>0 & choice(timing==1)==1) + sum(errme(timing==1)<0 & choice(timing==1)==2))/length(choice(timing==1)); 
        avgPred_l(iTRx,isub) = (sum(errme(timing==2)>0 & choice(timing==2)==1) + sum(errme(timing==2)<0 & choice(timing==2)==2))/length(choice(timing==2)); 
    end
    
    
    dec_dv = nan(length(Decoded_result{3}.est(1,:)), nTR); 
    for iTRx = 1:nTR
        tempse = circ_mean([Decoded_result{3}.est(iTRx,:); Decoded_result{4}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempse(tempse<0) = tempse(tempse<0)+180; 
        
        tempx = circ_mean([Decoded_result{9}.est(iTRx,:); Decoded_result{10}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        
        errme = tempx - tempse; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        dec_dv(timing==1,iTRx) = errme(timing==1); 

        tempx = circ_mean([Decoded_result{6}.est(iTRx,:); Decoded_result{7}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        
        errme = tempx - tempse; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        dec_dv(timing==2,iTRx) = errme(timing==2); 
%         avgPred_Near_e(iTRx,isub) = (sum(dec_dv(timing==1)>0 & choice(timing==1)==1 & abs(ref)<5) + sum(dec_dv(timing==1)<0 & choice(timing==1)==2 & abs(ref)<5))/length(choice(timing==1 & abs(ref)<5)); 
%         avgPred_Near_l(iTRx,isub) = (sum(dec_dv(timing==2)>0 & choice(timing==2)==1 & abs(ref)<5) + sum(dec_dv(timing==2)<0 & choice(timing==2)==2 & abs(ref)<5))/length(choice(timing==2 & abs(ref)<5)); 
%         avgPred_All_e(iTRx,isub) = (sum(dec_dv(timing==1)>0 & choice(timing==1)==1 & abs(ref)<5) + sum(dec_dv(timing==1)<0 & choice(timing==1)==2 & abs(ref)<5))/length(choice(timing==1 & abs(ref)<5)); 
%         avgPred_All_l(iTRx,isub) = (sum(dec_dv(timing==2)>0 & choice(timing==2)==1 & abs(ref)<5) + sum(dec_dv(timing==2)<0 & choice(timing==2)==2 & abs(ref)<5))/length(choice(timing==2 & abs(ref)<5)); 
    
    end
    
    clear Decoded_result
end


%% Figure drawing 

nearfar_color = [248 186 0; 170 121 66]/255; 

set(figure(1),'position',[1 1125 288 220]); clf; 

SP = subplot(1,1,1); cla; hold on; 
errorbar(refs, circ_mean(Behav_std_ref'*2*pi/180)*180/pi/2, circ_std(Behav_std_ref'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1), 'ko','capsize',0,'markerfacecolor','w','color',nearfar_color(1,:),'linewidth',1)
plot(refs, circ_mean(Behav_std_ref'*2*pi/180)*180/pi/2, 'ko','linewidth',1.3,'markerfacecolor',nearfar_color(1,:),'color',nearfar_color(1,:)); 

errorbar(refs(2:4), circ_mean(Behav_std_ref(2:4,:)'*2*pi/180)*180/pi/2, circ_std(Behav_std_ref(2:4,:)'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1), 'ko','capsize',0,'markerfacecolor','w','color',nearfar_color(2,:),'linewidth',1)
plot(refs(2:4), circ_mean(Behav_std_ref(2:4,:)'*2*pi/180)*180/pi/2, 'ko','linewidth',1.3,'markerfacecolor',nearfar_color(2,:),'color',nearfar_color(2,:)); 

xlabel('Reference (deg)'); 
ylabel('STD (deg)'); 


set(figure(2),'position',[1 491 1343 314]); clf; 
for iTR = 1:nTR
    SP = subplot(2,7,iTR); cla; hold on; 
    errorbar(refs, circ_mean(squeeze(std_ref_mne_e(iTR,:,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(std_ref_mne_e(iTR,:,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1), 'ko','capsize',0,'markerfacecolor','w','color',nearfar_color(1,:),'linewidth',1)
    plot(refs, circ_mean(squeeze(std_ref_mne_e(iTR,:,:))'*2*pi/180)*180/pi/2, 'ko','linewidth',1.3,'markerfacecolor',nearfar_color(1,:),'color',nearfar_color(1,:)); 
    
    
    errorbar(refs(2:4), circ_mean(squeeze(std_ref_mne_e(iTR,2:4,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(std_ref_mne_e(iTR,2:4,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1), 'ko','capsize',0,'markerfacecolor','w','color',nearfar_color(2,:),'linewidth',1)
    plot(refs(2:4), circ_mean(squeeze(std_ref_mne_e(iTR,2:4,:))'*2*pi/180)*180/pi/2, 'ko','linewidth',1.3,'markerfacecolor',nearfar_color(2,:),'color',nearfar_color(2,:)); 
    
    xlabel('Reference (deg)'); 
    ylabel('STD (deg)'); 
    
    title(['TR=' num2str(iTR)]); 
end



set(figure(3),'position',[1 103 1343 314]); clf; 
for iTR = 1:nTR
    SP = subplot(2,7,iTR); cla; hold on; 
    errorbar(refs, circ_mean(squeeze(std_ref_mne_l(iTR,:,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(std_ref_mne_l(iTR,:,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1), 'ko','capsize',0,'markerfacecolor','w','color',nearfar_color(1,:),'linewidth',1)
    plot(refs, circ_mean(squeeze(std_ref_mne_l(iTR,:,:))'*2*pi/180)*180/pi/2, 'ko','linewidth',1.3,'markerfacecolor',nearfar_color(1,:),'color',nearfar_color(1,:)); 
    
    
    errorbar(refs(2:4), circ_mean(squeeze(std_ref_mne_l(iTR,2:4,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(std_ref_mne_l(iTR,2:4,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1), 'ko','capsize',0,'markerfacecolor','w','color',nearfar_color(2,:),'linewidth',1)
    plot(refs(2:4), circ_mean(squeeze(std_ref_mne_l(iTR,2:4,:))'*2*pi/180)*180/pi/2, 'ko','linewidth',1.3,'markerfacecolor',nearfar_color(2,:),'color',nearfar_color(2,:)); 
    
    xlabel('Reference (deg)'); 
    ylabel('STD (deg)'); 
    
    title(['TR=' num2str(iTR)]); 
end

xx_e = nan(5,nTR); 
xx_l = nan(5,nTR); 
xx_near_e = nan(3,nTR); 
xx_near_l = nan(3,nTR);
xx_far_e = nan(2,nTR); 
xx_far_l = nan(2,nTR); 
for iTR = 1:nTR
    xx_e(:,iTR) = circ_mean(squeeze(std_ref_mne_e(iTR,:,:))'*2*pi/180)*180/pi/2; 
    xx_l(:,iTR) = circ_mean(squeeze(std_ref_mne_l(iTR,:,:))'*2*pi/180)*180/pi/2; 
    xx_near_e(:,iTR) = circ_mean(squeeze(std_ref_mne_e(iTR,2:4,:))'*2*pi/180)*180/pi/2; 
    xx_near_l(:,iTR) = circ_mean(squeeze(std_ref_mne_l(iTR,2:4,:))'*2*pi/180)*180/pi/2; 
    xx_far_e(:,iTR) = circ_mean(squeeze(std_ref_mne_e(iTR,[1, 5],:))'*2*pi/180)*180/pi/2; 
    xx_far_l(:,iTR) = circ_mean(squeeze(std_ref_mne_l(iTR,[1, 5],:))'*2*pi/180)*180/pi/2; 
end

set(figure(4),'position',[1 598 268 207]); clf; 
plot(nanmean(xx_e),'bo-'); hold on; 
plot(nanmean(xx_l),'ro-');
legend('Early', 'Late'); 
ylabel('Baseline STD'); 
xlabel('Time (TR)'); 

set(figure(5),'position',[1 598 268 207]); clf; 
% plot(nanmean(xx_near_e)-nanmean(xx_far_e),'bo-'); hold on; 
% plot(nanmean(xx_near_l)-nanmean(xx_far_l),'ro-');
% plot([0 14],[0 0],'k--'); 
plot(nanmean(xx_near_e),'bo-'); hold on; 
plot(nanmean(xx_far_e),'bo--'); 
plot(nanmean(xx_near_l),'ro-');
plot(nanmean(xx_far_l),'ro--');
% legend('Early', 'Late'); 
ylabel('Near, Far STD'); 

ylabel('Near - Far STD'); 
xlabel('Time (TR)'); 




%% Impact of trial-selection 
std_nf = nan(2,length(sub_list)); 
err_crit = 20; 

for isub = runSub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_estimated/VC_sub-' sub_list(isub,:) '_dec.mat'])
    [nTestTR, nTrials] = size(Decoded_result{1}.est); 
    
    % Behavioral errors
    errme = response - stimulus; 
    errme(errme>90) = errme(errme>90) -180; 
    errme(errme<-90) = errme(errme<-90) +180; 
    
    % Correct for the reference-dependent bias
    errme_corr = nan(size(errme)); 
    for ir = 1:length(refs)
        temp = errme(ref==refs(ir) & ~isnan(errme) & abs(errme)<err_crit) - circ_mean((errme(ref==refs(ir) & ~isnan(errme) & abs(errme)<err_crit))'*2*pi/180)*180/pi/2; 
        temp(temp>90) = temp(temp>90) -180; 
        temp(temp<-90) = temp(temp<-90) +180;
        errme_corr(ref==refs(ir) & ~isnan(errme) & abs(errme)<err_crit) = temp ; 
    end
    
    % Near vs. Far 
    vals = errme_corr(abs(ref)<5  & ~isnan(errme) & abs(errme)<err_crit); 
    std_nf(1,isub) = circ_std(vals'*2*pi/180)*180/pi/2; 
    vals = errme_corr(abs(ref)>5  & ~isnan(errme) & abs(errme)<err_crit); 
    std_nf(2,isub) = circ_std(vals'*2*pi/180)*180/pi/2; 
    
    nTrials = sum(~isnan(errme_corr))/length(errme_corr); 
end


set(figure(10),'position',[1 587 230 218]); clf; 
errorbar(1:2, nanmean(std_nf'), nanstd(std_nf')/sqrt(length(sub_list)-1)); 

[h,p] = ttest(std_nf(1,:), std_nf(2,:)); 

ranges = [5, 11]; 
set(figure(11),'position',[1 587 230 218]); clf; 
SP = subplot(1,1,1); cla; hold on; 
plot(std_nf(1,:), std_nf(2,:),'ko','markerfacecolor','w'); 

plot(ranges,ranges,'k--'); 
xlim(ranges) ; ylim(ranges) ;
xlabel('STD Near'); ylabel('STD Far'); 
title(['p = ' num2str(p)]); 




%% Trial-filtering 
std_fmri_nf_e = nan(2,nTR,length(sub_list)); 
std_fmri_nf_l = nan(2,nTR,length(sub_list)); 

for isub = runSub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_estimated/VC_sub-' sub_list(isub,:) '_dec.mat'])
    [nTestTR, nTrials] = size(Decoded_result{1}.est); 
    
    % Behavioral errors
    errme = response - stimulus; 
    errme(errme>90) = errme(errme>90) -180; 
    errme(errme<-90) = errme(errme<-90) +180; 
    
    % Correct for the reference-dependent bias
    errme_corr = nan(size(errme)); 
    for ir = 1:length(refs)
        temp = errme(ref==refs(ir) & ~isnan(errme) & abs(errme)<err_crit) - circ_mean((errme(ref==refs(ir) & ~isnan(errme) & abs(errme)<err_crit))'*2*pi/180)*180/pi/2; 
        temp(temp>90) = temp(temp>90) -180; 
        temp(temp<-90) = temp(temp<-90) +180; 
        errme_corr(ref==refs(ir) & ~isnan(errme) & abs(errme)<err_crit) = temp ; 
    end
    errme_behav = errme; 
    
    % Sensory decoded
    errme_sensory = nan(length(Decoded_result{3}.est(1,:)), nTR); 
    for iTRx = 1:nTR
        tempx = circ_mean([Decoded_result{3}.est(iTRx,:); Decoded_result{4}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        errme = tempx - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180; 
        errme_sensory(:,iTRx) = errme; 
    end
     % Memory decoded
    errme_memory = nan(length(Decoded_result{3}.est(1,:)), nTR); 
    for iTRx = 1:nTR
        tempx = circ_mean([Decoded_result{9}.est(iTRx,:); Decoded_result{10}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        errme = tempx - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        errme_memory(timing==1,iTRx) = errme(timing==1); 
        
        tempx = circ_mean([Decoded_result{6}.est(iTRx,:); Decoded_result{7}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        errme = tempx - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        errme_memory(timing==2,iTRx) = errme(timing==2); 
    end
    
    % Correct for the reference-dependent bias
    errme_memory_corr = nan(size(errme_memory)); 
    for iTR = 1:nTR
        for ir = 1:length(refs)
            temp = errme_sensory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR) - circ_mean((errme_sensory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1,iTR))*2*pi/180)*180/pi/2; 
            temp(temp>90) = temp(temp>90) -180; 
            temp(temp<-90) = temp(temp<-90) +180; 
            errme_memory_corr(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR) = temp ; 

            temp = errme_sensory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR) - circ_mean((errme_sensory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2,iTR))*2*pi/180)*180/pi/2; 
            temp(temp>90) = temp(temp>90) -180; 
            temp(temp<-90) = temp(temp<-90) +180; 
            errme_memory_corr(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR) = temp ; 
            
%             temp = errme_memory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR) - circ_mean((errme_memory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1,iTR))*2*pi/180)*180/pi/2; 
%             temp(temp>90) = temp(temp>90) -180; 
%             temp(temp<-90) = temp(temp<-90) +180; 
%             errme_memory_corr(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR) = temp ; 
% 
%             temp = errme_memory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR) - circ_mean((errme_memory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2,iTR))*2*pi/180)*180/pi/2; 
%             temp(temp>90) = temp(temp>90) -180; 
%             temp(temp<-90) = temp(temp<-90) +180; 
%             errme_memory_corr(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR) = temp ; 
        end
    end
    
    % Near vs. Far 
    for iTR = 1:nTR
        vals = errme_memory_corr(abs(ref)<5  & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR); 
        std_fmri_nf_e(1,iTR,isub) = circ_std(vals*2*pi/180)*180/pi/2; 
        vals = errme_memory_corr(abs(ref)>5  & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR); 
        std_fmri_nf_e(2,iTR,isub) = circ_std(vals*2*pi/180)*180/pi/2; 
        
        vals = errme_memory_corr(abs(ref)<5  & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR); 
        std_fmri_nf_l(1,iTR,isub) = circ_std(vals*2*pi/180)*180/pi/2; 
        vals = errme_memory_corr(abs(ref)>5  & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR); 
        std_fmri_nf_l(2,iTR,isub) = circ_std(vals*2*pi/180)*180/pi/2; 
    end
    
    
    clear Decoded_result
end

ranges = [20, 40]; 
set(figure(12),'position',[1 491 1426 387]); clf; 
for iTR = 1:nTR
    [h,p] = ttest(squeeze(std_fmri_nf_e(1,iTR,:)), squeeze(std_fmri_nf_e(2,iTR,:))); 

    SP = subplot(2,7,iTR); cla; hold on; 
    plot(squeeze(std_fmri_nf_e(1,iTR,:)), squeeze(std_fmri_nf_e(2,iTR,:)),'ko','markersize',6); 

    plot(ranges,ranges,'k--'); 
    xlim(ranges) ; ylim(ranges) ;
    xlabel('STD Near'); ylabel('STD Far'); 
    if nanmean(squeeze(std_fmri_nf_e(1,iTR,:)))-nanmean(squeeze(std_fmri_nf_e(2,iTR,:)))>0
        title(['Near > Far, p = ' num2str(p)]); 
    else
        title(['Near < Far, p = ' num2str(p)]); 
    end
end


ranges = [20, 40]; 
set(figure(13),'position',[1 1 1426 387]); clf; 
for iTR = 1:nTR
    [h,p] = ttest(squeeze(std_fmri_nf_l(1,iTR,:)), squeeze(std_fmri_nf_l(2,iTR,:))); 

    SP = subplot(2,7,iTR); cla; hold on; 
    plot(squeeze(std_fmri_nf_l(1,iTR,:)), squeeze(std_fmri_nf_l(2,iTR,:)),'ko','markersize',6); 

    plot(ranges,ranges,'k--'); 
    xlim(ranges) ; ylim(ranges) ;
    xlabel('STD Near'); ylabel('STD Far'); 
    if nanmean(squeeze(std_fmri_nf_l(1,iTR,:)))-nanmean(squeeze(std_fmri_nf_l(2,iTR,:)))>0
        title(['Near > Far, p = ' num2str(p)]); 
    else
        title(['Near < Far, p = ' num2str(p)]); 
    end
end

for iTR = 1:nTR
    [h,p,ci,stats] = ttest(squeeze(std_fmri_nf_e(1,iTR,:)), squeeze(std_fmri_nf_e(2,iTR,:))); 
    t_stats(iTR) = stats.tstat; 
    p_value(iTR) = p; 
end

set(figure(14),'position',[14 602 310 203]); clf; 
SP = subplot(1,1,1); cla; hold on; 
plot(t_stats,'ko-');
for iTR = 1:nTR
    if p_value(iTR) < 0.05
        plot(iTR, t_stats(iTR),'ko','markerfacecolor','k'); 
    else
        plot(iTR, t_stats(iTR),'ko','markerfacecolor','w'); 
    end
end
plot([0 14],[0 0], 'k--'); 
xlim([0,14]); ylim([-2 4]); 
xlabel('Time (TR)'); 
ylabel('t value (near - far std)'); 





%% Orientation-specific correction
err_crit = 20; 

for isub = runSub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_estimated/VC_sub-' sub_list(isub,:) '_dec.mat'])
    [nTestTR, nTrials] = size(Decoded_result{1}.est); 
    
    % Behavioral errors
    errme = response - stimulus; 
    errme(errme>90) = errme(errme>90) -180; 
    errme(errme<-90) = errme(errme<-90) +180; 
    
    % Correct for orientation-dependent bias
    errme_ori = nan(size(errme)); 
    for istim = 1:length(stimcond)
        temp = errme(stimulus==stimcond(istim) & ~isnan(errme) & abs(errme)<err_crit) - circ_mean((errme(stimulus==stimcond(istim) & ~isnan(errme) & abs(errme)<err_crit))'*2*pi/180)*180/pi/2; 
        temp(temp>90) = temp(temp>90) -180; 
        temp(temp<-90) = temp(temp<-90) +180; 
        errme_ori(stimulus==stimcond(istim) & ~isnan(errme) & abs(errme)<err_crit) = temp ; 
    end
    errme = errme_ori; 
    
    % Correct for the reference-dependent bias
    errme_corr = nan(size(errme)); 
    for ir = 1:length(refs)
        temp = errme(ref==refs(ir) & ~isnan(errme) & abs(errme)<err_crit) - circ_mean((errme(ref==refs(ir) & ~isnan(errme) & abs(errme)<err_crit))'*2*pi/180)*180/pi/2; 
        temp(temp>90) = temp(temp>90) -180; 
        temp(temp<-90) = temp(temp<-90) +180; 
        errme_corr(ref==refs(ir) & ~isnan(errme) & abs(errme)<err_crit) = temp ; 
    end
    errme_behav = errme_corr; 
    errme = errme_corr;
    

    for ir = 1:length(refs)
        Behav_mean_ref(ir,isub) = circ_mean((errme_ori(ref==refs(ir) & ref==refs(ir)  & ~isnan(errme_ori) & abs(errme_ori)<err_crit ))'*2*pi/180)*180/pi/2;
        Behav_std_ref(ir,isub) = circ_std((errme_ori(ref==refs(ir) & ref==refs(ir)  & ~isnan(errme_ori) & abs(errme_ori)<err_crit))'*2*pi/180)*180/pi/2;
    end
    
    
    % Sensory decoded
    errme_sensory = nan(length(Decoded_result{3}.est(1,:)), nTR); 
    for iTRx = 1:nTR
        tempx = circ_mean([Decoded_result{3}.est(iTRx,:); Decoded_result{4}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        errme = tempx - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180; 
        errme_sensory(:,iTRx) = errme; 
    end
     % Memory decoded
    errme_memory = nan(length(Decoded_result{3}.est(1,:)), nTR); 
    for iTRx = 1:nTR
        tempx = circ_mean([Decoded_result{9}.est(iTRx,:); Decoded_result{10}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        errme = tempx - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        errme_memory(timing==1,iTRx) = errme(timing==1); 
        
        tempx = circ_mean([Decoded_result{6}.est(iTRx,:); Decoded_result{7}.est(iTRx,:)]*2*pi/180)*180/pi/2;
        tempx(tempx<0) = tempx(tempx<0)+180; 
        errme = tempx - stimulus; 
        errme(errme>90) = errme(errme>90) -180; 
        errme(errme<-90) = errme(errme<-90) +180;  
        errme_memory(timing==2,iTRx) = errme(timing==2); 
    end
    
    % Correct for the reference-dependent bias
    errme_memory_corr = nan(size(errme_memory)); 
    for iTR = 1:nTR
        for ir = 1:length(refs)
%             temp = errme_sensory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR) - circ_mean((errme_sensory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1,iTR))*2*pi/180)*180/pi/2; 
%             temp(temp>90) = temp(temp>90) -180; 
%             temp(temp<-90) = temp(temp<-90) +180; 
%             errme_memory_corr(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR) = temp ; 
% 
%             temp = errme_sensory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR) - circ_mean((errme_sensory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2,iTR))*2*pi/180)*180/pi/2; 
%             temp(temp>90) = temp(temp>90) -180; 
%             temp(temp<-90) = temp(temp<-90) +180; 
%             errme_memory_corr(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR) = temp ; 
            
            temp = errme_memory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR) - circ_mean((errme_memory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1,iTR))*2*pi/180)*180/pi/2; 
            temp(temp>90) = temp(temp>90) -180; 
            temp(temp<-90) = temp(temp<-90) +180; 
            errme_memory_corr(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR) = temp ; 

            temp = errme_memory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR) - circ_mean((errme_memory(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2,iTR))*2*pi/180)*180/pi/2; 
            temp(temp>90) = temp(temp>90) -180; 
            temp(temp<-90) = temp(temp<-90) +180; 
            errme_memory_corr(ref==refs(ir) & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR) = temp ; 
        end
    end
    
    % Near vs. Far 
    for iTR = 1:nTR
        vals = errme_memory_corr(abs(ref)<5  & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR); 
        std_fmri_nf_e(1,iTR,isub) = circ_std(vals*2*pi/180)*180/pi/2; 
        vals = errme_memory_corr(abs(ref)>5  & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==1, iTR); 
        std_fmri_nf_e(2,iTR,isub) = circ_std(vals*2*pi/180)*180/pi/2; 
        
        vals = errme_memory_corr(abs(ref)<5  & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR); 
        std_fmri_nf_l(1,iTR,isub) = circ_std(vals*2*pi/180)*180/pi/2; 
        vals = errme_memory_corr(abs(ref)>5  & ~isnan(errme_behav) & abs(errme_behav)<err_crit & timing==2, iTR); 
        std_fmri_nf_l(2,iTR,isub) = circ_std(vals*2*pi/180)*180/pi/2; 
    end
    
    
    clear Decoded_result
end


set(figure(101),'position',[1 1125 288 220]); clf; 

SP = subplot(1,1,1); cla; hold on; 
errorbar(refs, circ_mean(Behav_std_ref'*2*pi/180)*180/pi/2, circ_std(Behav_std_ref'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1), 'ko','capsize',0,'markerfacecolor','w','color',nearfar_color(1,:),'linewidth',1)
plot(refs, circ_mean(Behav_std_ref'*2*pi/180)*180/pi/2, 'ko','linewidth',1.3,'markerfacecolor',nearfar_color(1,:),'color',nearfar_color(1,:)); 

errorbar(refs(2:4), circ_mean(Behav_std_ref(2:4,:)'*2*pi/180)*180/pi/2, circ_std(Behav_std_ref(2:4,:)'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1), 'ko','capsize',0,'markerfacecolor','w','color',nearfar_color(2,:),'linewidth',1)
plot(refs(2:4), circ_mean(Behav_std_ref(2:4,:)'*2*pi/180)*180/pi/2, 'ko','linewidth',1.3,'markerfacecolor',nearfar_color(2,:),'color',nearfar_color(2,:)); 

xlabel('Reference (deg)'); 
ylabel('STD (deg)'); 














