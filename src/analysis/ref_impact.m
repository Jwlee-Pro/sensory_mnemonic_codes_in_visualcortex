clear all; close all; clc;

% Load sub/roi-list 
load('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/sub_list.mat')
addpath('/Users/joonwon/Dropbox/0_PRESENTATION_FILES_JWL2020/BRL_daily/20211007_/WorkingMemoryDynamics-master/CircStat2012a')
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



xval = [xvec flip(xvec)]; 
ymean = circ_mean(squeeze(temp_choice_TR_e(:,1,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_TR_e(:,1,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'r','facealpha',0.5,'edgecolor','none'); 
    
ymean = circ_mean(squeeze(temp_choice_TR_e(:,2,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_TR_e(:,2,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd +flip(ymean+ystd)] ; 
patch(xval, yval, 'b','facealpha',0.5,'edgecolor','none'); 
plot([0 15],[0 0],'k--','linewidth',1.5); 
% ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 
title('Early Dm'); 