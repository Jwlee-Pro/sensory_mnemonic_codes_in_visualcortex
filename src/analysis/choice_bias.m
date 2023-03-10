

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


timing_m = []; 
ref_m = []; 
choice_m = []; 
dec_tr = []; 
dec_sen = []; 
dec_mne = []; 
for isub = runSub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_estimated/VC_sub-' sub_list(isub,:) '_dec.mat'])
    [nTestTR, nTrials] = size(Decoded_result{1}.est); 
    
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
        avgPred_Near_e(iTRx,isub) = (sum(dec_dv(timing==1)>0 & choice(timing==1)==1 & abs(ref)<5) + sum(dec_dv(timing==1)<0 & choice(timing==1)==2 & abs(ref)<5))/length(choice(timing==1 & abs(ref)<5)); 
        avgPred_Near_l(iTRx,isub) = (sum(dec_dv(timing==2)>0 & choice(timing==2)==1 & abs(ref)<5) + sum(dec_dv(timing==2)<0 & choice(timing==2)==2 & abs(ref)<5))/length(choice(timing==2 & abs(ref)<5)); 
        avgPred_All_e(iTRx,isub) = (sum(dec_dv(timing==1)>0 & choice(timing==1)==1 & abs(ref)<5) + sum(dec_dv(timing==1)<0 & choice(timing==1)==2 & abs(ref)<5))/length(choice(timing==1 & abs(ref)<5)); 
        avgPred_All_l(iTRx,isub) = (sum(dec_dv(timing==2)>0 & choice(timing==2)==1 & abs(ref)<5) + sum(dec_dv(timing==2)<0 & choice(timing==2)==2 & abs(ref)<5))/length(choice(timing==2 & abs(ref)<5)); 
    
    end
    
    clear Decoded_result
end

asdf
%% Figure drawing


% xvec = (1:nTR) * 2 -1; 
xvec = (1:nTR); 

set(figure(1),'position',[1 779 397 566]); clf; 
SP = subplot(3,2,1); cla; hold on; 
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

SP = subplot(3,2,2); cla; hold on; 
xval = [xvec flip(xvec)]; 
ymean = circ_mean(squeeze(temp_choice_TR_l(:,1,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_TR_l(:,1,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'r','facealpha',0.5,'edgecolor','none'); 
    
ymean = circ_mean(squeeze(temp_choice_TR_l(:,2,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_TR_l(:,2,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd +flip(ymean+ystd)] ; 
patch(xval, yval, 'b','facealpha',0.5,'edgecolor','none'); 
plot([0 15],[0 0],'k--','linewidth',1.5); 
title('Late Dm'); 
% ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 


% Sensory 
SP = subplot(3,2,3); cla; hold on; 
xval = [xvec flip(xvec)]; 
ymean = circ_mean(squeeze(temp_choice_sen_e(:,1,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_sen_e(:,1,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'r','facealpha',0.5,'edgecolor','none'); 
    
ymean = circ_mean(squeeze(temp_choice_sen_e(:,2,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_sen_e(:,2,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd +flip(ymean+ystd)] ; 
patch(xval, yval, 'b','facealpha',0.5,'edgecolor','none'); 
plot([0 15],[0 0],'k--','linewidth',1.5); 
% ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,4); cla; hold on; 
xval = [xvec flip(xvec)]; 
ymean = circ_mean(squeeze(temp_choice_sen_l(:,1,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_sen_l(:,1,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'r','facealpha',0.5,'edgecolor','none'); 
    
ymean = circ_mean(squeeze(temp_choice_sen_l(:,2,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_sen_l(:,2,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd +flip(ymean+ystd)] ; 
patch(xval, yval, 'b','facealpha',0.5,'edgecolor','none'); 
plot([0 15],[0 0],'k--','linewidth',1.5); 
title('Late Dm'); 
% ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 



% Mnemonic 
SP = subplot(3,2,5); cla; hold on; 
xval = [xvec flip(xvec)]; 
ymean = circ_mean(squeeze(temp_choice_mne_e(:,1,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_mne_e(:,1,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'r','facealpha',0.5,'edgecolor','none'); 
    
ymean = circ_mean(squeeze(temp_choice_mne_e(:,2,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_mne_e(:,2,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd +flip(ymean+ystd)] ; 
patch(xval, yval, 'b','facealpha',0.5,'edgecolor','none'); 
plot([0 15],[0 0],'k--','linewidth',1.5); 
% ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,6); cla; hold on; 
xval = [xvec flip(xvec)]; 
ymean = circ_mean(squeeze(temp_choice_mne_l(:,1,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_mne_l(:,1,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'r','facealpha',0.5,'edgecolor','none'); 
    
ymean = circ_mean(squeeze(temp_choice_mne_l(:,2,:))'*2*pi/180)*180/pi/2; 
ystd = circ_std(squeeze(temp_choice_mne_l(:,2,:))'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd +flip(ymean+ystd)] ; 
patch(xval, yval, 'b','facealpha',0.5,'edgecolor','none'); 
plot([0 15],[0 0],'k--','linewidth',1.5); 
title('Late Dm'); 
% ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 






set(figure(2),'position',[1 779 397 566]); clf; 
SP = subplot(3,2,1); cla; hold on; 
xval = [xvec flip(xvec)]; 
xx = squeeze(temp_choice_TR_e(:,1,:)-temp_choice_TR_e(:,2,:)); 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180;

ymean = circ_mean(xx'*2*pi/180)*180/pi/2; 
ystd = circ_std(xx'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'k','facealpha',0.5,'edgecolor','none'); 
    
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('\Delta decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,2); cla; hold on; 
xval = [xvec flip(xvec)]; 
xx = squeeze(temp_choice_TR_l(:,1,:)-temp_choice_TR_l(:,2,:)); 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180;

ymean = circ_mean(xx'*2*pi/180)*180/pi/2; 
ystd = circ_std(xx'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'k','facealpha',0.5,'edgecolor','none'); 
    
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
title('Late Dm'); 

SP = subplot(3,2,3); cla; hold on; 
xval = [xvec flip(xvec)]; 
xx = squeeze(temp_choice_sen_e(:,1,:)-temp_choice_sen_e(:,2,:)); 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180;

ymean = circ_mean(xx'*2*pi/180)*180/pi/2; 
ystd = circ_std(xx'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'k','facealpha',0.5,'edgecolor','none'); 
    
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('\Delta decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,4); cla; hold on; 
xval = [xvec flip(xvec)]; 
xx = squeeze(temp_choice_sen_l(:,1,:)-temp_choice_sen_l(:,2,:)); 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180;

ymean = circ_mean(xx'*2*pi/180)*180/pi/2; 
ystd = circ_std(xx'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'k','facealpha',0.5,'edgecolor','none'); 
    
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
title('Late Dm'); 


SP = subplot(3,2,5); cla; hold on; 
xval = [xvec flip(xvec)]; 
xx = squeeze(temp_choice_mne_e(:,1,:)-temp_choice_mne_e(:,2,:)); 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180;

ymean = circ_mean(xx'*2*pi/180)*180/pi/2; 
ystd = circ_std(xx'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'k','facealpha',0.5,'edgecolor','none'); 
    
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('\Delta decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,6); cla; hold on; 
xval = [xvec flip(xvec)]; 
xx = squeeze(temp_choice_mne_l(:,1,:)-temp_choice_mne_l(:,2,:)); 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180;

ymean = circ_mean(xx'*2*pi/180)*180/pi/2; 
ystd = circ_std(xx'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
yval = [ymean-ystd flip(ymean+ystd)] ; 
patch(xval, yval, 'k','facealpha',0.5,'edgecolor','none'); 
    
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
title('Late Dm'); 



%% Merged

set(figure(3),'position',[1 779 397 566]); clf; 

SP = subplot(3,2,1); cla; hold on; 
plot(rangeu, circ_mean(dec_tr(timing_m==1 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_tr(timing_m==1 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,2); cla; hold on; 
plot(rangeu, circ_mean(dec_tr(timing_m==2 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_tr(timing_m==2 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
title('Late Dm'); 
ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 

SP = subplot(3,2,3); cla; hold on; 
plot(rangeu, circ_mean(dec_sen(timing_m==1 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_sen(timing_m==1 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,4); cla; hold on; 
plot(rangeu, circ_mean(dec_sen(timing_m==2 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_sen(timing_m==2 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
title('Late Dm'); 
ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 

SP = subplot(3,2,5); cla; hold on; 
plot(rangeu, circ_mean(dec_mne(timing_m==1 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_mne(timing_m==1 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,6); cla; hold on; 
plot(rangeu, circ_mean(dec_mne(timing_m==2 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_mne(timing_m==2 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
title('Late Dm'); 
ylim([-20 20]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 





set(figure(4),'position',[1 779 397 566]); clf; 

SP = subplot(3,2,1); cla; hold on; 
xx = circ_mean(dec_tr(timing_m==1 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2 - ...
    circ_mean(dec_tr(timing_m==1 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2; 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180; 
plot(rangeu, xx, 'r.-','color','k','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,2); cla; hold on; 
xx = circ_mean(dec_tr(timing_m==2 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2 - ...
    circ_mean(dec_tr(timing_m==2 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2; 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180; 
plot(rangeu, xx, 'r.-','color','k','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
title('Late Dm'); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 


SP = subplot(3,2,3); cla; hold on; 
xx = circ_mean(dec_sen(timing_m==1 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2 - ...
    circ_mean(dec_sen(timing_m==1 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2; 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180; 
plot(rangeu, xx, 'r.-','color','k','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,4); cla; hold on; 
xx = circ_mean(dec_sen(timing_m==2 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2 - ...
    circ_mean(dec_sen(timing_m==2 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2; 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180; 
plot(rangeu, xx, 'r.-','color','k','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
title('Late Dm'); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 


SP = subplot(3,2,5); cla; hold on; 
xx = circ_mean(dec_mne(timing_m==1 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2 - ...
    circ_mean(dec_mne(timing_m==1 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2; 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180; 
plot(rangeu, xx, 'r.-','color','k','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 
ylabel('Decoded bias (deg)'); 
title('Early Dm'); 

SP = subplot(3,2,6); cla; hold on; 
xx = circ_mean(dec_mne(timing_m==2 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2 - ...
    circ_mean(dec_mne(timing_m==2 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2; 
xx(xx>90) = xx(xx>90) -180; 
xx(xx<-90) = xx(xx<-90) +180; 
plot(rangeu, xx, 'r.-','color','k','markersize',13,'linewidth',2);
plot([0 15],[0 0],'k--','linewidth',1.5); 
title('Late Dm'); 
ylim([-10 30]); 
xlim([3 15]); 
xlabel('Time (TR)'); 

%% 
% Reference conditions
refcond = [19 47 120; 147 189 234; 190 190 190; 238 139 129; 120 14 32]/255; 

rangeu = 1:(nTR); 

% Sensory code
set(figure(101),'position',[1 560 742 245]); clf; 
SP = subplot(1,2,1); cla; hold on; 
for ir = 1:5
    plot(rangeu, circ_mean(dec_sen(timing_m==1 & ref_m ==refs(ir),rangeu)*2*pi/180)*180/pi/2, 'r.-','color',refcond(ir,:),'markersize',13,'linewidth',2);
    plot([0 14],[0 0]+refs(ir),'k-','color',refcond(ir,:)); 
end
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Early Dm'); 
ylim([-30 30]);
SP = subplot(1,2,2); cla; hold on; 
for ir = 1:5
    plot(rangeu, circ_mean(dec_sen(timing_m==2 & ref_m ==refs(ir),rangeu)*2*pi/180)*180/pi/2, 'r.-','color',refcond(ir,:),'markersize',13,'linewidth',2);
    plot([0 14],[0 0]+refs(ir),'k-','color',refcond(ir,:)); 
end
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Late Dm'); 
ylim([-30 30]);

% Mnemonic code
set(figure(102),'position',[1 560 742 245]); clf; 
SP = subplot(1,2,1); cla; hold on; 
for ir = 1:5
    plot(rangeu, circ_mean(dec_mne(timing_m==1 & ref_m ==refs(ir),rangeu)*2*pi/180)*180/pi/2, 'r.-','color',refcond(ir,:),'markersize',13,'linewidth',2);
    plot([0 14],[0 0]+refs(ir),'k-','color',refcond(ir,:)); 
end
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Early Dm'); 
ylim([-30 30]);
SP = subplot(1,2,2); cla; hold on; 
for ir = 1:5
    plot(rangeu, circ_mean(dec_mne(timing_m==2 & ref_m ==refs(ir),rangeu)*2*pi/180)*180/pi/2, 'r.-','color',refcond(ir,:),'markersize',13,'linewidth',2);
    plot([0 14],[0 0]+refs(ir),'k-','color',refcond(ir,:)); 
end
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Late Dm'); 
ylim([-30 30]);


% Mixed code
set(figure(103),'position',[1 560 742 245]); clf; 
SP = subplot(1,2,1); cla; hold on; 
for ir = 1:5
    plot(rangeu, circ_mean(dec_tr(timing_m==1 & ref_m ==refs(ir),rangeu)*2*pi/180)*180/pi/2, 'r.-','color',refcond(ir,:),'markersize',13,'linewidth',2);
    plot([0 14],[0 0]+refs(ir),'k-','color',refcond(ir,:)); 
end
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Early Dm'); 
ylim([-30 30]);
SP = subplot(1,2,2); cla; hold on; 
for ir = 1:5
    plot(rangeu, circ_mean(dec_tr(timing_m==2 & ref_m ==refs(ir),rangeu)*2*pi/180)*180/pi/2, 'r.-','color',refcond(ir,:),'markersize',13,'linewidth',2);
    plot([0 14],[0 0]+refs(ir),'k-','color',refcond(ir,:)); 
end
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Late Dm'); 
ylim([-30 30]);


%% Choice
% Sensory code
set(figure(105),'position',[1 560 742 245]); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(rangeu, circ_mean(dec_sen(timing_m==1 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_sen(timing_m==1 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Early Dm'); 
ylim([-30 30]);
SP = subplot(1,2,2); cla; hold on; 
plot(rangeu, circ_mean(dec_sen(timing_m==2 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_sen(timing_m==2 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Late Dm'); 
ylim([-30 30]);


set(figure(106),'position',[1 241 742 245]); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(rangeu, circ_mean(dec_mne(timing_m==1 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_mne(timing_m==1 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Early Dm'); 
ylim([-30 30]); 
SP = subplot(1,2,2); cla; hold on; 
plot(rangeu, circ_mean(dec_mne(timing_m==2 & choice_m ==1 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_mne(timing_m==2 & choice_m ==2 & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Late Dm'); 
ylim([-30 30]);

set(figure(107),'position',[1 8 742 245]); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(rangeu, circ_mean(dec_tr(timing_m==1 & choice_m ==1  & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_tr(timing_m==1 & choice_m ==2  & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Early Dm'); 
ylim([-30 30]); 

SP = subplot(1,2,2); cla; hold on; 
plot(rangeu, circ_mean(dec_tr(timing_m==2 & choice_m ==1  & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_tr(timing_m==2 & choice_m ==2  & abs(ref_m)<5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Late Dm'); 
ylim([-30 30]);








set(figure(105),'position',[1 560 742 245]); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(rangeu, circ_mean(dec_sen(timing_m==1 & choice_m ==1 & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_sen(timing_m==1 & choice_m ==2 & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Early Dm'); 
ylim([-30 30]);
SP = subplot(1,2,2); cla; hold on; 
plot(rangeu, circ_mean(dec_sen(timing_m==2 & choice_m ==1 & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_sen(timing_m==2 & choice_m ==2 & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Late Dm'); 
ylim([-30 30]);


set(figure(106),'position',[1 241 742 245]); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(rangeu, circ_mean(dec_mne(timing_m==1 & choice_m ==1 & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_mne(timing_m==1 & choice_m ==2 & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Early Dm'); 
ylim([-30 30]); 
SP = subplot(1,2,2); cla; hold on; 
plot(rangeu, circ_mean(dec_mne(timing_m==2 & choice_m ==1 & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_mne(timing_m==2 & choice_m ==2 & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Late Dm'); 
ylim([-30 30]);

set(figure(107),'position',[1 8 742 245]); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(rangeu, circ_mean(dec_tr(timing_m==1 & choice_m ==1  & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_tr(timing_m==1 & choice_m ==2  & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Early Dm'); 
ylim([-30 30]); 

SP = subplot(1,2,2); cla; hold on; 
plot(rangeu, circ_mean(dec_tr(timing_m==2 & choice_m ==1  & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','r','markersize',13,'linewidth',2);
plot(rangeu, circ_mean(dec_tr(timing_m==2 & choice_m ==2  & abs(ref_m)>5,rangeu)*2*pi/180)*180/pi/2, 'r.-','color','b','markersize',13,'linewidth',2);
plot([0 14],[0 0],'k-','color',refcond(3,:)); 
xlabel('Time (TR)'); ylabel('Decoded bias'); 
title('Late Dm'); 
ylim([-30 30]);
% figure(50); clf; 
% plot(double(stimulus), Decoded_result{13}.est(13,:),'ko'); 
% 
% errme = Decoded_result{13}.est(13,:) - stimulus;
% errme(errme>90) = errme(errme>90) -180; 
% errme(errme<-90) = errme(errme<-90) +180;  
% 
% figure(51); clf; 
% SP = subplot(2,1,1); cla; hold on; 
% plot(stimulus,'ko','markerfacecolor','k'); 
% plot(Decoded_result{13}.est(13,:),'ro','linewidth',2); 
% xlim([0, length(stimulus)]); 
% SP = subplot(2,1,2); cla; hold on; 
% plot(errme,'bo','markerfacecolor','k'); 
% plot([0 length(stimulus)+1], [0 0],'k--'); 
% xlim([0, length(stimulus)]); 

%% 

asdf


figure(24); clf; 
SP = subplot(3,2,1); cla; hold on; 
plot(rangeu, circ_mean(Sen_dec_choice_1_e(:,rangeu)*2*pi/180)*180/pi/2,'ro-','markerfacecolor','w')
plot(rangeu, circ_mean(Sen_dec_choice_2_e(:,rangeu)*2*pi/180)*180/pi/2,'bo-','markerfacecolor','w')

SP = subplot(3,2,2); cla; hold on; 
plot(rangeu, circ_mean(Sen_dec_choice_1_l(:,rangeu)*2*pi/180)*180/pi/2,'ro-','markerfacecolor','w')
plot(rangeu, circ_mean(Sen_dec_choice_2_l(:,rangeu)*2*pi/180)*180/pi/2,'bo-','markerfacecolor','w')

SP = subplot(3,2,3); cla; hold on; 
plot(rangeu, circ_mean(Mne_dec_choice_1_e(:,rangeu)*2*pi/180)*180/pi/2,'ro-','markerfacecolor','w')
plot(rangeu, circ_mean(Mne_dec_choice_2_e(:,rangeu)*2*pi/180)*180/pi/2,'bo-','markerfacecolor','w')

SP = subplot(3,2,4); cla; hold on; 
plot(rangeu, circ_mean(Mne_dec_choice_1_l(:,rangeu)*2*pi/180)*180/pi/2,'ro-','markerfacecolor','w')
plot(rangeu, circ_mean(Mne_dec_choice_2_l(:,rangeu)*2*pi/180)*180/pi/2,'bo-','markerfacecolor','w')

SP = subplot(3,2,5); cla; hold on; 
plot(rangeu, circ_mean(TR_dec_choice_1_e(:,rangeu)*2*pi/180)*180/pi/2,'ro-','markerfacecolor','w')
plot(rangeu, circ_mean(TR_dec_choice_2_e(:,rangeu)*2*pi/180)*180/pi/2,'bo-','markerfacecolor','w')

SP = subplot(3,2,6); cla; hold on; 
plot(rangeu, circ_mean(TR_dec_choice_1_l(:,rangeu)*2*pi/180)*180/pi/2,'ro-','markerfacecolor','w')
plot(rangeu, circ_mean(TR_dec_choice_2_l(:,rangeu)*2*pi/180)*180/pi/2,'bo-','markerfacecolor','w')

%%




figure(13) ;clf; 
SP = subplot(1,1,1); cla; hold on; 
plot(nanmean(avgPred_e'),'b-'); 
plot(nanmean(avgPred_l'),'r-'); 

% asdf

figure(4); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(circ_mean(temp_mean_e'*2*pi/180)*180/pi/2)
SP = subplot(1,2,2); cla; hold on; 
plot(circ_mean(temp_mean_l'*2*pi/180)*180/pi/2)


set(figure(5),'position',[440 623 865 175]); clf; 
SP = subplot(1,3,1); cla; hold on; 
plot(circ_mean(temp_sen_mean_e'*2*pi/180)*180/pi/2,'b-')
plot(circ_mean(temp_sen_mean_l'*2*pi/180)*180/pi/2,'r-')
title('sensory code'); 
SP = subplot(1,3,2); cla; hold on; 
plot(circ_mean(temp_mne_mean_e'*2*pi/180)*180/pi/2,'b-')
plot(circ_mean(temp_mne_mean_l'*2*pi/180)*180/pi/2,'r-')
title('mnemonic code'); 
SP = subplot(1,3,3); cla; hold on; 
plot(circ_mean(temp_TR_mean_e'*2*pi/180)*180/pi/2,'b-')
plot(circ_mean(temp_TR_mean_l'*2*pi/180)*180/pi/2,'r-')
title('mixed code'); 


set(figure(6),'position',[440 322 560 476]); clf; 
SP = subplot(3,2,1); cla; hold on; 
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_sen_e(4:end,1,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_sen_e(4:end,1,:))'*2*pi/180)*180/pi/2/sqrt(49),'ro-')
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_sen_e(4:end,2,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_sen_e(4:end,2,:))'*2*pi/180)*180/pi/2/sqrt(49),'bo-')
% plot(4:14, circ_mean(squeeze(temp_choice_sen_e(4:end,1,:))'*2*pi/180)*180/pi/2,'r-')
% plot(4:14, circ_mean(squeeze(temp_choice_sen_e(4:end,2,:))'*2*pi/180)*180/pi/2,'b-')
title('sensory code/ early dm'); 

SP = subplot(3,2,2); cla; hold on; 
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_sen_l(4:end,1,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_sen_l(4:end,1,:))'*2*pi/180)*180/pi/2/sqrt(49),'ro-')
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_sen_l(4:end,2,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_sen_l(4:end,2,:))'*2*pi/180)*180/pi/2/sqrt(49),'bo-')
% plot((4:14)*2, circ_mean(squeeze(temp_choice_sen_l(4:end,1,:))'*2*pi/180)*180/pi/2,'r-')
% plot((4:14)*2, circ_mean(squeeze(temp_choice_sen_l(4:end,2,:))'*2*pi/180)*180/pi/2,'b-')
title('sensory code/ late dm'); 

SP = subplot(3,2,3); cla; hold on; 
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_mne_e(4:end,1,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_mne_e(4:end,1,:))'*2*pi/180)*180/pi/2/sqrt(49),'ro-')
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_mne_e(4:end,2,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_mne_e(4:end,2,:))'*2*pi/180)*180/pi/2/sqrt(49),'bo-')

% plot((4:14)*2, circ_mean(squeeze(temp_choice_mne_e(4:end,1,:))'*2*pi/180)*180/pi/2,'r-')
% plot((4:14)*2, circ_mean(squeeze(temp_choice_mne_e(4:end,2,:))'*2*pi/180)*180/pi/2,'b-')
title('mnemonic code/ early dm'); 

SP = subplot(3,2,4); cla; hold on; 
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_mne_l(4:end,1,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_mne_l(4:end,1,:))'*2*pi/180)*180/pi/2/sqrt(49),'ro-')
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_mne_l(4:end,2,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_mne_l(4:end,2,:))'*2*pi/180)*180/pi/2/sqrt(49),'bo-')

% plot((4:14)*2, circ_mean(squeeze(temp_choice_mne_l(4:end,1,:))'*2*pi/180)*180/pi/2,'r-')
% plot((4:14)*2, circ_mean(squeeze(temp_choice_mne_l(4:end,2,:))'*2*pi/180)*180/pi/2,'b-')
title('mnemonic code/ late dm'); 

SP = subplot(3,2,5); cla; hold on; 
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_TR_e(4:end,1,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_TR_e(4:end,1,:))'*2*pi/180)*180/pi/2/sqrt(49),'ro-')
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_TR_e(4:end,2,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_TR_e(4:end,2,:))'*2*pi/180)*180/pi/2/sqrt(49),'bo-')
% plot((4:14)*2, circ_mean(squeeze(temp_choice_TR_e(4:end,1,:))'*2*pi/180)*180/pi/2,'r-')
% plot((4:14)*2, circ_mean(squeeze(temp_choice_TR_e(4:end,2,:))'*2*pi/180)*180/pi/2,'b-')
title('mixed code/ early dm'); 

SP = subplot(3,2,6); cla; hold on; 
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_TR_l(4:end,1,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_TR_l(4:end,1,:))'*2*pi/180)*180/pi/2/sqrt(49),'ro-')
errorbar((4:14)*2, circ_mean(squeeze(temp_choice_TR_l(4:end,2,:))'*2*pi/180)*180/pi/2, circ_std(squeeze(temp_choice_TR_l(4:end,2,:))'*2*pi/180)*180/pi/2/sqrt(49),'bo-')
% plot((4:14)*2, circ_mean(squeeze(temp_choice_TR_l(4:end,1,:))'*2*pi/180)*180/pi/2,'r-')
% plot((4:14)*2, circ_mean(squeeze(temp_choice_TR_l(4:end,2,:))'*2*pi/180)*180/pi/2,'b-')
title('mixed code/ late dm'); 





set(figure(7),'position',[440 322 560 476]); clf; 
SP = subplot(3,2,1); cla; hold on; 
plot(circ_mean(squeeze(tempx_choice_sen_e(:,1,:))'*2*pi/180)*180/pi/2,'r-')
plot(circ_mean(squeeze(tempx_choice_sen_e(:,2,:))'*2*pi/180)*180/pi/2,'b-')
title('sensory code/ early dm'); 

SP = subplot(3,2,2); cla; hold on; 
plot(circ_mean(squeeze(tempx_choice_sen_l(:,1,:))'*2*pi/180)*180/pi/2,'r-')
plot(circ_mean(squeeze(tempx_choice_sen_l(:,2,:))'*2*pi/180)*180/pi/2,'b-')
title('sensory code/ late dm'); 

SP = subplot(3,2,3); cla; hold on; 
plot(circ_mean(squeeze(tempx_choice_mne_e(:,1,:))'*2*pi/180)*180/pi/2,'r-')
plot(circ_mean(squeeze(tempx_choice_mne_e(:,2,:))'*2*pi/180)*180/pi/2,'b-')
title('mnemonic code/ early dm'); 

SP = subplot(3,2,4); cla; hold on; 
plot(circ_mean(squeeze(tempx_choice_mne_l(:,1,:))'*2*pi/180)*180/pi/2,'r-')
plot(circ_mean(squeeze(tempx_choice_mne_l(:,2,:))'*2*pi/180)*180/pi/2,'b-')
title('mnemonic code/ late dm'); 

SP = subplot(3,2,5); cla; hold on; 
plot(circ_mean(squeeze(tempx_choice_TR_e(:,1,:))'*2*pi/180)*180/pi/2,'r-')
plot(circ_mean(squeeze(tempx_choice_TR_e(:,2,:))'*2*pi/180)*180/pi/2,'b-')
title('mixed code/ early dm'); 

SP = subplot(3,2,6); cla; hold on; 
plot(circ_mean(squeeze(tempx_choice_TR_l(:,1,:))'*2*pi/180)*180/pi/2,'r-')
plot(circ_mean(squeeze(tempx_choice_TR_l(:,2,:))'*2*pi/180)*180/pi/2,'b-')
title('mixed code/ late dm'); 


figure(14); clf; 
SP = subplot(1,2,1); cla; hold on; 
plot(squeeze(temp_choice_TR_l(:,1,:)),'bo')
plot(squeeze(temp_choice_TR_l(:,2,:)),'ro')


for iTRx = 1:nTR
    for istim = 1:2
        tempx = circ_mean(squeeze(aDecErr_sen_e(iTRx,istim,:))*2*pi/180)*180/pi/2; 
%         tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecErr_sen_e_fin(iTRx,istim) = tempx; 
        aDecSTD_sen_e_fin(iTRx,istim) = circ_std(squeeze(aDecErr_sen_e(iTRx,istim,:))*2*pi/180)*180/pi/2/sqrt(49); 
        
        tempx = circ_mean(squeeze(aDecErr_sen_l(iTRx,istim,:))*2*pi/180)*180/pi/2; 
%         tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecErr_sen_l_fin(iTRx,istim) = tempx; 
        aDecSTD_sen_l_fin(iTRx,istim) = circ_std(squeeze(aDecErr_sen_l(iTRx,istim,:))*2*pi/180)*180/pi/2/sqrt(49); 

        tempx = circ_mean(squeeze(aDecErr_mne_e(iTRx,istim,:))*2*pi/180)*180/pi/2; 
%         tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecErr_mne_e_fin(iTRx,istim) = tempx; 
        aDecSTD_mne_e_fin(iTRx,istim) = circ_std(squeeze(aDecErr_mne_e(iTRx,istim,:))*2*pi/180)*180/pi/2/sqrt(49); 

        tempx = circ_mean(squeeze(aDecErr_mne_l(iTRx,istim,:))*2*pi/180)*180/pi/2; 
%         tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecErr_mne_l_fin(iTRx,istim) = tempx; 
        aDecSTD_mne_l_fin(iTRx,istim) = circ_std(squeeze(aDecErr_TR_l(iTRx,istim,:))*2*pi/180)*180/pi/2/sqrt(49); 
        
        
        
        tempx = circ_mean(squeeze(aDecErr_TR_e(iTRx,istim,:))*pi/180)*180/pi; 
%         tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecErr_TR_e_fin(iTRx,istim) = tempx; 
        aDecSTD_TR_e_fin(iTRx,istim) = circ_std(squeeze(aDecErr_TR_e(iTRx,istim,:))*pi/180)*180/pi/sqrt(49); 

        tempx = circ_mean(squeeze(aDecErr_TR_l(iTRx,istim,:))*pi/180)*180/pi; 
%         tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecErr_TR_l_fin(iTRx,istim) = tempx; 
        aDecSTD_TR_l_fin(iTRx,istim) = circ_std(squeeze(aDecErr_TR_l(iTRx,istim,:))*pi/180)*180/pi/sqrt(49); 
    
    end
end


figure(1); clf; 
SP = subplot(1,2,1); cla; hold on; 
errorbar(1:nTR, aDecErr_sen_e_fin(:,1), aDecSTD_sen_e_fin(:,1), 'bo-','markerfacecolor','w'); 
errorbar(1:nTR, aDecErr_sen_e_fin(:,2), aDecSTD_sen_e_fin(:,2), 'ro-','markerfacecolor','w'); 

SP = subplot(1,2,2); cla; hold on; 
errorbar(1:nTR, aDecErr_sen_l_fin(:,1), aDecSTD_sen_l_fin(:,1), 'bo-','markerfacecolor','w'); 
errorbar(1:nTR, aDecErr_sen_l_fin(:,2), aDecSTD_sen_l_fin(:,2), 'ro-','markerfacecolor','w'); 


figure(2); clf; 
SP = subplot(1,2,1); cla; hold on; 
errorbar(1:nTR, aDecErr_mne_e_fin(:,1), aDecSTD_mne_e_fin(:,1), 'bo-','markerfacecolor','w'); 
errorbar(1:nTR, aDecErr_mne_e_fin(:,2), aDecSTD_mne_e_fin(:,2), 'ro-','markerfacecolor','w'); 

SP = subplot(1,2,2); cla; hold on; 
errorbar(1:nTR, aDecErr_mne_l_fin(:,1), aDecSTD_mne_l_fin(:,1), 'bo-','markerfacecolor','w'); 
errorbar(1:nTR, aDecErr_mne_l_fin(:,2), aDecSTD_mne_l_fin(:,2), 'ro-','markerfacecolor','w'); 



figure(3); clf; 
SP = subplot(1,2,1); cla; hold on; 
errorbar(1:nTR, aDecErr_TR_e_fin(:,1), aDecSTD_TR_e_fin(:,1), 'bo-','markerfacecolor','w'); 
errorbar(1:nTR, aDecErr_TR_e_fin(:,2), aDecSTD_TR_e_fin(:,2), 'ro-','markerfacecolor','w'); 

SP = subplot(1,2,2); cla; hold on; 
errorbar(1:nTR, aDecErr_TR_l_fin(:,1), aDecSTD_TR_l_fin(:,1), 'bo-','markerfacecolor','w'); 
errorbar(1:nTR, aDecErr_TR_l_fin(:,2), aDecSTD_TR_l_fin(:,2), 'ro-','markerfacecolor','w'); 




for iTRx = 1:nTR
    for istim = 1:length(stimcond)
        tempx = circ_mean(squeeze(aDecAngle_sen_e(iTRx,istim,:))*2*pi/180)*180/pi/2; 
        aDecAngle_sen_e_fin(iTRx,istim) = tempx; 
        
        tempx = circ_mean(squeeze(aDecAngle_sen_l(iTRx,istim,:))*2*pi/180)*180/pi/2; 
%         tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecAngle_sen_l_fin(iTRx,istim) = tempx; 
        
        tempx = circ_mean(squeeze(aDecAngle_mne_e(iTRx,istim,:))*2*pi/180)*180/pi/2; 
%         tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecAngle_mne_e_fin(iTRx,istim) = tempx; 
        
        tempx = circ_mean(squeeze(aDecAngle_mne_l(iTRx,istim,:))*2*pi/180)*180/pi/2; 
%         if sign(tempx) > 0
%             tempx(tempx<0) = tempx(tempx<0)+180; 
%         else
%             tempx(tempx>0) = tempx(tempx>0)-180; 
%         end
        aDecAngle_mne_l_fin(iTRx,istim) = tempx; 
    end
end




J = jet(24)*0.7 + 0.3; 

set(figure(1), 'position',[1 579 723 226]); clf; 
subplot(1,2,1); cla; hold on; 
for ii = 1:length(stimcond)
    plot(aDecAngle_sen_e_fin(:,ii), 'color', J(ii,:), 'linewidth',1.5)
end
xlabel(' Time (TR)' ); 
ylabel(' Decoded orientation (deg)' ); 
title('Early Dm'); 
ylim([-10 190]); xlim([0 15]); 
yticks(linspace(0,180,5)); 

subplot(1,2,2); cla; hold on; 
for ii = 1:length(stimcond)
    plot(aDecAngle_sen_l_fin(:,ii), 'color', J(ii,:), 'linewidth',1.5)
end
xlabel(' Time (TR)' ); 
ylabel(' Decoded orientation (deg)' ); 
title('Late Dm'); 
ylim([-10 190]); xlim([0 15]); 
yticks(linspace(0,180,5)); 


set(figure(2), 'position',[3 564 866 241]); clf; 
for iTR = 1:nTR
    subplot(2,7,iTR); cla; hold on; 
    for ii = 1:length(stimcond)
        plot(stimcond(ii), aDecAngle_sen_e_fin(iTR,ii), 'ko','markerfacecolor', J(ii,:),'markersize',4)
    end
    plot([0 180], [0 180],'k-'); 
    ylim([0 180]); xlim([0 180]); 
    title(['TR=' num2str(iTR)]); 
end


set(figure(3), 'position',[3 249 866 241]); clf; 
for iTR = 1:nTR
    subplot(2,7,iTR); cla; hold on; 
    for ii = 1:length(stimcond)
        plot(stimcond(ii), aDecAngle_sen_l_fin(iTR,ii), 'ko','markerfacecolor', J(ii,:),'markersize',4)
    end
    plot([0 180], [0 180],'k-'); 
    ylim([0 180]); xlim([0 180]); 
    title(['TR=' num2str(iTR)]); 
end





set(figure(11), 'position',[1 579 723 226]); clf; 
subplot(1,2,1); cla; hold on; 
for ii = 1:length(stimcond)
    plot(aDecAngle_mne_e_fin(:,ii), 'color', J(ii,:), 'linewidth',1.5)
end
xlabel(' Time (TR)' ); 
ylabel(' Decoded orientation (deg)' ); 
title('Early Dm'); 
ylim([-10 190]); xlim([0 15]); 
yticks(linspace(0,180,5)); 

subplot(1,2,2); cla; hold on; 
for ii = 1:length(stimcond)
    plot(aDecAngle_mne_l_fin(:,ii), 'color', J(ii,:), 'linewidth',1.5)
end
xlabel(' Time (TR)' ); 
ylabel(' Decoded orientation (deg)' ); 
title('Late Dm'); 
ylim([-10 190]); 
xlim([0 15]); 
yticks(linspace(0,180,5)); 


set(figure(12), 'position',[3 564 866 241]); clf; 
for iTR = 1:nTR
    subplot(2,7,iTR); cla; hold on; 
    for ii = 1:length(stimcond)
        plot(stimcond(ii), aDecAngle_mne_e_fin(iTR,ii), 'ko','markerfacecolor', J(ii,:),'markersize',4)
    end
    plot([0 180], [0 180],'k-'); 
    ylim([0 180]); xlim([0 180]); 
    title(['TR=' num2str(iTR)]); 
end


set(figure(13), 'position',[3 249 866 241]); clf; 
for iTR = 1:nTR
    subplot(2,7,iTR); cla; hold on; 
    for ii = 1:length(stimcond)
        plot(stimcond(ii), aDecAngle_mne_l_fin(iTR,ii), 'ko','markerfacecolor', J(ii,:),'markersize',4)
    end
    plot([0 180], [0 180],'k-'); 
    ylim([0 180]); xlim([0 180]); 
    title(['TR=' num2str(iTR)]); 
end




figure(); clf; hold on; 
for ii = 1:length(stimcond)
    plot(stimcond(ii),1, 'ko','markerfacecolor', J(ii,:));
end
xticks(linspace(0,180,5)); 
xlim([-5 185])
view(270,90)




%%
set(figure(20), 'position',[1 44 804 761]); clf; 
SP = subplot(2,2,1); cla; 
hold on; 
for ii = 1:length(stimcond)
    plot([0 15*cos(stimcond(ii)*2*pi/180)], [0, 15*sin(stimcond(ii)*2*pi/180)], 'k-', 'color', J(ii,:), 'linewidth',2.4)
end
xys = nan(nTR, length(stimcond),2); 
for iTR = 1:nTR
    for ii = 1:length(stimcond)
        xys(iTR,ii,:) = [iTR*cos(aDecAngle_sen_e_fin(iTR,ii)*2*pi/180), iTR*sin(aDecAngle_sen_e_fin(iTR,ii)*2*pi/180)]; 
    end
end
for ii = 1:length(stimcond)
    plot(squeeze(xys(:,ii,1)), squeeze(xys(:,ii,2)), 'ko-', 'color', J(ii,:), 'linewidth',1.5,'markerfacecolor', J(ii,:),'markersize',4)
end
title('Early Dm, sensory code (3-4TR)'); 

SP = subplot(2,2,2); cla; 
hold on; 
for ii = 1:length(stimcond)
    plot([0 15*cos(stimcond(ii)*2*pi/180)], [0, 15*sin(stimcond(ii)*2*pi/180)], 'k-', 'color', J(ii,:), 'linewidth',2.4)
end
xys = nan(nTR, length(stimcond),2); 
for iTR = 1:nTR
    for ii = 1:length(stimcond)
        xys(iTR,ii,:) = [iTR*cos(aDecAngle_sen_l_fin(iTR,ii)*2*pi/180), iTR*sin(aDecAngle_sen_l_fin(iTR,ii)*2*pi/180)]; 
    end
end
for ii = 1:length(stimcond)
    plot(squeeze(xys(:,ii,1)), squeeze(xys(:,ii,2)), 'ko-', 'color', J(ii,:), 'linewidth',1.5,'markerfacecolor', J(ii,:),'markersize',4)
end
title('Late Dm, sensory code (3-4TR)'); 

SP = subplot(2,2,3); cla; 
hold on; 
for ii = 1:length(stimcond)
    plot([0 15*cos(stimcond(ii)*2*pi/180)], [0, 15*sin(stimcond(ii)*2*pi/180)], 'k-', 'color', J(ii,:), 'linewidth',2.4)
end
xys = nan(nTR, length(stimcond),2); 
for iTR = 1:nTR
    for ii = 1:length(stimcond)
        xys(iTR,ii,:) = [iTR*cos(aDecAngle_mne_e_fin(iTR,ii)*2*pi/180), iTR*sin(aDecAngle_mne_e_fin(iTR,ii)*2*pi/180)]; 
    end
end
for ii = 1:length(stimcond)
    plot(squeeze(xys(:,ii,1)), squeeze(xys(:,ii,2)), 'ko-', 'color', J(ii,:), 'linewidth',1.5,'markerfacecolor', J(ii,:),'markersize',4)
end
title('Early Dm, mnemonic code (9-10TR)'); 

SP = subplot(2,2,4); cla; 
hold on; 
for ii = 1:length(stimcond)
    plot([0 15*cos(stimcond(ii)*2*pi/180)], [0, 15*sin(stimcond(ii)*2*pi/180)], 'k-', 'color', J(ii,:), 'linewidth',2.4)
end
xys = nan(nTR, length(stimcond),2); 
for iTR = 1:nTR
    for ii = 1:length(stimcond)
        xys(iTR,ii,:) = [iTR*cos(aDecAngle_mne_l_fin(iTR,ii)*2*pi/180), iTR*sin(aDecAngle_mne_l_fin(iTR,ii)*2*pi/180)]; 
    end
end
for ii = 1:length(stimcond)
    plot(squeeze(xys(:,ii,1)), squeeze(xys(:,ii,2)), 'ko-', 'color', J(ii,:), 'linewidth',1.5,'markerfacecolor', J(ii,:),'markersize',4)
end
title('Late Dm, mnemonic code (6-7TR)'); 






for iTRx = 1:nTR
    for istim = 1:length(stimcond)
        tempx = circ_mean(squeeze(aDecAngle_TR_e(iTRx,istim,:))*2*pi/180)*180/pi/2; 
        tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecAngle_TR_e_fin(iTRx,istim) = tempx; 
        
        tempx = circ_mean(squeeze(aDecAngle_TR_l(iTRx,istim,:))*2*pi/180)*180/pi/2; 
        tempx(tempx<0) = tempx(tempx<0)+180; 
        aDecAngle_TR_l_fin(iTRx,istim) = tempx; 
    end
end


set(figure(21), 'position',[1 44 804 761]); clf; 
SP = subplot(2,2,1); cla; 
hold on; 
for ii = 1:length(stimcond)
    plot([15*cos(stimcond(ii)*2*pi/180) 30*cos(stimcond(ii)*2*pi/180)], [15*sin(stimcond(ii)*2*pi/180), 30*sin(stimcond(ii)*2*pi/180)], 'k-', 'color', J(ii,:), 'linewidth',2.4)
end
xys = nan(nTR, length(stimcond),2); 
for iTR = 1:nTR
    for ii = 1:length(stimcond)
        xys(iTR,ii,:) = [(iTR+15)*cos(aDecAngle_TR_e_fin(iTR,ii)*2*pi/180), (iTR+15)*sin(aDecAngle_TR_e_fin(iTR,ii)*2*pi/180)]; 
    end
end
for ii = 1:length(stimcond)
    plot(squeeze(xys(:,ii,1)), squeeze(xys(:,ii,2)), 'ko-', 'color', J(ii,:), 'linewidth',1.5,'markerfacecolor', J(ii,:),'markersize',4)
end
title('Early Dm, same TRs'); 

SP = subplot(2,2,2); cla; 
hold on; 
for ii = 1:length(stimcond)
    plot([0 15*cos(stimcond(ii)*2*pi/180)], [0, 15*sin(stimcond(ii)*2*pi/180)], 'k-', 'color', J(ii,:), 'linewidth',2.4)
end
xys = nan(nTR, length(stimcond),2); 
for iTR = 1:nTR
    for ii = 1:length(stimcond)
        xys(iTR,ii,:) = [iTR*cos(aDecAngle_TR_l_fin(iTR,ii)*2*pi/180), iTR*sin(aDecAngle_TR_l_fin(iTR,ii)*2*pi/180)]; 
    end
end
for ii = 1:length(stimcond)
    plot(squeeze(xys(:,ii,1)), squeeze(xys(:,ii,2)), 'ko-', 'color', J(ii,:), 'linewidth',1.5,'markerfacecolor', J(ii,:),'markersize',4)
end
title('Late Dm, same TRs'); 




set(figure(22), 'position',[3 249 866 241]); clf; 
for iTR = 1:nTR
    subplot(2,7,iTR); cla; hold on; 
    for ii = 1:length(stimcond)
        plot(stimcond(ii), aDecAngle_TR_e_fin(iTR,ii), 'ko','markerfacecolor', J(ii,:),'markersize',4)
    end
    plot([0 180], [0 180],'k-'); 
    ylim([0 180]); xlim([0 180]); 
    title(['TR=' num2str(iTR)]); 
end



set(figure(23), 'position',[3 249 866 241]); clf; 
for iTR = 1:nTR
    subplot(2,7,iTR); cla; hold on; 
    for ii = 1:length(stimcond)
        plot(stimcond(ii), aDecAngle_TR_l_fin(iTR,ii), 'ko','markerfacecolor', J(ii,:),'markersize',4)
    end
    plot([0 180], [0 180],'k-'); 
    ylim([0 180]); xlim([0 180]); 
    title(['TR=' num2str(iTR)]); 
end

