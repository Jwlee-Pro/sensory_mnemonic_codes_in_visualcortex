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


% plotting timecourse of estimates
% boostrap 95% confidence interval
nB=10000;
matBE=[];
matBL=[];
for iB=1:nB
    tv= floor(rand(1,50)*50)+1;
    matBE(:,:,iB)=mean(matfMRI.ori.mean_e(:,:,tv),3);
    matBL(:,:,iB)=mean(matfMRI.ori.mean_l(:,:,tv),3);
end
matBE_mean=mean(matfMRI.ori.mean_e,3);
matBE_ste=std(matBE,[],3);
matBL_mean=mean(matfMRI.ori.mean_l,3);
matBL_ste=std(matBL,[],3);

matBEL_mean=(matBE_mean+matBL_mean)/2;
matBEL_ste=sqrt((matBE_ste.^2+matBL_ste.^2)/2);

cmap_ori=hsv(24)*.75;
figure(201); clf;


axixT=1:2:1+2*13;
tValid=5:13;
tL=find(stimcond>90);
stimcond_a=stimcond;
stimcond_a(tL)=stimcond(tL)-180;

SP=subplot(1,3,1)
hold on;
for iOri=1:length(stimcond)
     if(rem(iOri,3)==1)
        errorbar(axisT(tValid)',matBE_mean(tValid,iOri)+stimcond_a(iOri),matBE_ste(tValid,iOri),'LineWidth',2,'Color',cmap_ori(iOri,:));
    else
        plot(axisT(tValid)',matBE_mean(tValid,iOri)+stimcond_a(iOri),'-','LineWidth',2,'Color',cmap_ori(iOri,:));
   
    end
end
xlim([5 28]); ylim([-100 100]);
set(line([5 28],[0 0]),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]+45),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]-45),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]+90),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]-90),'LineStyle','--','Color',[0 0 0]+.5)

SP=subplot(1,3,2)
hold on;
for iOri=1:length(stimcond)
     if(rem(iOri,3)==1)
        errorbar(axisT(tValid)',matBL_mean(tValid,iOri)+stimcond_a(iOri),matBL_ste(tValid,iOri),'LineWidth',2,'Color',cmap_ori(iOri,:));
    else
        plot(axisT(tValid)',matBL_mean(tValid,iOri)+stimcond_a(iOri),'-','LineWidth',2,'Color',cmap_ori(iOri,:));
   
    end
end
xlim([5 28]); ylim([-100 100]);
set(line([5 28],[0 0]),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]+45),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]-45),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]+90),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]-90),'LineStyle','--','Color',[0 0 0]+.5)

SP=subplot(1,3,3); cla;
hold on;
for iOri=1:length(stimcond)
    if(rem(iOri,3)==1)
        errorbar(axisT(tValid)',matBEL_mean(tValid,iOri)+stimcond_a(iOri),matBEL_ste(tValid,iOri),'LineWidth',2,'Color',cmap_ori(iOri,:));
    else
        plot(axisT(tValid)',matBEL_mean(tValid,iOri)+stimcond_a(iOri),'-','LineWidth',2,'Color',cmap_ori(iOri,:));
   
    end
end
xlim([5 28]); ylim([-100 100]);
set(line([5 28],[0 0]),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]+45),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]-45),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]+90),'LineStyle','--','Color',[0 0 0]+.5)
set(line([5 28],[0 0]-90),'LineStyle','--','Color',[0 0 0]+.5)









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
% temp1 = matfMRI.ori.mean(5:9,:,:);
% temp2 = matfMRI.ori.mean(10:13,:,:);

temp1 = matfMRI.ori.mean(5:7,:,:);
temp2 = matfMRI.ori.mean(6:8,:,:);
temp3 = matfMRI.ori.mean(7:9,:,:);
temp4 = matfMRI.ori.mean(8:10,:,:);
temp5 = matfMRI.ori.mean(9:11,:,:);
temp6 = matfMRI.ori.mean(10:12,:,:);
temp7 = matfMRI.ori.mean(11:13,:,:);

e=[];
for isub = 1:length(sub_list)
    e{1}(:,isub) = circ_m(temp1(:,:,isub)); 
    e{2}(:,isub) = circ_m(temp2(:,:,isub)); 
    e{3}(:,isub) = circ_m(temp3(:,:,isub)); 

    e{4}(:,isub) = circ_m(temp4(:,:,isub)); 
    e{5}(:,isub) = circ_m(temp5(:,:,isub)); 
    e{6}(:,isub) = circ_m(temp6(:,:,isub)); 

    e{7}(:,isub) = circ_m(temp7(:,:,isub)); 

end

SP = subplot(1,2,1); cla; hold on; 
cmap=ones(7,3)-pink(7);
rb=[];
for iE=1:7
    ee=e{iE};
    errorbar(stimHalf,circ_m(ori_cardinal(ee)'),circ_s(ori_cardinal(ee)')/sqrt(length(sub_list)-1) ,'o-','capsize',1,'color',cmap(iE,:));

    pre = ori_cardinal(ee); 

    rb(iE,:)=nanmean([pre(8:12,:); -pre(2:6,:)],1);
end
% errorbar(stimHalf,circ_m(ori_cardinal(e1)'),circ_s(ori_cardinal(e1)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1,'color',[0 0 0]+0.7);
% errorbar(stimHalf,circ_m(ori_cardinal(e2)'),circ_s(ori_cardinal(e2)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1,'color',[0 0 0]+0.4); 
% errorbar(stimHalf,circ_m(ori_cardinal(e3)'),circ_s(ori_cardinal(e3)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1,'color',[0 0 0]); 
% plot([-45 45],[0 0]+circ_m(circ_m(ori_cardinal(e1)')'),'k--'); 
xlabel('orientation from cardinal (deg)'); ylabel('bias (deg)');
xlim([-45 45]); xticks(linspace(-45,45,3));
% legend('5-7TR','8-10TR','11-13TR');

figure(202); clf; hold on;
set(line([5 28],[0 0]),'LineWidth',2,'Color',[0 0 0]+.75)
axisT_bin=[6:1:12]*2-1;
eb=std(rb,[],2)/sqrt(49);
errorbar(axisT_bin, mean(rb,2)',eb,'Color',[0 0 0]);
plot(axisT_bin, mean(rb,2)','wo','MarkerSize',10,'MarkerFaceColor',[0 0 0]+.5,'LineWidth',2);
plot(axisT_bin(7), mean(rb(7,:)),'wo','MarkerSize',10,'MarkerFaceColor',[0 0 0],'LineWidth',2);
xlim([5 28]); ylim([-7 7])
[H,P]=ttest(rb(7,:));
text(axisT_bin(7)+1, mean(rb(7,:))+1, ['p = ',num2str(P)])



% early stage (TR5~9) vs. late stage (TR10~13)
% temp1 = matfMRI.ori.std(5:9,:,:);
% temp2 = matfMRI.ori.std(10:13,:,:);

temp1 = matfMRI.ori.std(5:7,:,:);
temp2 = matfMRI.ori.std(8:10,:,:);
temp3 = matfMRI.ori.std(11:13,:,:);
for isub = 1:length(sub_list)
    k1(:,isub) = circ_m(temp1(:,:,isub)); 
    k2(:,isub) = circ_m(temp2(:,:,isub)); 
    k3(:,isub) = circ_m(temp3(:,:,isub));
end

SP = subplot(1,2,2); cla; hold on; 
errorbar(stimHalf,circ_m(ori_cardinal(k1)'),circ_s(ori_cardinal(k1)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1,'color',[0 0 0]+0.7);
errorbar(stimHalf,circ_m(ori_cardinal(k2)'),circ_s(ori_cardinal(k2)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1,'color',[0 0 0]+0.4); 
errorbar(stimHalf,circ_m(ori_cardinal(k3)'),circ_s(ori_cardinal(k3)')/sqrt(length(sub_list)-1) ,'ko-','capsize',1,'color',[0 0 0]);
% plot([-45 45],[0 0]+circ_m(circ_m(ori_cardinal(e1)')'),'k--'); 
xlabel('orientation from cardinal (deg)'); ylabel('variability (deg)');
xlim([-45 45]); xticks(linspace(-45,45,3));




set(figure(15),'position',[1 129 398 176]); clf; 
SP = subplot(1,2,1); cla; hold on; 
pre = ori_cardinal(e{1}); 
peri= ori_cardinal(e{5}); 
post= ori_cardinal(e{7}); 
% plot(circ_m([pre(8:9,:); -pre(5:6,:)]), circ_m([post(8:9,:); -post(5:6,:)]),'ko');
plot(nanmean([pre(8:12,:); -pre(2:6,:)],1), nanmean([peri(8:12,:); -peri(2:6,:)],1),'bo'); hold on;
plot(nanmean([pre(8:12,:); -pre(2:6,:)],1), nanmean([post(8:12,:); -post(2:6,:)],1),'ro'); 
[h,p,ci,stats] = ttest(nanmean([pre(8:12,:); -pre(2:6,:)],1)- nanmean([post(8:12,:); -post(2:6,:)],1))

% plot(nanmean([pre(8:10,:); -pre(4:6,:)],1), nanmean([post(8:10,:); -post(4:6,:)],1),'ko');
% [h,p,ci,stats] = ttest(nanmean([pre(8:10,:); -pre(4:6,:)],1)- nanmean([post(8:10,:); -post(4:6,:)],1))


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
xval = [axisT(5:13) flip(axisT(5:13))]; 
temp = near_e; 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 

meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 
yval = [-stdval(5:13)+meanval(5:13) +flip(stdval(5:13)+meanval(5:13))] ; 
patch(xval, yval, nearfar_color(1,:),'facealpha',0.2,'edgecolor','none'); 
plot(axisT(5:13),meanval(5:13),'color',nearfar_color(1,:),'linewidth',1.5)

temp = far_e; 
stdval  = circ_std(temp'*2*pi/180)*180/pi/2/sqrt(length(sub_list)-1); 
meanval = circ_mean(temp'*2*pi/180)*180/pi/2; 

yval = [-stdval(5:13)+meanval(5:13) +flip(stdval(5:13)+meanval(5:13))] ; 
patch(xval, yval, nearfar_color(2,:),'facealpha',0.2,'edgecolor','none'); 
plot(axisT(5:13),meanval(5:13),'color',nearfar_color(2,:),'linewidth',1.5)

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


% Individual analysis
div_fmri = nanmean(near_e(6:11,:),1) ./ nanmean(far_e(6:11,:),1); 
div_beh  = nanmean(matBeh.ref.std(2:4,:)) ./ nanmean(matBeh.ref.std([1 5],:)); 

figure(19); clf; 
SP = subplot(1,1,1); cla; hold on; 
plot(div_beh, div_fmri, 'ko'); 
[coef, pval] = corr(div_beh', div_fmri')
xlabel('DIV behavior'); ylabel('DIV fmri (early Dm)');
title(['r=' num2str(coef) ', p=' num2str(pval) ]); 




% Near vs Far
near_e
vectP=[];
for iBin=1:7

    nn=mean(near_e(iBin+[4 5 6],:));
    BNE_mean(iBin)=circ_m(nn')';
    BNE_ste(iBin)=circ_s(nn')'/sqrt(49);

    ff=mean(far_e(iBin+[4 5 6],:));
    BFE_mean(iBin)=circ_m(ff')';
    BFE_ste(iBin)=circ_s(ff')'/sqrt(49);

    [H,P]=ttest(nn-ff);
    vectP(iBin)=P;

end

figure(203); clf; hold on;
% set(line([5 28],[0 0]),'LineWidth',2,'Color',[0 0 0]+.75)
axisT_bin=[6:1:12]*2-1;

errorbar(axisT_bin, BFE_mean,BFE_ste,'Color',[0 0 0]+.5);
plot(axisT_bin, BFE_mean,'wo','MarkerSize',10,'MarkerFaceColor',[0 0 0]+.5,'LineWidth',2);

errorbar(axisT_bin, BNE_mean,BNE_ste,'Color',[0 0 0]);
plot(axisT_bin, BNE_mean,'wo','MarkerSize',10,'MarkerFaceColor',[0 0 0],'LineWidth',2);

xlim([5 28]); ylim([-7 7])

plot(axisT_bin(7), mean(rb(7,:)),'wo','MarkerSize',10,'MarkerFaceColor',[0 0 0],'LineWidth',2);

[H,P]=ttest(rb(7,:));
text(axisT_bin(7)+1, mean(rb(7,:))+1, ['p = ',num2str(P)])




% sliding bin choice plotting
figure(204); clf; hold on;
eCCW = squeeze(matfMRI.choice.mean_e(:,1,:)); 
eCW = squeeze(matfMRI.choice.mean_e(:,2,:)); 

lCCW = squeeze(matfMRI.choice.mean_l(:,1,:)); 
lCW = squeeze(matfMRI.choice.mean_l(:,2,:)); 


vectP=[];
for iBin=1:7

    ccw=mean(eCCW(iBin+[4 5 6],:));
    BCE_mean(iBin)=circ_m(ccw')';
    BCE_ste(iBin)=circ_s(ccw')'/sqrt(49);

    cw=mean(eCW(iBin+[4 5 6],:));
    BWE_mean(iBin)=circ_m(cw')';
    BWE_ste(iBin)=circ_s(cw')'/sqrt(49);

    [H,P]=ttest(ccw-cw);
    vectP(iBin)=P;

end

SP=subplot(1,2,1); cla; hold on;

set(line([5 28],[0 0]),'LineWidth',2,'Color',[0 0 0]+.75)

errorbar(axisT_bin, BCE_mean,BCE_ste,'Color','m');
plot(axisT_bin, BCE_mean,'wo','MarkerSize',10,'MarkerFaceColor','m','LineWidth',2);
tL=find(vectP<.05);
plot(axisT_bin(tL), BCE_mean(tL),'wo','MarkerSize',10,'MarkerFaceColor','r','LineWidth',2);

errorbar(axisT_bin, BWE_mean,BWE_ste,'Color','c');
plot(axisT_bin, BWE_mean,'wo','MarkerSize',10,'MarkerFaceColor','c','LineWidth',2);
tL=find(vectP<.05);
plot(axisT_bin(tL), BWE_mean(tL),'wo','MarkerSize',10,'MarkerFaceColor','b','LineWidth',2);

xlim([5 28]); ylim([-1 1]*12)


% Late trial

vectP=[];
for iBin=1:7

    ccw=mean(lCCW(iBin+[4 5 6],:));
    BCL_mean(iBin)=circ_m(ccw')';
    BCL_ste(iBin)=circ_s(ccw')'/sqrt(49);

    cw=mean(lCW(iBin+[4 5 6],:));
    BWL_mean(iBin)=circ_m(cw')';
    BWL_ste(iBin)=circ_s(cw')'/sqrt(49);

    [H,P]=ttest(ccw-cw);
    vectP(iBin)=P;

end

SP=subplot(1,2,2); cla; hold on;

set(line([5 28],[0 0]),'LineWidth',2,'Color',[0 0 0]+.75)

errorbar(axisT_bin, BCL_mean,BCL_ste,'Color','m');
plot(axisT_bin, BCL_mean,'wo','MarkerSize',10,'MarkerFaceColor','m','LineWidth',2);
tL=find(vectP<.05);
plot(axisT_bin(tL), BCL_mean(tL),'wo','MarkerSize',10,'MarkerFaceColor','r','LineWidth',2);

errorbar(axisT_bin, BWL_mean,BWL_ste,'Color','c');
plot(axisT_bin, BWL_mean,'wo','MarkerSize',10,'MarkerFaceColor','c','LineWidth',2);
tL=find(vectP<.05);
plot(axisT_bin(tL), BWL_mean(tL),'wo','MarkerSize',10,'MarkerFaceColor','b','LineWidth',2);

xlim([5 28]); ylim([-1 1]*12)