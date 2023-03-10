

clear all; close all; clc;



% Load sub/roi-list 
load('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/sub_list.mat')


%% Channel response 
nTR = 14; 
refs = [-21, -4, 0, 4, 21]; 
stimcond = 0:7.5:172.5; 

refChannel_e = nan(nTR,5,120,length(sub_list)); 
refChannel_l = nan(nTR,5,120,length(sub_list)); 

fidelity_decx = cell(1,length(sub_list)); 
fidelityest_decx = cell(1,length(sub_list)); 
for isub = 1:length(sub_list)
    isub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded/VC_sub-' sub_list(isub,:) '.mat'])
    
    [nx,ny,nTrials,nChannels] = size(Early); 
    s_pre = linspace(0, 2*pi, nChannels +1)'; s_pre(end) = []; 
    
    DecMat = nan(size(Early)); 
    DecMat(:,:,Timing==1,:) = Early(:,:,Timing==1,:);
    DecMat(:,:,Timing==2,:) = Late(:,:,Timing==2,:);
    
    Decoded_result = cell(1,14); 
    for iTR1 = 1:14
        Decoded_result{iTR1}.est = nan(14,nTrials);
        Decoded_result{iTR1}.channelResp = squeeze(DecMat(iTR1,1:14,:,:)); 
        for iTR2 = 1:14
            pop_vec = squeeze(DecMat(iTR1,iTR2,:,:))*exp(1i*s_pre);
            Decoded_result{iTR1}.est(iTR2,:) = mod(angle(pop_vec)/pi*90, 180); 
        end
    end
    
    stimulus = stim; 
    
    % Shift to stimulus 
    for iTR = 1:nTR
        lkhCurr = nan(length(stimulus) , 120); 

        encodingTrained1 = squeeze((Decoded_result{9}.channelResp(iTR,:,:) + Decoded_result{10}.channelResp(iTR,:,:))/2) ; 
        encodingTrained2 = squeeze((Decoded_result{6}.channelResp(iTR,:,:) + Decoded_result{7}.channelResp(iTR,:,:))/2) ; 
        encodingTrained = nan(size(encodingTrained1)); 
        encodingTrained(Timing==1,:) = encodingTrained1(Timing==1,:); 
        encodingTrained(Timing==2,:) = encodingTrained2(Timing==2,:); 
        for iT = 1:length(stimulus) 
            lkhCurr(iT,:) = circshift(squeeze(encodingTrained(iT,:)), (90-stimulus(iT))/1.5); 
        end
        
        for ir = 1:length(refs)
            refChannel_e(iTR,ir,:,isub) = nanmean(lkhCurr(ref==refs(ir) & Timing==1,:)); 
            refChannel_l(iTR,ir,:,isub) = nanmean(lkhCurr(ref==refs(ir) & Timing==2,:)); 
        end
    end
    
end


%% figure plotting
xvals = linspace(-90,90,120); 
refcond = [19 47 120; 147 189 234; 190 190 190; 238 139 129; 120 14 32]/255; 

set(figure(2),'position',[1 491 1343 314]); clf; 
for iTR = 1:nTR
    SP = subplot(2,7,iTR); cla; hold on; 
    for ir = 1:5
        plot(xvals, nanmean(squeeze(refChannel_e(iTR,ir,:,:))'*2*pi/180)*180/pi/2, 'k-','color',refcond(ir,:),'linewidth',1.3); 
    end
    
    xlabel('Reference (deg)'); 
    ylabel('STD (deg)'); 
    
    title(['TR=' num2str(iTR)]); 
end
figure()
hold on; 
for ir = 1:5
    plot(ir, 0, 'ko','markerfacecolor', refcond(ir,:))
end

set(figure(3),'position',[1 103 1343 314]); clf; 
for iTR = 1:nTR
    SP = subplot(2,7,iTR); cla; hold on; 
    for ir = 1:5
        plot(xvals, nanmean(squeeze(refChannel_l(iTR,ir,:,:))'*2*pi/180)*180/pi/2, 'k-','color',refcond(ir,:),'linewidth',1.3); 
    end
    
    xlabel('Reference (deg)'); 
    ylabel('STD (deg)'); 
    
    title(['TR=' num2str(iTR)]); 
end

%% Re-centering & calculate fidelity
recenter_e = nan(nTR,5,120); 
recenter_l = nan(nTR,5,120); 
fidelity_e = nan(nTR,5); 
fidelity_l = nan(nTR,5); 
for iTR = 1:nTR
    for ir = 1:5
        pop_vec = nanmean(squeeze(refChannel_e(iTR,ir,:,:))')*exp(1i*s_pre);
        delta_value = mod(angle(pop_vec)/pi*90, 180);
        temp = circshift(nanmean(squeeze(refChannel_e(iTR,ir,:,:))'), round((90-delta_value)/1.5));
        recenter_e(iTR,ir,:) = temp; 
        fidelity_e(iTR,ir) = fidelity(temp', 90*2*pi/180);
        
        pop_vec = nanmean(squeeze(refChannel_l(iTR,ir,:,:))')*exp(1i*s_pre);
        delta_value = mod(angle(pop_vec)/pi*90, 180);
        temp = circshift(nanmean(squeeze(refChannel_l(iTR,ir,:,:))'), round((90-delta_value)/1.5));
        recenter_l(iTR,ir,:) = temp; 
        fidelity_l(iTR,ir) = fidelity(temp', 90*2*pi/180);
    end
end


set(figure(4),'position',[1 103 1343 314]); clf; 
for iTR = 1:nTR
    SP = subplot(2,7,iTR); cla; hold on; 
    plot(refs, fidelity_e(iTR,:),'bo-'); 
    plot(refs, fidelity_l(iTR,:),'ro-'); 
    xlabel('Reference (deg)'); 
    ylabel('Fidelity (deg)'); 
    
    title(['TR=' num2str(iTR)]); 
end
