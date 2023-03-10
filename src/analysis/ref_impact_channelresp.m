

clear all; close all; clc;



% Load sub/roi-list 
load('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/sub_list.mat')


%% Channel response 
nTR = 14; 
refs = [-21, -4, 0, 4, 21]; 
stimcond = 0:7.5:172.5; 

refChannel_e = nan(nTR,120,length(sub_list)); 
refChannel_l = nan(nTR,120,length(sub_list)); 

fidelity_decx = cell(1,length(sub_list)); 
fidelityest_decx = cell(1,length(sub_list)); 
for isub = 1:length(sub_list)
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
            refChannel_e(iTR,:,isub) = nanmean(lkhCurr(ref==refs(ir) & Timing==1,:)); 
            refChannel_l(iTR,:,isub) = nanmean(lkhCurr(ref==refs(ir) & Timing==2,:)); 
        end
    end
    
end