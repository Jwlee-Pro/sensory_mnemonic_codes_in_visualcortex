

clear all; close all; clc;



% Load sub/roi-list 
load('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/sub_list.mat')


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
        for iTR2 = 1:14
            pop_vec = squeeze(DecMat(iTR1,iTR2,:,:))*exp(1i*s_pre);
            Decoded_result{iTR1}.est(iTR2,:) = mod(angle(pop_vec)/pi*90, 180); 
        end
    end
    
    % Save files
    stimulus = stim; 
    timing = Timing; 
    response = esti; 
    save(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_estimated/VC_sub-' sub_list(isub,:) '_dec.mat'],'stimulus','timing','response','ref','choice','Decoded_result')
end
