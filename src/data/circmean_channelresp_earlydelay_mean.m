

clear all; close all; clc;



% Load sub/roi-list 
load('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/sub_list.mat')


for isub = 1:length(sub_list) 
    isub
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded/VC_sub-' sub_list(isub,:) '.mat'])
    
    [nx,ny,nTrials,nChannels] = size(Early); 
    s_pre = linspace(0, 2*pi, nChannels +1)'; s_pre(end) = []; 
    
    % early decision (delay period) 
    clear Decoded_result
    DecMat = squeeze(nanmean(Early(9:10,1:14,:,:),1)); 
    DecMat_avg = squeeze(nanmean(DecMat(5:13,:,:),1)); 
    
    pop_vec = DecMat_avg*exp(1i*s_pre);
    Decoded_result.est = mod(angle(pop_vec)/pi*90, 180); 
    
    % Save files
    stimulus = stim; 
    timing = Timing; 
    response = esti; 

    pop_vec = squeeze(nanmean(DecMat(5:9,:,:),1))*exp(1i*s_pre);
    Decoded_result.est_firstHalf = mod(angle(pop_vec)/pi*90, 180); 
    
    pop_vec = squeeze(nanmean(DecMat(10:13,:,:),1))*exp(1i*s_pre);
    Decoded_result.est_secondHalf = mod(angle(pop_vec)/pi*90, 180); 
    
    save(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_earlydelay/VC_sub-' sub_list(isub,:) '_mean_dec.mat'],'stimulus','timing','response','ref','choice','Decoded_result')
end
