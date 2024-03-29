clear all; close all; clc;

% Normalizing each time point with baseline channel response (could be replaced with BOLD)
normalize = 1; 

% Load sub/roi-list 
load('/Volumes/ROOT/CSNL_temp/JWL/Analysis_2021DecSummary/sub_list.mat')

% Time points to average (5~13TRs)
rangeT = 5:13; 

% Decoding channel responses
for isub = 1:length(sub_list) 
    fprintf([num2str(isub) '\n'])
    load(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded/VC_sub-' sub_list(isub,:) '.mat'])
    
    [nx,ny,nTrials,nChannels] = size(Early); 
    s_pre = linspace(0, 2*pi, nChannels +1)'; s_pre(end) = []; 
    
    % Decoded channel response trained with "early decision's delay period" 
    clear Decoded_result
    DecMat = squeeze(nanmean(Early(9:10,rangeT,:,:),1)); 
    
    if normalize == 1
        for ii = 1:size(DecMat,1)
            temp = DecMat(ii,Timing==1,:); 
            DecMat(ii,Timing==1,:) = (temp - nanmean(temp(:)))/nanstd(temp(:)); 
            temp = DecMat(ii,Timing==2,:); 
            DecMat(ii,Timing==2,:) = (temp - nanmean(temp(:)))/nanstd(temp(:)); 
        end
    end
    
    % Decode channel response 
    DecMat_avg = squeeze(nanmean(DecMat,1));
    pop_vec = DecMat_avg*exp(1i*s_pre);
    Decoded_result.est = mod(angle(pop_vec)/pi*90, 180); 
    
    % Save files
    stimulus = stim; 
    timing = Timing; 
    response = esti; 
    
    % Split first & second half 
    pop_vec = squeeze(nanmean(DecMat(1:5,:,:),1))*exp(1i*s_pre);
    Decoded_result.est_firstHalf = mod(angle(pop_vec)/pi*90, 180); 
    
    pop_vec = squeeze(nanmean(DecMat(6:9,:,:),1))*exp(1i*s_pre);
    Decoded_result.est_secondHalf = mod(angle(pop_vec)/pi*90, 180); 
    
    % Save data 
    if normalize == 1
        save(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_earlydelay/VC_sub-' sub_list(isub,:) '_meannorm_dec.mat'],'stimulus','timing','response','ref','choice','Decoded_result')
    else
        save(['/Volumes/ROOT/CSNL_temp/JWL/sensory_mnemonic_codes_in_visualcortex/data/decoded_earlydelay/VC_sub-' sub_list(isub,:) '_mean_dec.mat'],'stimulus','timing','response','ref','choice','Decoded_result')
    end
end
