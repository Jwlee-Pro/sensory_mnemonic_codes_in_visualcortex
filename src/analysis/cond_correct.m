function err_corr2 = cond_correct(err_behav, err_m, ref_m, stim_m, stim_corr, ref_corr)

stimcond = 0:7.5:179; 
refcond = [-21, -4, 0, 4, 21]; 

lapse_crit = 30; 

err_corr = err_m ; 
if stim_corr == 'c'
    for istim = 1:length(stimcond)
        ind1 = find(stim_m==stimcond(istim) & ~isnan(err_behav) & abs(err_behav)<lapse_crit); 
        if length(ind1)>1
            temp = err_m(ind1) - circ_m(err_m(ind1)'); 
        else
            temp = err_m(ind1); 
        end
        temp(temp>90) = temp(temp>90) -180; 
        temp(temp<-90) = temp(temp<-90) +180;  
        err_corr(ind1) = temp; 
    end
end

err_corr2 = err_corr ;
if ref_corr == 'c'
    for ir = 1:length(refcond)
        ind1 = find(ref_m==refcond(ir) & ~isnan(err_behav) & abs(err_behav)<lapse_crit); 
        if length(ind1)>1
            temp = err_corr(ind1) - circ_m(err_corr(ind1)'); 
        else
            temp = err_corr(ind1); 
        end
        temp(temp>90) = temp(temp>90) -180; 
        temp(temp<-90) = temp(temp<-90) +180;  
        err_corr2(ind1) = temp; 
    end
end
