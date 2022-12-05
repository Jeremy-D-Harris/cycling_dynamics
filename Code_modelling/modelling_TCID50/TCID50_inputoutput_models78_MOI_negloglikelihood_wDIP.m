function f = TCID50_inputoutput_models78_MOI_negloglikelihood_wDIP(y, newdata_TCID50, TCID50_noise, WT_mois, DI_mois, pars, eps,Ncells,which_model)

% n_mios = length(WT_mois(:,1));
n_replicates = length(WT_mois(1,:));
n_treatments = length(WT_mois(:,1));
% m = pars(1); K = pars(2);


logL = 0;
for nn = 1:n_replicates
    
    for kk=1:n_treatments
        
        these_data = newdata_TCID50(kk,nn);
        this_WT_moi = WT_mois(kk,nn); this_DI_moi = DI_mois(kk,nn);
        
        if isnan(these_data)~=1
            if isnan(this_WT_moi)~=1
                
                TCID50_model = Get_TCID50_models78_MOI_wDIP(this_WT_moi,this_DI_moi, pars, y, eps,Ncells,which_model);
                
                logL = logL + log(normpdf(log10(these_data), log10(TCID50_model), TCID50_noise(kk)));
            end
        end
        
    end
    
    
end
f = -logL;
