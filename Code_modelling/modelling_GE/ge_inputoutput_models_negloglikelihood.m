function f = ge_inputoutput_models_negloglikelihood(x, data_hpi_18hrs, effective_moi, Ncells, GE_noise, which_model)

n_reps = length(effective_moi(1,:));
n_treats = length(effective_moi(:,1));
logL = 0;
for ii = 1:n_reps
    GE_model = Get_GE_models(effective_moi(:,ii), x, Ncells, which_model);
    these_data = data_hpi_18hrs(:,ii);
    for kk = 1:n_treats
        logL = logL + log(normpdf(log10(these_data(kk)), log10(GE_model(kk)), GE_noise(kk)));
    end
end
f = -logL;