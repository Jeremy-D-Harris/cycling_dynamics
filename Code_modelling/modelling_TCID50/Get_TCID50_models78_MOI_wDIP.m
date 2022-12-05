function TCID50_output = Get_TCID50_models78_MOI_wDIP(WT_moi,DI_moi,pars,y,eps,Ncells,which_model)

max_N_W = 1000; max_N_D = 1000;
N_W=max_N_W; N_D = max_N_D;

m_val = pars(1);
K_val = pars(2);

psiW = WT_moi; psiD = DI_moi;

num_WT_vector = 0:N_W;
prob_psiW_vector = poisspdf(num_WT_vector,psiW);
prob_above_N_W = 1 - poisscdf(N_W,psiW);
prob_psiW_vector(1,end) = prob_above_N_W;


num_DI_vector = 0:N_D;
prob_psiD_vector = poisspdf(num_DI_vector,psiD);
prob_above_N_D = 1 - poisscdf(N_D,psiD);
prob_psiD_vector(1,end) = prob_above_N_D;

% joint probability matrix, Poisson(i WT,j DI) = Poisson(i WT)*Poisson(j DI)
prob_psiW_psiD_matrix = transpose(prob_psiW_vector)*prob_psiD_vector;

phiW_ij_matrix = zeros(size(prob_psiW_psiD_matrix));
nu_ij_matrix = zeros(size(prob_psiW_psiD_matrix));

sum_W_ij = 0;

if which_model==7
    
    for i=1:length(num_WT_vector)
        
        nu_ij_matrix(i,:) = (num_WT_vector(i)*ones(size(num_DI_vector))+eps*num_DI_vector)./(K_val+(num_WT_vector(i)*ones(size(num_DI_vector))+eps*num_DI_vector));
        
        phiW_ij_matrix(i,1) = 1;
        phiW_ij_matrix(i,2:end) = y(1);
        
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(i,:).*nu_ij_matrix(i,:).*prob_psiW_psiD_matrix(i,:));
        
    end
    
else
    
    for i=1:length(num_WT_vector)
        
        nu_ij_matrix(i,:) = (num_WT_vector(i)*ones(size(num_DI_vector))+eps*num_DI_vector)./(K_val+(num_WT_vector(i)*ones(size(num_DI_vector))+eps*num_DI_vector));
        if i ==1
            phiW_ij_matrix(i,2:end) = (num_WT_vector(i)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end))./(num_WT_vector(i)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end)+y(1)*(1-eps)*num_DI_vector(1,2:end));
            %         phiW_ij_matrix(i,2:end) = (num_WT_vector(i)*ones(size(num_DI_vector(1,2:end))))./(num_WT_vector(i)*ones(size(num_DI_vector(1,2:end)))+y(1)*num_DI_vector(1,2:end));
            phiW_ij_matrix(1,1) = 1;
        else
            phiW_ij_matrix(i,:) = (num_WT_vector(i)*ones(size(num_DI_vector))+eps*num_DI_vector)./(num_WT_vector(i)*ones(size(num_DI_vector))+eps*num_DI_vector+y(1)*(1-eps)*num_DI_vector);
        end
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(i,:).*nu_ij_matrix(i,:).*prob_psiW_psiD_matrix(i,:));
        
    end
    
    
end

% y(2) = f, a scalar multiple to convert from WT to TCID50 measurement
TCID50_output=y(2)*m_val*Ncells*sum_W_ij;



