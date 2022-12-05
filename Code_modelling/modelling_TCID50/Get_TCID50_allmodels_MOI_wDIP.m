function TCID50_output = Get_TCID50_allmodels_MOI_wDIP(WT_moi,DI_moi,pars,y,eps,Ncells,which_model)

max_N_W = 1000; max_N_D = 1000;
N_W=max_N_W; N_D = max_N_D;

if ismember(which_model,[1,2])
    
    nu_val = pars(1);
    
elseif ismember(which_model,[3,4])
    
    m_val = pars(1);
    
elseif ismember(which_model,[5,6])
    
    m_val = pars(1);
    K_val = pars(2);
    n_val = pars(3);
    
else
    
    m_val = pars(1);
    K_val = pars(2);
    
end

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

%% go through models 1-8
if which_model==1
    
    %% models 1-2: constant viral yield
    
    for ii=1:length(num_WT_vector)
        
        if ii==1
            nu_ij_matrix(ii,:)=0;
            
            ind_effectiveDI = find(num_DI_vector>1/eps);
            nu_ij_matrix(ii,ind_effectiveDI) = nu_val;
            
        else
            nu_ij_matrix(ii,:) = nu_val;
        end
        
        phiW_ij_matrix(ii,1) = 1;
        phiW_ij_matrix(ii,2:end) = y(1);
        
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(ii,:).*nu_ij_matrix(ii,:).*prob_psiW_psiD_matrix(ii,:));
        
    end
    
elseif which_model==2
    
    for ii=1:length(num_WT_vector)
        
        if ii==1
            nu_ij_matrix(ii,:)=0;
            
            ind_effectiveDI = find(num_DI_vector>1/eps);
            nu_ij_matrix(ii,ind_effectiveDI) = nu_val;
            
        else
            nu_ij_matrix(ii,:) = nu_val;
        end
        
        if ii ==1
            phiW_ij_matrix(ii,2:end) = (num_WT_vector(ii)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end))./(num_WT_vector(ii)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end)+y(1)*(1-eps)*num_DI_vector(1,2:end));
            %         phiW_ij_matrix(i,2:end) = (num_WT_vector(i)*ones(size(num_DI_vector(1,2:end))))./(num_WT_vector(i)*ones(size(num_DI_vector(1,2:end)))+y(1)*num_DI_vector(1,2:end));
            phiW_ij_matrix(1,1) = 1;
        else
            phiW_ij_matrix(ii,:) = (num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector)./(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector+y(1)*(1-eps)*num_DI_vector);
        end
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(ii,:).*nu_ij_matrix(ii,:).*prob_psiW_psiD_matrix(ii,:));
        
    end
    
    
elseif which_model==3
    
    %% models 3-4: Linear viral yield
    
    for ii=1:length(num_WT_vector)
        
        nu_ij_matrix(ii,:) = m_val*(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector);
        
        phiW_ij_matrix(ii,1) = 1;
        phiW_ij_matrix(ii,2:end) = y(1);
        
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(ii,:).*nu_ij_matrix(ii,:).*prob_psiW_psiD_matrix(ii,:));
        
    end
    
elseif which_model==4
    
    for ii=1:length(num_WT_vector)
        
        nu_ij_matrix(ii,:) = m_val*(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector);
        
        if ii ==1
            phiW_ij_matrix(ii,2:end) = (num_WT_vector(ii)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end))./(num_WT_vector(ii)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end)+y(1)*(1-eps)*num_DI_vector(1,2:end));
            %         phiW_ij_matrix(i,2:end) = (num_WT_vector(i)*ones(size(num_DI_vector(1,2:end))))./(num_WT_vector(i)*ones(size(num_DI_vector(1,2:end)))+y(1)*num_DI_vector(1,2:end));
            phiW_ij_matrix(1,1) = 1;
        else
            phiW_ij_matrix(ii,:) = (num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector)./(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector+y(1)*(1-eps)*num_DI_vector);
        end
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(ii,:).*nu_ij_matrix(ii,:).*prob_psiW_psiD_matrix(ii,:));
        
    end
    
    
elseif which_model==5
    
    %% models 5-6: Hill viral yield
    
    for ii=1:length(num_WT_vector)
        
        nu_ij_matrix(ii,:) = m_val*(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector).^n_val/(K_val+(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector).^n_val);
        
        phiW_ij_matrix(ii,1) = 1;
        phiW_ij_matrix(ii,2:end) = y(1);
        
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(ii,:).*nu_ij_matrix(ii,:).*prob_psiW_psiD_matrix(ii,:));
        
    end
    
elseif which_model==6
    
    for ii=1:length(num_WT_vector)
        
        nu_ij_matrix(ii,:) = m_val*(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector).^n_val/(K_val+(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector).^n_val);
        
        if ii ==1
            phiW_ij_matrix(ii,2:end) = (num_WT_vector(ii)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end))./(num_WT_vector(ii)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end)+y(1)*(1-eps)*num_DI_vector(1,2:end));
            %         phiW_ij_matrix(i,2:end) = (num_WT_vector(i)*ones(size(num_DI_vector(1,2:end))))./(num_WT_vector(i)*ones(size(num_DI_vector(1,2:end)))+y(1)*num_DI_vector(1,2:end));
            phiW_ij_matrix(1,1) = 1;
        else
            phiW_ij_matrix(ii,:) = (num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector)./(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector+y(1)*(1-eps)*num_DI_vector);
        end
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(ii,:).*nu_ij_matrix(ii,:).*prob_psiW_psiD_matrix(ii,:));
        
    end
    
    
elseif which_model==7
    %% models 7-8: MM viral yield
    for ii=1:length(num_WT_vector)
        
        nu_ij_matrix(ii,:) = m_val*(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector)./(K_val+(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector));
        
        phiW_ij_matrix(ii,1) = 1;
        phiW_ij_matrix(ii,2:end) = y(1);
        
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(ii,:).*nu_ij_matrix(ii,:).*prob_psiW_psiD_matrix(ii,:));
        
    end
    
else
    
    for ii=1:length(num_WT_vector)
        
        nu_ij_matrix(ii,:) = m_val*(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector)./(K_val+(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector));
        if ii ==1
            phiW_ij_matrix(ii,2:end) = (num_WT_vector(ii)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end))./(num_WT_vector(ii)*ones(size(num_DI_vector(1,2:end)))+eps*num_DI_vector(1,2:end)+y(1)*(1-eps)*num_DI_vector(1,2:end));
            %         phiW_ij_matrix(i,2:end) = (num_WT_vector(i)*ones(size(num_DI_vector(1,2:end))))./(num_WT_vector(i)*ones(size(num_DI_vector(1,2:end)))+y(1)*num_DI_vector(1,2:end));
            phiW_ij_matrix(1,1) = 1;
        else
            phiW_ij_matrix(ii,:) = (num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector)./(num_WT_vector(ii)*ones(size(num_DI_vector))+eps*num_DI_vector+y(1)*(1-eps)*num_DI_vector);
        end
        sum_W_ij = sum_W_ij + sum(phiW_ij_matrix(ii,:).*nu_ij_matrix(ii,:).*prob_psiW_psiD_matrix(ii,:));
        
    end
    
    
end

% y(2) = f, a scalar multiple to convert from WT to TCID50 measurement
TCID50_output=y(2)*Ncells*sum_W_ij;



