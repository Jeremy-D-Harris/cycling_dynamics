function [W, D, TCID50, HAU] = simulate_passagestudy_models_deterministic_consolidate(params,which_model)

% consider consolidating further but good for now - 05/19/22

max_N_W = 5000; max_N_D = 5000;
N_W=max_N_W; N_D = max_N_D;
num_WT_vector = 1:N_W;
num_DI_vector = 1:N_D;

c_vals = params.c_vals;
tp = params.num_pass;
mu_val = params.mu;
f_val = params.f;

% initial conditions
W_in = params.W_init;
D_in = params.D_init;


% get total viral yield
if ismember(which_model,[1,2])
    
    %input-independent
    nu_val = params.nu;
    nu_i_vector = nu_val*ones(size(num_WT_vector));
    
elseif ismember(which_model,[3,4])
    
    % linear
    m_val = params.m;
    nu_i_vector = m_val*num_WT_vector;
    
elseif ismember(which_model,[5,6])
    
    m_val = params.m;
    K_val = params.K;
    n_val = params.n;
    nu_i_vector = m_val*(num_WT_vector.^n_val)./(K_val^n_val+(num_WT_vector.^n_val));
    
else
    
    m_val = params.m;
    K_val = params.K;
    nu_i_vector = m_val*(num_WT_vector)./(K_val+num_WT_vector);
    
end


% get proportion output function
if ismember(which_model,[1,3,5,7])
    
    % constant proportion WT when WT+DI coinfect
    phi_val = params.phi;
    
    phiW=phi_val;
    phiD=1-phi_val;
    
    phiW_matrix = phiW*ones(N_W,N_D);
    phiD_matrix = phiD*ones(N_W,N_D);
    
else
    
    % relative fitness of DI vs. WT
    lambda_val = params.lambda;
    
    phiW_matrix = zeros(N_W,N_D);
    phiD_matrix = zeros(N_W,N_D);
    for i = 1:N_W
        
        phiW_matrix(i,:) = (num_WT_vector(i)*ones(1,N_D))./(num_WT_vector(i)*ones(1,N_D)+lambda_val*num_DI_vector);
        phiD_matrix(i,:) = lambda_val*num_DI_vector./(num_WT_vector(i)*ones(1,N_D)+lambda_val*num_DI_vector);
    end
    
end


for k=1:tp
    
    % mean MOI's
    psiW = W_in/c_vals(k); psiD = D_in/c_vals(k);
    
    % set up WT vector
    prob_psiW_vector = poisspdf(num_WT_vector,psiW);
    prob_above_N_W = 1 - poisscdf(N_W,psiW);
    prob_psiW_vector(1,N_W) = prob_above_N_W;
    
    % set up DI vector
    prob_psiD_vector = poisspdf(num_DI_vector,psiD);
    prob_above_N_D = 1 - poisscdf(N_D,psiD);
    prob_psiD_vector(1,N_D) = prob_above_N_D;
    
    % joint probability matrix, Poisson(i WT,j DI) = Poisson(i WT)*Poisson(j DI)
    prob_psiW_psiD_matrix = transpose(prob_psiW_vector)*prob_psiD_vector;
    
    % probability of no DI and at least one WT
    alpha = sum(nu_i_vector.*prob_psiW_vector);
    alpha = exp(-psiD)*alpha;
    
    prob_psiW_psiD_phiW_matrix  = prob_psiW_psiD_matrix.*phiW_matrix;
    sum_W_ij = sum(nu_i_vector*prob_psiW_psiD_phiW_matrix);
    
    prob_psiW_psiD_phiD_matrix  = prob_psiW_psiD_matrix.*phiD_matrix;
    sum_D_ij = sum(nu_i_vector*prob_psiW_psiD_phiD_matrix);
    
    W(k,1)=c_vals(k).*(alpha+sum_W_ij)*(1-mu_val)/2;
    D(k,1)=c_vals(k).*((alpha+sum_W_ij)*mu_val+sum_D_ij)/2;
    
    W_in = W(k,1); D_in = D(k,1);
    
end

% make W and D variables row vectors
W = W';
D = D';

% simulate TCID50 and HAU from W and D dynamics
TCID50 = GetTCID50(W, params,f_val);
HAU = GetHAU(W, D, params);
