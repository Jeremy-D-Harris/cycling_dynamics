function GE_output = Get_GE_MM_wDIP(WT_moi, DIP_moi, m_fit, K_fit,eps, Ncells)

max_N_W = 1000; N_W=max_N_W;
max_N_D = 1000; N_D=max_N_D;

num_WT_vector = 0:N_W;
num_DIP_vector = 0:N_D;

GE_output = zeros(size(WT_moi));

% for i =1:length(WT_moi)

psiW = WT_moi;

% set up WT vector
prob_psiW_vector = poisspdf(num_WT_vector,psiW);
prob_above_N_W = 1 - poisscdf(N_W,psiW);
prob_psiW_vector(1,N_W) = prob_above_N_W;


for k = 1:length(DIP_moi)
    
    psiD = DIP_moi(1,k);
    
    % set up DIP vector
    prob_psiD_vector = poisspdf(num_DIP_vector,psiD);
    prob_above_N_D = 1 - poisscdf(N_D,psiD);
    prob_psiD_vector(1,N_D) = prob_above_N_D;
    
    for j = 1:(N_D+1)
        
        % total cellular output as a function of WT only
        nu_i_vector_j = (num_WT_vector+eps*(j-1)*ones(size(num_WT_vector)))./(K_fit+num_WT_vector+eps*(j-1)*ones(size(num_WT_vector)));
        
        % probability of no DI and at least one WT
        alpha(j) = sum(nu_i_vector_j.*prob_psiW_vector)*prob_psiD_vector(1,j);
        
    end
    
    % GE output as a function of WT and DIP MOI
    GE_output(1,k)=m_fit*Ncells*sum(alpha);
    
end

% end


