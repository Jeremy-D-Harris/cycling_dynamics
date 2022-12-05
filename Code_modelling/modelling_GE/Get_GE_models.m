function GE_output = Get_GE_models(moi, x, Ncells, which_model)

GE_output = zeros(size(moi));

switch which_model
    
    % input-independent model
    case 1
    
    nu = x(1);
    GE_output = Ncells*nu*(1-exp(-moi));
    
        
    % linear model
    case 2
    
    m = x(1);
    GE_output = Ncells*m*moi;
    
        
    % MM model
    case 3
        
        m_fit = x(1); K_fit = x(2);
        
        max_N_W = 1000; N_W=max_N_W;
        
        num_WT_vector = 1:N_W;
        
        % total cellular output as a function of WT only
        nu_i_vector = num_WT_vector./(K_fit+num_WT_vector);
        
        for k =1:length(moi)
            
            psiW = moi(k,1);
            
            % set up WT vector
            prob_psiW_vector = poisspdf(num_WT_vector,psiW);
            prob_above_N_W = 1 - poisscdf(N_W,psiW);
            prob_psiW_vector(1,N_W) = prob_above_N_W;
            
            % probability of no DI and at least one WT
            alpha = sum(nu_i_vector.*prob_psiW_vector);
            
            % GE output as a function of WT MOI
            GE_output(k,1)=m_fit*Ncells*alpha;
            
        end
        
        
    % Hill model
    case 4
        
        m_fit = x(1); K_fit = x(2); n_fit = x(3);
        
        max_N_W = 1000; N_W=max_N_W;
        
        num_WT_vector = 1:N_W;
        
        % total cellular output as a function of WT only
        nu_i_vector = (num_WT_vector.^n_fit)./(K_fit^n_fit+(num_WT_vector.^n_fit)); % updated the model to reflect the form: m*x^n/(K^n+x^n)
        
        for k =1:length(moi)
            
            psiW = moi(k,1);
            
            % set up WT vector
            prob_psiW_vector = poisspdf(num_WT_vector,psiW);
            prob_above_N_W = 1 - poisscdf(N_W,psiW);
            prob_psiW_vector(1,N_W) = prob_above_N_W;
            
            % probability of no DI and at least one WT
            alpha = sum(nu_i_vector.*prob_psiW_vector);
            
            % GE output as a function of WT MOI
            GE_output(k,1)=m_fit*Ncells*alpha;
            
        end
        
end



