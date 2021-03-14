function k = mem_update(SF, SCR, M_F, M_CR, k)
    %MEM_UPDATE memory update of scale factor F and crossover rate CR
    
    % input:
        % SF        -- memory save of scale factor F
        % SCR       -- memory save of crossover rate CR
        % M_F       -- mean of Gaussian distribution
        % M_CR      -- mean of Cauchy distribution
        
    % output : 
        % k         -- index of memory save
    ter_value = -1;             % terminal value 
    
    SF_size = length(SF);
    SCR_size = length(SCR);
    
    if SF_size ~= 0 && SCR_size ~= 0
        if MCR(k) == ter_value || max(SCR) == 0
            MCR(k) = ter_value;
        else
            MCR(k) = 
        end 
    
    else
    end
end


