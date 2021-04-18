function [restart_index_fr, restart_index_fo] = stop_trigger(last_bsf,bsf,pop_fr_struct,pop_fo_struct)
%STOP_TRIGGER check if cma meets the criteria of restart
% input:
    % pop_struct            -- struct of population information
    % bsf                   -- the current best so far solution
    
% output:
    % valid                 -- valid == 1 -> restart; valid == 0 -> no restart
    
    Tolfun = 1e-12;
    % the "best so far solution" in the last 10 + floor(30 * problem_size / lambda) generation 
    
    %% condition 1: 
    if abs(last_bsf(end-1) - bsf(end-1)) < Tolfun && last_bsf(end-1) ~= inf
       restart_index_fo = 1; restart_index_fr = 1;
       return;
    end
    
    %% condition 2 , condition 3, condition 4
    restart_index_fr = check_subpopulation(pop_fr_struct);
    restart_index_fo = check_subpopulation(pop_fo_struct);
    
end

function logic = check_subpopulation(pop_struct)
 %% condition 2:
    TolX = 1e-12 * 0.3; % sigma_0 = 0.3
    % number of elements in the standard deviation matrix that not smaller than TolX
    n_sd = numel(find(pop_struct.C >= TolX));
    % number of elements in the sigma*pc vector that not smaller than TolX
    n_sc = numel(find(pop_struct.sigma * pop_struct.pc >= TolX));
    n_sd_matrix = numel(pop_struct.C);
    n_sc_vector = numel(pop_struct.pc);
    
    if n_sd > ceil(n_sd_matrix * 0.1) && n_sc > ceil(n_sc_vector * 0.3) 
        logic = 1;
        return;
    end
    
    %% condition 3:
    if any(pop_struct.xmean == pop_struct.xmean + 0.1 * pop_struct.sigma * pop_struct.D * pop_struct.B)
        logic = 1;
    end
    
    %% condition 4:
    if max(pop_struct.D) > 1e14*min(pop_struct.D)
        logic = 1;
    end
end


