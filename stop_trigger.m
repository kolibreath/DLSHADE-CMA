function [restart_index_fr, restart_index_fo,bsf_unchanged_counter] = stop_trigger(bsf_unchanged_counter,pop_fr_struct,pop_fo_struct)
%STOP_TRIGGER check if cma meets the criteria of restart
% input:
    % pop_struct            -- struct of population information
    % bsf                   -- the current best so far solution
    
% output:
    % valid                 -- valid == 1 -> restart; valid == 0 -> no restart
    
    % the "best so far solution" in the last 10 + floor(30 * problem_size / lambda) generation 
    
    %% condition 1: 
    if bsf_unchanged_counter > 50
       restart_index_fo = 1; 
       restart_index_fr = 1;
       bsf_unchanged_counter = 0;
       return;
    end
    
    %% condition 2 , condition 3, condition 4
    restart_index_fr = check_subpopulation(pop_fr_struct);
    restart_index_fo = check_subpopulation(pop_fo_struct);
    
end

function logic = check_subpopulation(pop_struct)
    logic = 0;
   %% condition 2:
    TolX = 1e-21 * 0.3; % sigma_0 = 0.3
    n_sd = all(all(abs(pop_struct.C) <= TolX));
    n_sc = all(abs(pop_struct.sigma * pop_struct.pc) <= TolX);
   
    if n_sd && n_sc 
        logic = 1;
        return;
    end
    
    %% condition 3:
    if numel(find(pop_struct.xmean == pop_struct.xmean + 0.1 * ...
            pop_struct.sigma * pop_struct.D' * pop_struct.B)) >= 2
        logic = 1;
        return;
    end
    
    %% condition 4:
    if max(pop_struct.D) > 1e14*min(pop_struct.D)
        logic = 1;
        return;
    end
end