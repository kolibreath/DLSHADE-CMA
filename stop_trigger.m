function logic = stop_trigger(last_bsf,bsf,pop_struct)
%STOP_TRIGGER check if cma meets the criteria of restart
% input:
    % pop_struct            -- struct of population information
    % bsf                   -- the current best so far solution
    
% output:
    % valid                 -- valid == 1 -> restart; valid == 0 -> no restart
    
    Tolfun = 1e-12;
    % the "best so far solution" in the last 10 + floor(30 * problem_size / lambda) generation 
    logic = 0;
    %% condition 1: 
    if abs(last_bsf(end-1) - bsf(end-1)) < Tolfun
       logic = 1;
       return;
    end
    %% condition 2:
    TolX = 1e-12 * 0.3; % sigma_0 = 0.3
end

