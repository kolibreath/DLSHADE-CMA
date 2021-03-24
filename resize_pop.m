function pop_struct = resize_pop(max_popsize,min_popsize,pop_struct,max_nfes,nfes)
%RESIZE_POP resize population using LSPR strategy
% input:
    % max_popsize           -- maximum popsize of pop
    % min_popsize           -- minimum popsize of pop
    % pop_struct            -- population struct
    % max_nfes              -- maximum of function evaluation
    % nfes                  -- current funciton evaluation
    % archive               -- archive of eliminated individuals
% output:
    % pop                   -- updated population after resizing
% Steps
    % 1) output the number of individual should be deleted and denoted as reducation_ind_num
    % 2) sort population based on proportion of feasibility solutions, in the
    % early stages of evolution, population moves in the direction of feasible
    % regions, and in the later stages of evolution, individuals converge in
    % feasible regions and look for the best solution by deleting worst
    % solution based on the rank of fitness
    
% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

    %% for resizing the population size
    pop = pop_struct.pop;
    plan_popsize = round((((min_popsize - max_popsize) / max_nfes) * nfes) + max_popsize);
    popsize = pop_struct.popsize;
    if popsize > plan_popsize
       % prevent popsize after resize is smaller than min_popsize
       reduction_ind_num = min(popsize - plan_popsize, popsize - min_popsize);
       
       k = 0; 
       while k < reduction_ind_num
           worst_ind = popsize - k;
           pop(worst_ind) = [];
           k = k+1;
       end 

       pop_struct.popsize = popsize - reduction_ind_num;

    end
    pop_struct.pop = pop;
end

