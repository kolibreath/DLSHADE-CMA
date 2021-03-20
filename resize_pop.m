function [pop] = resize_pop(max_popsize,min_popsize,pop,max_nfes,nfes)
%RESIZE_POP resize population using LSPR strategy
% input:
    % max_popsize           -- maximum popsize of pop
    % min_popsize           -- minimum popsize of pop
    % pop                   -- population (pop_fr, pop_ec)
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
%   此处显示详细说明
  %% for resizing the population size
    plan_popsize = round((((min_popsize - max_popsize) / max_nfes) * nfes) + max_popsize);
    popsize = pop.popsize;
    if popsize > plan_popsize
       % prevent popsize after resize is smaller than min_popsize
       reduction_ind_num = min(popsize - plan_popsize, popsize - min_popsize);
       popsize = popsize - reduction_ind_num;
       pop.popsize = popsize;

       % find PFS
       % TODO 使用这样的方式是否可行?
       feasible_nums = length(find(pop(end) == 0));
       pfs = feasible_nums / popsize;
       
       if rand > pfs % early stage of evolution
           conv = pop(:, end);
           [~, indBest] = sort(conv, 'ascend');
       else          % later stage of evolution
           fitness = pop(:, end);
           [~, indBest] = sort(fitness, 'ascend');
       end
       
       for r = 1:reduction_ind_num
            worst_ind = indBest(end);
            pop(worst_ind, :) = [];
       end

    end
end

