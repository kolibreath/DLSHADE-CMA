function sorted_pop = sort_fr(pop)
%SORT_FR sort population according to feasibility rule

% Steps:
    % 1) feasible individuals are better than all infeasible ones
    % 2) in all infeasible individuals, less conv is preferred
    % less function value is considered better individuals

% Version 1.4 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

    feasible_index = find(pop(:, end) == 0); % indicies
    % feasible_index will be a zero vector, if there is no feasible individuals
    if ~any(feasible_index) 
        sorted_feasible_indices = [];
    else
        feasible_individuals = pop(feasible_index,:);
        % sort based on fitness (function value)
        [~, sorted_feasible_indices] = sort(feasible_individuals(:,end-1), 'ascend');    
    end
    [popsize,~] = size(pop);
    all_index = 1:popsize;
    all_index(feasible_index) = [];
    infeasible_index = all_index;
    infeasible_individuals = pop(infeasible_index,:); 
    if ~any(infeasible_individuals)
        sorted_infeasible_indices = [];
    else     
        [~, sorted_infeasible_indices] = sort(infeasible_individuals(:,end-1), 'ascend');
    end
    
    sorted_pop = [pop(sorted_feasible_indices,:);pop(sorted_infeasible_indices,:)];
    
end

