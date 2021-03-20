function sorted_pop = sort_fr(pop)
%SORT_FR sort population according to feasibility rule

% Steps:
    % 1) feasible individuals are better than all infeasible ones
    % 2) in all infeasible individuals, less conv is preferred
    % less function value is considered better individuals
    
    feasible_index = find(pop(:, end) == 0); % indicies
    feasible_individuals = pop(feasible_index);
    % sort based on fitness (function value)
    [~, sorted_feasible_indices] = sort(feasible_individuals(end-1), 'ascend');    
    
    infeasible_individuals = pop( ~feasible_index); 
    [~, sorted_infeasible_indices] = sort(infeasible_individuals(end), 'ascend');
    
    sorted_pop = [pop(sorted_feasible_indices);pop(sorted_infeasible_indices)];
    
end

