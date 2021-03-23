function sorted_ec = sort_ec(pop, epsilon)
%SORT_EC sort pop based on epsilon constriant
% input: 
    % pop           -- population
    % epsilon       -- epsilon constraint value
% output:
    % sorted_ec     -- sorted population using epsilon_constriant method
% Steps
    % x1 is better than x2 when
    % 1) f(x1) < f(x2) if conV(x1), conV(x2) < epsilon
    % 2) f(x1) < f(x2) if conV(x1) == conV(x2)
    % 3) conV(x1) < conV(x2) otherwise

% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18
    
    %TODO further check here!
    conv = pop(:, end);
    
    [popsize,columns] = size(pop);

    % surely conv of feasible individuals is less than infeasible ones
    feasible_index = find(conv == 0);
    if ~any(feasible_index)
        feasible_individuals = [];
    else
        feasible_individuals = pop(feasible_index,:);
        % based on fitness
        [~, sorted_index] = sort(feasible_individuals(:,end - 1), 'ascend');
        feasible_individuals = feasible_individuals(sorted_index,:);
    end

    % find individuals that conv is less than epsilon value, but exclude (absolutely) feasible ones
    epsilon_feasible_index = find(xor(conv <= epsilon, conv == 0));
    if ~any(epsilon_feasible_index)
        epsilon_feasible_individuals = [];
    else
        epsilon_feasible_individuals = pop(epsilon_feasible_index,:);
        % based on fitness
        [~, sorted_index] = sort(epsilon_feasible_individuals(:,end - 1), 'ascend');
        epsilon_feasible_individuals = epsilon_feasible_individuals(sorted_index,:);
    end
    

    % find conv more than or equal to epsilon but conv equal
    % 1) find all conv > epsilon
    epsilon_infeasible_index = find(conv > epsilon);
    % 2) find conv > epsilon and equal 
    epsilon_infeasible_equal_individuals = pop(epsilon_infeasible_index,:);
    epsilon_infeasible_equal_index = [];
    % if epsilon_infeasible_equal_individual is not empty
    if numel(epsilon_infeasible_equal_individuals) ~= 0
        [~, i, ~] = unique(epsilon_infeasible_equal_individuals(end), 'first');
        index_of_duplicates = find(not(ismember(1:numel(epsilon_infeasible_equal_individuals(end)), i)));
        epsilon_infeasible_equal_individuals = []; 
        % if there exist any non-zero element
        if any(index_of_duplicates) 
            % duplicates of constraint violation value
            % eg A = [1 1 2 3 4 3], duplicate_conv = [1 3]
            duplicate_conv = conv(index_of_duplicates);
            % sort duplicate_conv, individuals with less conv will be placed before others
            duplicate_conv = sort(duplicate_conv, 'ascend');
            % pass feasible and epsilon_feasible index and make following code run faster
            conv_copy = conv;
            delete_index = find(conv < epsilon); % feasible and epsilon-feasible
            conv_copy(delete_index) = inf; 
            len = length(conv_copy);
            
            % Radix sort: 1) ascend conv 2) ascend fitness
            for dupli = duplicate_conv
                temp = [];
                for i = 1:len
                    if conv_copy(i) == inf; continue; end  % pass feasible and epsilon-feasible for higher speed
                    if dupli == conv_copy(i)
                        temp = [temp; pop(i,:)];
                        epsilon_infeasible_equal_index = [epsilon_infeasible_equal_index; i]; 
                    end
                end
                % sort temp based on fitness
                temp = sortrows(temp, columns-1);
                epsilon_infeasible_equal_individuals = [epsilon_infeasible_equal_individuals;temp];
            end
        end
    end
    
    % otherwise
    selected_index = [feasible_index; epsilon_feasible_index; epsilon_infeasible_equal_index];
    % check unique
    uniq = length(selected_index) - length(unique(selected_index));
    if uniq ~= 0
        ME = MException("duplicates in selected index");
        throw(ME);
    end
    
    pop(selected_index,:) = [];
    % based on conv
    pop = sortrows(pop,columns);

    % combine all the index from the best to worse:
    % absolutely feasible, epsilon_feasible, epsilon_infeasible but have
    % equal conv, otherwise
    sorted_ec = [feasible_individuals;epsilon_feasible_individuals;epsilon_infeasible_equal_individuals;pop];
end
