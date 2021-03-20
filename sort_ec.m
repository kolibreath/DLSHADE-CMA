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

    % surely conv of feasible individuals is less than infeasible ones
    feasible_index = find(conv == 0);
    feasible_individuals = pop(feasible_index);
    % based on fitness
    [feasible_individuals, sorted_feasible_index] = sort(feasible_individuals(end - 1), 'ascend');

    % find individuals that conv is less than epsilon value, but exclude (absolutely) feasible ones
    epsilon_feasible_index = find(xor(conv < epsilon, conv == 0));
    epsilon_feasible_individuals = pop(epsilon_feasible_index);
    % based on fitness
    [epsilon_feasible_individuals, sorted_epsilon_feasible_index] = sort(epsilon_feasible_individuals(end - 1), 'ascend');

    % find conv more than or equal to epsilon but conv equal
    % 1) find all conv >= epsilon
    epsilon_infeasible_index = find(conv >= epsilon);
    % 2) find conv >= epsilon and equal
    epsilon_infeasible_equal_individuals = pop(epsilon_infeasible_index);
    [~, i, ~] = unique(epsilon_infeasible_equal_individuals(end), 'first');
    index_of_duplicates = find(not(ismember(1:numel(epsilon_infeasible_equal_individuals(end)), i)));
    % duplicates of constraint violation value
    % eg A = [1 1 2 3 4 3], duplicate_conv = [1 3]
    duplicate_conv = conv(index_of_duplicates);
    % sort duplicate_conv, individuals with less conv will be placed before others
    duplicate_conv = sort(duplicate_conv, 'ascend');
    % pass feasible and epsilon_feasible index and make following code run faster
    conv_copy = conv;
    delete_index = find(conv < epsilon); % feasible and epsilon-feasible
    conv_copy(delete_index) = -1; 
    len = length(conv_copy);

    epsilon_infeasible_equal_individuals = [];
    epsilon_infeasible_equal_index = [];
    for dupli = duplicate_conv
        for i = 1:len
            if conv_copy(i) == -1; continue; end
            if dupli == conv_copy(i)
                epsilon_infeasible_equal_individuals = [epsilon_infeasible_equal_individuals; pop(i)];
                epsilon_infeasible_equal_index = [epsilon_infeasible_equal_index; i]; 
            end
        end
        % based on fitness
        [epsilon_infeasible_indviduals, ~] = sort(epsilon_infeasible_equal_individuals(end - 1), 'ascend');
    end

    % otherwise
    pop_copy = pop;
    selected_index = [feasible_index; epsilon_feasible_index; epsilon_infeasible_equal_index];
    pop_copy(selected_index) = [];
    otherwise_individuals = pop_copy;
    % based on conv
    [~, sorted_otherwise_index] = sort(otherwise_individuals(end), 'ascend');

    % combine all the index from the best to worse:
    % absolutely feasible, epsilon_feasible, epsilon_infeasible but have
    % equal conv, otherwise
    sorted_ec = [feasible_individuals;epsilon_feasible_individuals;epsilon_infeasible_equal_individuals;otherwise_individuals];
end
