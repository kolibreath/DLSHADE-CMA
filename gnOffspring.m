function ui = gnOffspring(pop, mu, popsize, f, cr)
%GNOFFSPRING generate offspring by DE/current-to-pbetter*/
% input:
    % pop                   -- population of individuals (pop_fr or pop_ec)
    % mu:                   -- p_best_rate * popsize
    % popsize               -- size of population
    % f                     -- scale factor (vector of popsize * 1)
    % cr                    -- crossover rate (vector of popsize * 1)
% output:
    % ui                    -- individual after mutation and crossover
    
    %% DE/current-to-pbest mutation and crossover
    % pop will be sorted according to feasibilty proportion 
    violated_num = length(find(pop(end) >= 0));
    pfs = violated_num / popsize;     % proportion of violated numbers
    
    [~, columns] = size(pop);  % columns = N + 2 (fitness and conV)
    % TODO 在搜索前期，pfs比较小，需要多选择conv进行排序， 在搜索后期，根据fitness排序
    % TODO 可以采取不同的mutation 策略
    if rand  > pfs % sorted by conV
        pop = sortrows(pop, columns);
    else
        pop = sortrows(pop, columns - 1);
    end
    
    %TODO 添加ranking 根据迭代的情况采用 linear rank selection 或者是exponential
    %top p%
    pNP = max(round(mu), 2); %% choose at least two best solutions
    randindex = ceil(rand(1, popsize) .* pNP); %% select from [1, 2, 3, ..., pNP]
    randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
    pbest = pop(randindex, :); %% randomly choose one of the top 100p% solutions TODO 这个应该是排序之后的选择，检查！

    % mutation
    r0 = [1:popsize];
    popAll = [pop; archive.pop];
    [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);
            
    vi = pop + f(: , ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
    vi = boundConstraint(vi, pop, lu);

    % crossover
    mask = rand(popsize, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
    rows = (1:popsize)'; cols = floor(rand(popsize, 1) * problem_size) + 1; % choose one position where the element of ui doesn't come from the parent
    jrand = sub2ind([popsize problem_size], rows, cols); mask(jrand) = false;
    ui = vi; 
    ui(mask) = pop(mask);
end

