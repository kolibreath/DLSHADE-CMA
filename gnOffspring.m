function [ui,r0] = gnOffspring(pop_struct,lu,archive,nfes,max_nfes,f,cr)
%GNOFFSPRING generate offspring by DE/current-to-pbetter*/
% input:
    % pop_struct                   -- struct of population
    % lu                           -- lower and upper bounds 
    % archive                      -- archvie stores defeated parents
    % f                            -- generated scale factor f
    % cr                           -- generated crossover rate cr
% output:
    % ui                           -- individual after mutation and crossover
    % r0                           -- base vector index

% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

    %% DE/current-to-pbest* mutation and crossover
    pop = pop_struct.pop;
    popsize = pop_struct.popsize;
    problem_size = pop_struct.problem_size;
    
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
    
    %% mutation
    % select lambda parent as base vector 
    lambda = pop_struct.lambda;
    r0 = ceil(rand(1,lambda) * popsize);
    popAll = [pop; archive.pop];
    [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);

    pop = pop(:,1:problem_size);
    pbetter = zeros(lambda,problem_size);
    %TODO 设置这里pbest 概率 为 0.3
    pbest_rate = 0.3;
    for k = 1: lambda
        % 可能有部分来自best
        if rand < 0.5
            pbetter(k,:) = (pop_struct.xmean' + pop_struct.sigma ...
              * pop_struct.B * (pop_struct.D .* randn(problem_size, 1)))';
        else
            temp = ceil(pbest_rate * popsize);
            pbetter(k,:) = pop(ceil(temp * rand),:);
        end
    end
    
    popAll = popAll(:,1:problem_size);
    
    base_vectors = pop(r0,:);
    f_w = gnFw(nfes,max_nfes,f);
    vi = base_vectors + f_w(: , ones(1, problem_size)) .* (pbetter - base_vectors)...
       + f(: , ones(1, problem_size)).* (pop(r1, :) - popAll(r2, :));
    vi = boundConstraint(vi,pop,r0,lu);
    % crossover
    mask = rand(lambda, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
    rows = (1:lambda)'; 
    cols = floor(rand(lambda, 1) * problem_size) + 1; % choose one position where the element of ui doesn't come from the parent
    jrand = sub2ind([lambda problem_size], rows, cols); mask(jrand) = false;
    ui = vi; 
    ui(mask) = base_vectors(mask);
    
end

function f_w = gnFw(nfes,max_nfes,f)
    if nfes <= floor(0.2 * max_nfes)
        f_w = 0.7*f;
    elseif nfes > floor(0.2 * max_nfes) && nfes <= floor(0.4 * max_nfes)
        f_w = 0.8*f;
    else
        % TODO 改了这里 1.2
        f_w = 0.5*f;
    end
end


