function [ui,r0] = gnOffspring(pop_struct,archive,nfes,max_nfes,f,cr,pattern)
%GNOFFSPRING generate offspring by DE/current-to-pbetter*/
% input:
    % pop_struct                   -- struct of population
    % lu                           -- lower and upper bounds 
    % archive                      -- archvie stores defeated parents
    % f                            -- generated scale factor f
    % cr                           -- generated crossover rate cr
    % pattern                      -- pattern == 1: Global Search Stage; pattern == 2: Local Search Stage
% output:
    % ui                           -- individual after mutation and crossover
    % r0                           -- base vector index

% Version 1.9 Author: Shi Zeyuan 734780178@qq.com Date: 2021/4/19

    pop = pop_struct.pop;
    popsize = pop_struct.popsize;
    problem_size = pop_struct.problem_size;
    pop = pop(:, 1:problem_size);
    popAll = [pop; archive.pop(:,1:problem_size)];
    global pbest_rate;
    global mutation_strategy_de;
    global mutation_strategy_cma;
    global lu ;
    lambda = pop_struct.lambda;
    
    violated_num = length(find(pop(end) >= 0));
    pfs = violated_num / popsize; % proportion of violated numbers

    [~, columns] = size(pop);
    % TODO 检查排序
    if rand >= pfs
        [pop,~] = sortrows(pop, columns-1);
    else
        [pop,~] = sortrows(pop, columns);
    end

    % size of offspring population
    % pattern == 1, lambda = popsize
    % pattern == 2, lambda = popsize * 2
    %% mutation
    if pattern == 1
        %% DE mutation strategy without Covariance Matrix Adaptation
        % DE/rand/1
        if rand < mutation_strategy_de
            r0 = ceil(rand(1,popsize) * popsize);
            [r1, r2] = gnR1R2(popsize,size(popAll,1),r0); 
            base_vectors = pop(r0,:);
            vi = base_vectors + f(:, ones(1, problem_size)) .* (pop(r1, :) - popAll(r2, :));
            vi = boundConstraint(vi, pop, r0, lu);
        else
            % DE/best/1
            temp = ceil(popsize * pbest_rate);
            r0 = ceil(rand(1,popsize) * temp);
            [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);
            base_vectors = pop(r0,:);
            vi = base_vectors + f(:, ones(1,problem_size)) .* (pop(r1, :) - popAll(r2, :));
            vi = boundConstraint(vi, pop, r0, lu);
        end
    else  % pattern == 2
        %% DE mutation strategy with Covariance Matrix Adaptation 
        % pop will be sorted according to feasibilty proportion

        % select lambda parent as base vector
        lambda = pop_struct.lambda;
        r0 = ceil(rand(1, lambda) * popsize);
        [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);
        pbetter = zeros(lambda, problem_size);

        for k = 1:lambda
            % 可能有部分来自best
            if rand < mutation_strategy_cma
                pbetter(k, :) = (pop_struct.xmean' + pop_struct.sigma ...
                    * pop_struct.B * (pop_struct.D .* randn(problem_size, 1)))';
            else
                temp = ceil(pbest_rate * popsize);
                pbetter(k, :) = pop(randi(temp), :);
            end
        end

        base_vectors = pop(r0, :);
        f_w = gnFw(nfes, max_nfes, f);
        vi = base_vectors + f_w(:, ones(1, problem_size)) .* (pbetter - base_vectors) ...
            + f(:, ones(1, problem_size)) .* (pop(r1, :) - popAll(r2, :));
        vi = boundConstraint(vi, pop, r0, lu);
    end
  
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


