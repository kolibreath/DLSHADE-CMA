function [pop_array,nfes] = repair_population(pop,problem_size,nfes,func)
% input:
    % idx                       -- indices of clusters
    % pop                       -- global population 
% output:
    % pop_array                 -- array of repaired pop_structs 

    [cluster_size,centers] = find_cluster_center(pop);
    
    lambda = 4 + floor(3 * log(problem_size));
    popsize = floor(lambda / 2);
    B = eye(problem_size, problem_size);
    D = ones(problem_size, 1);
    C = B * diag(D.^2) * B';
    invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
    eigeneval = nfes; % track update of B and D
    sigma = 0.3;

    pop_array = cell(1,cluster_size);
    for i = 1 : cluster_size
        %% 直接对每一个聚类中的所有个体进行协方差矩阵的估计和计算会出现问题：
        % 因为这些个体可以求出协方差矩阵，但是这些个体的特征向量会出现问题
        %% TODO 想出更好的思路进行局部搜索
        %% 目前通过对当前的xmean 作为均值中心进行搜索
        xmean = centers(i,1:problem_size);
        
        pop_cluster = zeros(popsize,problem_size);
        k = 1;
        while k <= popsize
            pop_cluster(k,:) =  (xmean' + sigma * B * (D .* randn(problem_size, 1)))';
            k = k + 1;
        end
        pop_cluster = evalpop(pop_cluster,func);
        nfes = nfes + popsize;
        pop_struct = assem_pop(pop_cluster, popsize, lambda, problem_size, C, D, B, ...
                     invsqrtC, eigeneval, xmean, sigma);
        pop_array{i} = pop_struct;
    end
end

% 基本上是分成四个种群
% best_fitness best_conv mean_conv mean_fit 如果当best_fitness 和 best_conv 对应的个体是一个就会分成三个钟群
function [cluster_size, centers] = find_cluster_center(pop)
    [popsize, columns] = size(pop);
    [~,fit_index] = sortrows(pop,columns-1);
    [~,conv_index] = sortrows(pop,columns);
    if fit_index(1) == conv_index(1) 
        cluster_size = 3;
        centers = zeros(3,columns);
        centers(1) = pop(:, pop(fit_index(1)));
        % 找mean_fit mean_fit 如果使用两个向量的中间值会重新计算，没有必要
        centers(2) = pop(:, pop(fit_index(ceil(popsize / 2))));
        centers(3) = pop(:, pop(conv_index(ceil(popsize / 2))));
    else 
        cluster_size = 4;
        centers = zeros(4, columns);
        centers(1) = pop(:, pop(fit_index(1)));
        centers(2) = pop(:, pop(conv_index(1)));
        centers(3) = pop(:, pop(fit_index(ceil(popsize / 2))));
        centers(4) = pop(:, pop(conv_index(ceil(popsize / 2))));
    end

end