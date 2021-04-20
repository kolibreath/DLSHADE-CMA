function [pop_array,nfes] = repair_population(idx,cluster_size,pop,problem_size,nfes,func)
% input:
    % idx                       -- indices of clusters
    % pop                       -- global population 
% output:
    % pop_array                 -- array of repaired pop_structs 

    lambda = 4 + floor(3 * log(problem_size));
    popsize = floor(lambda / 2);
    B = eye(problem_size, problem_size);
    D = ones(problem_size, 1);
    C = B * diag(D.^2) * B';
    invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
    eigeneval = 0; % track update of B and D
    sigma = 0.3;

    pop_array = cell(1,cluster_size);
    for i = 1 : cluster_size
        %% 直接对每一个聚类中的所有个体进行协方差矩阵的估计和计算会出现问题：
        % 因为这些个体可以求出协方差矩阵，但是这些个体的特征向量会出现问题
        %% TODO 想出更好的思路进行局部搜索
        %% 目前通过对当前的xmean 作为均值中心进行搜索
        pop_cluster = pop(find(idx == i),:);
        xmean = mean(pop_cluster(:,1:problem_size));
        
        %% 重新初始化pop_cluster
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