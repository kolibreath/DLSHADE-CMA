function [pop_array,nfes,cluster_size] = create_cma_pop(pop,problem_size,nfes,func,lu)
% input:
    % idx                       -- indices of clusters
    % pop                       -- global population 
% output:
    % pop_array                 -- array of repaired pop_structs 

    idx = min_max_distance(pop,problem_size);
    [cluster_size,centers] = find_cluster_center(pop(idx,:));

    pop_array = cell(1,cluster_size);
    for i = 1 : cluster_size
        %% 直接对每一个聚类中的所有个体进行协方差矩阵的估计和计算会出现问题：
        % 因为这些个体可以求出协方差矩阵，但是这些个体的特征向量会出现问题
        %% TODO 想出更好的思路进行局部搜索
        %% 目前通过对当前的xmean 作为均值中心进行搜索
         xmean = centers(i,1:problem_size);
        [pop_struct,nfes] = initialize_cma_pop(xmean,0.3,problem_size,nfes,func,lu);
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
        centers(1,:) = pop(fit_index(1),:);
        % 找mean_fit mean_fit 如果使用两个向量的中间值会重新计算，没有必要
        centers(2,:) = pop(fit_index(ceil(popsize / 2)),:);
        centers(3,:) = pop(conv_index(ceil(popsize / 2)),:);
    else 
        cluster_size = 4;
        centers = zeros(4, columns);
        centers(1,:) = pop(fit_index(1),:);
        centers(2,:) = pop(conv_index(1),:);
        centers(3,:) = pop(fit_index(ceil(popsize / 2)),:);
        centers(4,:) = pop(conv_index(ceil(popsize / 2)),:);
    end

end

%% 通过最大最小距离确定最多4*lambda 大小的种群
function idx = min_max_distance(pop,problem_size)
    % 初始点设定为最优fitness的点
    [popsize, columns] = size(pop);
    [~,fit_index] = sortrows(pop,columns-1);
    best_index = fit_index(1);
    lambda = 4 + floor(3 * log(problem_size));
    max_n = 4 * lambda;
    best_individual = pop(best_index,:);

    idx = zeros(1,max_n);
    % 遍历所有其他个体找到里自身距离最小的
    for k = 1 : max_n
        pop(best_index,:) = []; % 删除当前best
        [popsize,~] = size(pop);
        idx(k) = best_index;
        diff = (pop(:,1:problem_size) - repmat(best_individual(1:problem_size),popsize,1)) .^2;
        distance = sum(diff,2);
        best_index = min(distance);
        best_individual = pop(best_index,:);
    end
end