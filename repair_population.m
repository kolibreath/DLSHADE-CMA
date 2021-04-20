function [pop_array,nfes] = repair_population(idx,cluster_size,pop,problem_size,mu,lambda,nfes,func)
% input:
    % idx                       -- indices of clusters
    % pop                       -- global population 
% output:
    % pop_array                 -- array of repaired pop_structs 

    pop_array = [];
    for i = 1 : cluster_size
        pop_cluster = pop(find(idx == i),:);
        [popsize,columns] = size(pop_cluster);
        if popsize < mu     % add some individual to pop_cluster
            % 估计当前popsize 的个体的协方差矩阵，通过协方差矩阵和均值中心再生n个个体
            n = mu - popsize;
            C = cov(pop_cluster(:,1:problem_size));
            sigma = 0.3;
            [B,D] = eig(C);
            B = B';
            xmean = mean(pop_cluster(:,1:problem_size));
            % 生成新的n个个体
            ui = zeros(n,problem_size);
            for k = 1 <= n
                ui(k, :) = (xmean' + sigma * B * (D * randn(problem_size, 1)))';
            end
            ui = evalpop(ui,func);
            nfes = nfes + n;
            pop_cluster = [pop_cluster;ui];
        elseif popsize > mu % remove some individuals from pop_cluster
            % sort based on fitness
            [~,sorted_index] =  sortrows(pop_cluster,columns-1);
            pop_cluster(sorted_index(mu+1:end),:) = [];  
        end

        % 完成修复这个聚类之后 计算协方差矩阵
        C = cov(pop_cluster(:,1:problem_size));
        C = triu(C) + triu(C, 1)';
        [B, D] = eig(C);
        B = B'; 
        % 
        D = sqrt(diag(D));
        invsqrtC = B * diag(D.^ - 1) * B';
        sigma = 0.3;
        eigeneval = nfes; %% TODO 这里要检查eigeneval 什么时候会触发
        xmean = mean(pop_cluster(:,1:problem_size));
        pop_struct = assem_pop(pop_cluster, mu, lambda, problem_size, C, D, B, ...
                     invsqrtC, eigeneval, xmean, sigma);
        pop_array = [pop_array;pop_struct];
    end
end