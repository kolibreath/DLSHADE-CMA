function [pop_struct,nfes,restart_record] = restart_trigger(pop_struct,nfes,func,restart_record)
%RESTART_TRIGGER 此处显示有关此函数的摘要
%   此处显示详细说明
   global sigma_restart;
   global restart_pos;
%% 只要是重新初始化 必然初始化为原来的两倍
%% condition1：如果正定性被破坏 重新初始化
   [popsize,columns] = size(pop_struct.pop);
   problem_size = pop_struct.problem_size;
   % 如果出现非正定的情况 重新初始化种群中的所有 或者是进行某一维度的初始化
   if ~isreal(pop_struct.C) || ~isreal(pop_struct.D) || ~isreal(pop_struct.B)
       restart_record = [restart_record;nfes];
       if rand < restart_pos
           [~,best_fit] = sortrows(pop_struct.pop,columns-1);
            xmean = pop_struct.pop(best_fit,:);% 使用最好的结果作为xmean
            xmean = xmean(1:problem_size);
            lambda = pop_struct.lambda  * 2; 
%             disp("restart");
           [pop_struct,nfes] = initialize_cma_pop(xmean,sigma_restart,lambda,problem_size,nfes,func);       
       else
          [pop_struct,nfes] = dim_restart(pop_struct,func,nfes);
       end
       return;
   end
end

function [pop_struct,nfes] = dim_restart(pop_struct, func, nfes)
%DIM_RESTART 此处显示有关此函数的摘要
%   此处显示详细说明
    pop = pop_struct.pop; % 求列的标准差
    popsize = pop_struct.popsize;
    problem_size = pop_struct.problem_size;
    columns = size(pop, 2);
    fitness = pop(:, columns-1);
    pop = pop(:,1:problem_size);
    standard_deviation = std(pop);
    global deviation_tolerance;
    global dim_restart_sigma;
    index =find(standard_deviation < deviation_tolerance);
    if numel(index) == 0 
        return 
    end
    
    % 如果存在这样的个体，随机选择一个个体的对应的index的维度进行重新初始化
    % 在原来的pop的基础上进行初始化
    for i = 1 : numel(index)
        dim_index = index(i);
        % 找出当前的最好的个体
        % 将最好的个体的这一维度的结果作为均值，正态分布初始化这一个维度中的任意一个值    
        rand_index = randi(popsize);
        [~, best_index] = sort(fitness,'ascend');
        mean_value = pop(best_index(1),dim_index);
        pop(rand_index,dim_index) = normrnd(mean_value, dim_restart_sigma);
    end
    % TODO 这里浪费了nfes 可以修改！！！
    [pop,nfes] = evalpop(pop,func,nfes);
    pop_struct.pop = pop;
end



