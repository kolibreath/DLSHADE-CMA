function delta_k = delta_k(par_conv, par_fit, off_conv, off_fit)
%DELTA_K calculating delta_k with conV and fitness value of an individual
%compared to original proposed LSHADE, output of delta_k is no longer
%focusing solely on improvement of fitness value. But when using constraint
%handling techique, a "better" individual may not beat its parent in both
%fitness and conV.
% TODO handle documentation
% input:    
    % par_conv          -- parent 
    % off          -- offspring
% output:
    % delta_k      -- combining improvement from both par and off
% Version 1.2 Author: Shi Zeyuan 734780178@qq.com

%% 
    delta_fitness = max(par_fit - off_fit, 0);
    delta_conv    = max(par_conv - off_conv, 0);
    delta_k = [delta_fitness, delta_conv];
end

function dif_val = weights_lshade(delta_k)
% input:
    % delta_k           -- vector (succ_num * 2) of delta_fitness and
    % delta_par
    
% ouput:
    % dif_val
    
    %TODO check here 原始的lshade中使用的是部分占总体的比值，我这里使用了标准化之后再计算总体占全体的比值
    % delta_k(:,1) = delta_fitness, delta_k(:,2) = delta_conv
    max_delta_fit = max(delta_k(:,1));
    min_delta_fit = min(delta_k(:,2));
    
    max_delta_conv = max(delta_k(:,1));
    min_delta_conv = min(delta_k(:,2));
    
    % normalization
    delta_k(:,1) = (max_delta_fit - delta_k(:,1)) ./ (max_delta_fit - min_delta_fit);
    delta_k(:,2) = (max_delta_conv - delta_k(:,2)) ./ (max_delta_conv - min_delta_conv);
    
    dif_val = sqrt(delta_k(:,1) .^ 2 + delta_k(:,2) .^ 2);
    dif_val = dif_val ./ sum(dif_val);
    
end