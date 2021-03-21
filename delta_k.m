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
% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

%% 
    delta_fitness = max(par_fit - off_fit, 0);
    delta_conv    = max(par_conv - off_conv, 0);
    delta_k = [delta_fitness, delta_conv];
end
