function [pop, archive, suc_f, suc_cr, delta_k] = replace_record(pop,k,cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k)                                                               
% REPLACE_RECORD replace defeated parent into the archive and record successful f and cr
% input: 
    % pop       -- population
    % k         -- k_th parent in population is replaced
    % cur_off   -- current offspring
    % cur_par   -- current parent
    % archive   -- archive saving defeating parents
    % suc_f     -- storing successful scale factor f
    % suc_cr    -- storing successful crossover rate cr
    % delta_k   -- vector recording improvement in fitness and conv
% output:
    % pop       -- population
    % archive   -- archive saving defeating parents
    % suc_f     -- storing successful scale factor f
    % suc_cr    -- storing successful crossover rate cr
    % delta_k   -- updated delta_k (append a new row of fitness and conv improvement)
    
% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

%%
    % find those from offspring population whose fitness is better than
    % their parent but defeated in Deb's feasibility rule selection 
% scale factor F and crossover rate cr
  pop(k, :) = cur_off;
  archive.pop = [archive.pop; cur_par];
  
  % record successful f and cr
  % TODO 使用最好的子代替换
  suc_f = [suc_f; f(k)];
  suc_cr = [suc_cr;cr(k)];
  
  % delta_k(par_conv, par_fit, off_conv, off_fiWt)
  delta_fitness = max(cur_par(end-1) - cur_off(end-1), 0);
  delta_k =  [delta_k;delta_fitness];
end