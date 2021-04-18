function [pop_struct, archive_fr, archive, suc_f, suc_cr,delta_k] = update_pop_fr(pop_struct,ui,base_index,archive,f,cr)
%COMPARE_FR update population using FROFI
% input:
    % pop_struct   -- struct of population matrix, problem_size etc.
    % ui           -- offspring population of pop
    % base_index   -- indices of base vector 
    % archive      -- archive stores defeated parents
    % f            -- scale factor used by ui
    % cr           -- crossover rate used by ui
% output:
    % pop_struct   -- struct saves updated population 
    % archive_fr   -- defeated offspring in ui_fr
    % archive      -- updated archive
    % suc_f        -- successful scale factor by ui
    % suc_cr       -- successful crossover rate by ui
    % both suc_f and suc_cr are column vector (suc_num * 1)
    % delta_k      -- improvement of successful offspring compared to
    % parent considering both fitness and conv
    
% Steps: 
%  1) if x1 and x2 are both feasible and the one with better fitness will
%  be selected. 
%  2) if x1 is feasible and x2 is not, x1 will be selected.
%  3) if x1 and x2 are both infeasible, the one with less conV will be
%  seleted
%  4) According to FROFI, ones from ui not selected by feasibilty rule but
%  better in fitness will be store in an archive which is used to replace
%  some individuals of pop.

%  NOTE: If the selected one is from pop, it will go to archive (shared by 
%  subpopulations), otherwise the one from ui will be stored at archive_fr

% Version 1.4 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/25
    
%% 
    suc_cr = [];
    suc_f = [];
    archive_fr = [];
    
    delta_k = [];
    lambda = pop_struct.lambda;
    pop = pop_struct.pop;
    
    % compare between offspring individual and parent
    for k = 1 : lambda
       cur_par = pop(base_index(k),:);
       cur_off = ui(k,:);
       
       par_fit = cur_par(end - 1);
       off_fit = cur_off(end - 1);
       
       par_conv = cur_par(end);
       off_conv = cur_off(end);
       %% between two infeasible solutions, the one with smaller conV wins
       if par_conv ~=0 && off_conv ~= 0
          if par_conv > off_conv
             %parent is defeated
             [pop,archive,suc_f,suc_cr, delta_k] = replace_record(pop,base_index(k),cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k);
          else
              % ui is defeated 
             archive_fr = [archive_fr; cur_off];
          end
       %% if both feasible, better fitness is preferred
       elseif par_conv == 0 && off_conv == 0
           if par_fit > off_fit 
              %parent is defeated
              [pop,archive,suc_f,suc_cr, delta_k] = replace_record(pop,base_index(k),cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k);
           else
               % ui is defeated 
               archive_fr = [archive_fr; cur_off];
           end
       %% feasible one is better than its infeasible counterpart
       else 
           if off_conv == 0 
             %parent is defeated
             [pop,archive,suc_f,suc_cr, delta_k] = replace_record(pop,base_index(k),cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k);
           else 
                % ui is defeated 
              archive_fr = [archive_fr; cur_off];
           end
       end
    end
   
    pop_struct.pop = pop;  % rewrite updated population into pop_struct
end
    
