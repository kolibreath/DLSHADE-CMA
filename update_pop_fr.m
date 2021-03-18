function [pop, archive_fr, archive, suc_f, suc_cr] = update_pop_fr(pop,ui,archive,f,cr)
%COMPARE_FR implementing Deb's feasibility rule 
% input:
    % pop          -- population
    % ui           -- offspring population of pop
    % archive      -- archive stores defeated parents
    % f            -- scale factor used by ui
    % cr           -- crossover rate used by ui
% output:
    % pop          -- population 
    % archive_fr   -- defeated offspring in ui_fr
    % archive      -- updated archive
    % suc_f        -- successful scale factor by ui
    % suc_cr       -- successful crossover rate by ui
    % both suc_f and suc_cr are column vector (suc_num * 1)
    
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

% Version 1.2 Author: Shi Zeyuan 734780178@qq.com
    
    % find those from offspring population whose fitness is better than
    % their parent but defeated in Deb's feasibility rule selection 
    
    record_index = find(ui(end) > pop(end) & ui(end-1) < pop(end-1));
    
    suc_cr = [];
    suc_f = [];
    archive_fr = [];
    
    for k = 1 : pop.popsize
       cur_par = pop(k);
       cur_off = ui(k);
       
       par_fit = cur_par(end - 1);
       off_fit = cur_off(end - 1);
       
       par_conv = cur_par(end);
       off_conv = cur_off(end);
       
       % TODO extract common logic
       % between two infeasible solutions, the one with smaller smaller conV wins
       if par_conv ~=0 && off_conv ~= 0
          if par_conv > off_conv
              % replace and store in archive
              pop(k, :) = cur_off;
              archive = [archive; cur_par];
              % record successful f and cr
              suc_f = [suc_f; f(k)];
              suc_cr = [suc_cr;cr(k)];
          else
              % ui is defeated 
              archive_fr = [archive_fr; ui];
          end
       % if both feasible, better fitness is preferred
       elseif par_conv == 0 && off_conv == 0
           
           if par_fit > off_fit 
               % replace and store in archive
               pop(k, :) = cur_off;
               archive = [archive; cur_par];
               % record successful f and cr
               suc_f = [suc_f; f(k)];
               suc_cr = [suc_cr;cr(k)];
           else
               % ui is defeated 
               archive_fr = [archive_fr; ui];
           end
       % feasible one is better than its infeasible counterpart
       else 
           if off_conv == 0 
              % replace and store in archive
              pop(k, :) = cur_off;
              archive = [archive; cur_par];
              % record successful f and cr
              suc_f = [suc_f; f(k)];
              suc_cr = [suc_cr;cr(k)];
           else 
               % ui is defeated 
               archive_fr = [archive_fr; ui];
           end
       end
    end
   
    % applying FROFI replacement strategy
    

