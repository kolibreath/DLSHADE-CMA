function [pop_struct,archive_ec,archive,suc_f,suc_cr, delta_k] = update_pop_ec(pop_struct,ui,archive,f,cr,epsilon)
%UPDATE_POP_EC update population using epsilon constraint handing (original)
% input:
    % pop_struct   -- population struct
    % ui           -- offspring population of pop
    % archive      -- archive stores defeated parents
    % f            -- scale factor used by ui
    % cr           -- crossover rate used by ui
% output:
    % pop_struct   -- updated population has been already stored here
    % archive_ec   -- defeated offspring in ui_ec
    % archive      -- updated archive
    % suc_f        -- successful scale factor by ui
    % suc_cr       -- successful crossover rate by ui
    % both suc_f and suc_cr are column vector (suc_num * 1)
    % delta_fit      -- improvement of successful offspring compared to
    % parent considering both fitness and conv

% Steps
    % x1 is better than x2 when 
    % 1) f(x1) < f(x2) if conV(x1), conV(x2) < epsilon
    % 2) f(x1) < f(x2) if conV(x1) == conV(x2)
    % 3) conV(x1) < conV(x2) otherwise

% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18
%%
    delta_k = [];
    archive_ec = [];
    popsize = pop_struct.popsize;
    pop = pop_struct.pop;
    
    suc_f = [];
    suc_cr = [];
    
    for k = 1 : popsize
        cur_off = ui(k, :);
        cur_par = pop(k, :);
        
        fit_off = cur_off(end-1);
        fit_par = cur_par(end-1);
        
        conv_off = cur_off(end);
        conv_par = cur_par(end);
        % f(x1) < f(x2) if conV(x1), conV(x2) < epsilon
        if conv_off < epsilon && conv_par < epsilon
            if fit_off < fit_par
                % defeated parent
                 [pop,archive,suc_f,suc_cr,delta_k] = replace_record(pop,k,cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k);
            else
                % put it into archive_ec
                archive_ec = [archive_ec; cur_off];
            end
        % f(x1) < f(x2) if conV(x1) == conV(x2)
        elseif conv_off == conv_par
            if fit_off < fit_par
                % defeated parent
                [pop,archive,suc_f,suc_cr,delta_k] = replace_record(pop,k,cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k);
            else
                % put it into archive_ec
                archive_ec = [archive_ec; cur_off];
            end
        else
            if conv_off < conv_par
                 % defeated parent
                 [pop,archive,suc_f,suc_cr,delta_k] = replace_record(pop,k,cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k);
            else
                % put it into archive_ec
                archive_ec = [archive_ec; cur_off];
            end
        end
    end
   
    pop_struct.pop = pop;
end


