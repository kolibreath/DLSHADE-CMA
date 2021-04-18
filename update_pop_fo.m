function [pop_struct, archive_fo, archive, suc_f, suc_cr, delta_k] = update_pop_fo(pop_struct,ui,base_index, archive, f, cr)
% UPDATE_POP_FO update population by comparing parent population and offspring population
% input:
    % pop_struct   -- struct of population matrix, problem_size etc.
    % ui           -- offspring population of pop
    % base_index   -- indices of parents which generate ui
    % archive      -- archive stores defeated parents
    % f            -- scale factor used by ui
    % cr           -- crossover rate used by ui
% output:
    % pop_struct   -- struct saves updated population
    % archive_fo   -- defeated offspring in ui_fr
    % archive      -- updated archive
    % suc_f        -- successful scale factor by ui
    % suc_cr       -- successful crossover rate by ui
    % both suc_f and suc_cr are column vector (suc_num * 1)
    % delta_k      -- improvement of successful offspring compared to
    % parent considering both fitness and conv

% Steps:
    % 1) rank individuals in pop_fo based on solely fitness
    % 2) better fitness wins
    % Version 1.4 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

    %%
    suc_cr = [];
    suc_f = [];
    archive_fo = [];
    delta_k = [];
    popsize = pop_struct.popsize;
    lambda = pop_struct.lambda; 
    pop = pop_struct.pop;

    % compare between offspring individual and parent
    for k = 1:lambda
        cur_par = pop(base_index(k), :);
        cur_off = ui(k, :);

        par_fit = cur_par(end - 1);
        off_fit = cur_off(end - 1);

        % TODO 这里可以修改，因为这个地方选入suc_f 和 suc_cr 中的F和CR可能都是来自同一个父代的子代，这些子代没有进行选拔比较
        % individual which has better fitness survive
        if par_fit > off_fit
            % parent is defeated
            [pop,archive,suc_f,suc_cr,delta_k] = replace_record(pop, base_index(k), cur_off, archive, cur_par, suc_f, suc_cr, f, cr, delta_k);
        else
            % offspring is defeated
            archive_fo = [archive_fo; cur_off];
        end
    end
    pop_struct.pop = pop; % rewrite updated population into pop_struct
end