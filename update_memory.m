function [memory_sf,memory_scr,memory_pos] = update_memory(suc_f,suc_cr,memory_sf,memory_scr,memory_size,memory_pos,delta_k)
%UPDATE_MEMORY update memory for scale factor F and crossover rate cr
% input:
    % suc_f                 -- successful f
    % suc_cr                -- successful cr
    % memory_sf             -- memory for f
    % memory_scr            -- memory for cr
    % memory_pos            -- position in memory
    % delta_k               -- compute dif_val combining information from
    % constraint violation and fitness
% output:
    % memory_sf             -- memory for f
    % memory_scr            -- memory for cr
    % memory_pos            -- position in memory
    num_success_params = numel(suc_cr);
    if num_success_params > 0
        dif_val = weights_lshade(delta_k);

        %% for updating the memory of scaling factor
        memory_sf(memory_pos) = (dif_val' * (suc_f.^2)) / (dif_val' * suc_f);

        %% for updating the memory of crossover rate
        %TODO 为什么会存在suc_cr 等于0 的情况?
        if max(suc_cr) == 0 || memory_scr(memory_pos) == -1
            memory_scr(memory_pos) = -1;
        else
            memory_scr(memory_pos) = (dif_val' * (suc_cr.^2)) / (dif_val' * suc_cr);
        end

        memory_pos = memory_pos + 1;

        if memory_pos > memory_size
            memory_pos = 1;
        end

    end
end

