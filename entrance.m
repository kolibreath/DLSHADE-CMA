%%%%%%%%%%%%%%%%%%%
%% Version 1.9 
%% 1. Global Search Stage
%% 2. Local Search Stage
%%%%%%%%%%%%%%%%%%%


clc;
clear all;

format long;
format compact;

global initial_flag; % for CEC2017 test function evaluation

rand('seed', sum(100 * clock));

for func = 1:1

    optimum = func * 100.0;
    %% PARAMETER SETTINGS FOR PROBLEM SIZE
    Dimension_size = [10, 30, 50, 100];
    fprintf('\n-------------------------------------------------------\n')

    %% for each problem size
    for dim_num = 1:1

        problem_size = Dimension_size(dim_num);
        initial_flag = 0;
        fprintf('Function = %d, Dimension size = %d\n', func, problem_size)

        %% for each run:
        for run_id = 2:2

            outcome = [];

            bsf_solution = zeros(problem_size + 4, 1);
            bsf_solution(end - 1) = inf;
            bsf_solution(end) = inf;

            last_bsf_solution = bsf_solution;

            lu = decision_range(func, problem_size)'; % 2 * problem_size matrix
            max_nfes = 10000 * problem_size;
            nfes = 0;

            %% initialize population
            global_popsize = problem_size * 10; 
            global_pop = repmat(lu(1, :), global_popsize, 1) +  ...
            rand(global_popsize, problem_size) .* (repmat(lu(2, :) - lu(1, :), global_popsize, 1));
            global_pop = evalpop(global_pop,func);
            
            global_pop_struct = []; 
            global_pop_struct.popsize = global_popsize; 
            global_pop_struct.problem_size = problem_size;
            global_pop_struct.pop = global_pop;
            global_pop_struct.lambda = global_popsize;
            
              %% PARAMETER SETTINGS FOR SHADE
            p_best_rate = 0.11;
            arc_rate = 1.4; % archive for saving defeated parents
            memory_size = 5; % memory for successful F and CR

            memory_sf = 0.5 .* ones(memory_size, 1);
            memory_scr = 0.5 .* ones(memory_size, 1);
            memory_pos = 1;

            archive.NP = arc_rate * global_popsize; % the maximum size of the archive
            archive.pop = zeros(0, problem_size); % the solutions stored in te archive
            archive.funvalues = zeros(0, 1); % the function value of the archived solutions

            % pos_entropy = []; % record for positional entropy
            % 使用位置熵的概念先看看各个函数的性质 然后定量设置一个值
            while nfes < floor(max_nfes * 0.75)
                %% 这里的检查一下权重设置有无问题
%                 memory_weight = ones(memory_size,1) / memory_size;
                memory_weight = [0.2 0.4 0.6 0.8 1.0];
                %% generate f and cr for subpopulations respectively
                [f, cr] = gnFCR(global_pop_struct, memory_size, memory_sf, memory_scr, memory_weight);
                % Note: ui is un-evaluated matrix (lambda * problem_size)
                [ui, base] = gnOffspring(global_pop_struct, lu, archive, nfes, max_nfes, f, cr,1);
                
                %% evaluate offspring populations of subpopulations
                ui = evalpop(ui, func);
                nfes = nfes + global_pop_struct.popsize;

                % updated population stored in structs
                % evaluate population using feasibility rule
                [global_pop_struct, ~, archive, suc_f, suc_cr, delta_k] = update_pop_fr(global_pop_struct, ui, base, archive, f, cr);
                %% update f and cr memory
                [memory_sf, memory_scr, memory_pos] = update_memory(suc_f, suc_cr, memory_sf, memory_scr, memory_size, memory_pos, delta_k);

                %% resize archive
                archive.NP = round(arc_rate * global_pop_struct.popsize);

                if size(archive.pop, 1) > archive.NP
                    rndpos = randperm(size(archive.pop, 1));
                    rndpos = rndpos(1:archive.NP);
                    archive.pop = archive.pop(rndpos, :);
                end

                %% update best so far solution
                bsf_solution = find_bsf(global_pop_struct, global_pop_struct, bsf_solution);
            end % end of while
            
            % 计算位置熵 通过fitness 和 conv 排序之后的位置之间的差值
            [~,columns] = size(global_pop_struct.pop);
            [~, sortedfit_index] = sortrows(global_pop_struct.pop, columns - 1);
            [~, sortedcon_index] = sortrows(global_pop_struct.pop, columns);
            pop_entropy = sum(abs(sortedfit_index - sortedcon_index));

            % fprintf('run= %d, fitness = %d\n, conv = %d\n', run_id, bsf_solution(end - 1), bsf_solution(end));
            fprintf('func= %d, dimension = %d, entropy = %d\n', func, problem_size, pop_entropy);
        end %% end 1 run

        fprintf("---------------------------------------------------\n");
        fprintf('fitness = %d\n, conv = %d\n', bsf_solution(end - 1), bsf_solution(end));
        plot(outcome);

    end %% end of iterate one problem size

    fprintf('\n')
end %% end 1 function run
