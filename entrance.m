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
    for dim_num = 2:2

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
            global_popsize = problem_size * 5; 
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
            %%  ----------------------------- Global Search Stage --------------------------------
            while nfes < floor(max_nfes * 0.75)
                %% generate f and cr for subpopulations respectively
                [f, cr] = gnFCR(global_pop_struct, memory_size, memory_sf, memory_scr);
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
                bsf_solution = find_bsf(global_pop_struct, bsf_solution);
            end % end of while
            
            %% ------------------------------ Local Search Stage --------------------------------
            %% INITIALIZATION FOR CMA
            lambda = 4 + floor(3 * log(problem_size));
            popsize = floor(lambda / 2);
            B = eye(problem_size, problem_size);
            D = ones(problem_size, 1);
            C = B * diag(D.^2) * B';
            invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
            eigeneval = 0; % track update of B and D
            xmean = rand(1, problem_size) .* (lu(2) - lu(1));
            sigma = 0.3;

            %% PARAMETER SETTINGS FOR COVARIANCE ADAPTATION MATIRX (CMA)
            cma = assem_cma(problem_size, lambda);
            sigma_lu = [1e-20, min((lu(2) - lu(1)) / 2)];
            %% kmeans
            cluster_number = 8;
            idx = kmeans(global_pop_struct.pop(:,1:problem_size), cluster_number);
            [pop_array, nfes] = repair(idx, cluster_number, global_pop_struct.pop, ... 
                                global_pop_struct.problem_size, mu, lambda, nfes, func);
            % 先完成一个没有子种群交换的版本
            % TODO 一定要修改evalpop 修改成会自动改变nfes的版本... 太傻比了
            while nfes < max_nfes
                % 每个子种群自行迭代 更新shade中的memory archive等等
                for i = 1 : cluster_number
                    pop_struct = pop_array(i);
                    [f, cr] = gnFCR(pop_struct, memory_size, memory_sf, memory_scr);
                    [ui, base] = gnOffspring(pop_struct, lu, archive, nfes, max_nfes, f, cr,2);
                    ui = evalpop(ui, func);
                    nfes = nfes + pop_struct.lambda;

                    [pop_struct, archive, archive, suc_f, suc_cr, delta_k] = ...
                     update_pop_fr(pop_struct, ui, base, archive, f, cr);

                    [memory_sf, memory_scr, memory_pos] = update_memory(suc_f, suc_cr, memory_sf, ... 
                    memory_scr, memory_size, memory_pos, delta_k);

                    %% resize the population size of pop_fo and pop_fr
                    archive.NP = round(arc_rate * lambda);

                    if size(archive.pop, 1) > archive.NP
                        rndpos = randperm(size(archive.pop, 1));
                        rndpos = rndpos(1:archive.NP);
                        archive.pop = archive.pop(rndpos, :);
                    end

                    % CMA parameters update (populations in pop_fr and pop_fo are sorted)
                    [pop_struct] = update_cma(pop_struct, nfes, sigma_lu);
                    bsf_solution = find_bsf(pop_struct, bsf_solution);
                end 
            end
            % fprintf('run= %d, fitness = %d\n, conv = %d\n', run_id, bsf_solution(end - 1), bsf_solution(end));
           
        end %% end 1 run

        fprintf("---------------------------------------------------\n");
        fprintf('fitness = %d\n, conv = %d\n', bsf_solution(end - 1), bsf_solution(end));
        plot(outcome);

    end %% end of iterate one problem size

    fprintf('\n')
end %% end 1 function run
