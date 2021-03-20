%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB/Octave source code of L-SHADE which is an improved version of SHADE 1.1.
%% Note that this source code is transferred from the C++ source code version.
%% About L-SHADE, please see following papers:
%%
%% * Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.
%%
%% For this package, we downloaded JADE's source code from Dr. Q. Zhang's website (http://dces.essex.ac.uk/staff/qzhang) and modified it.
%%
%% Update
%% 9/Oct/2014: incompatibilities-fix between Octave and Matlab and bug-fix of population size reduction procedure (thanks to Dr. Elsayed)
%%%%%%%%%%%%%%%%%%%

clc;
clear all;

format long;
format compact;

Dimension_size = [10, 30, 50, 100];

rand('seed', sum(100 * clock));

val_2_reach = 10^(-8);
max_region = 100.0;
min_region = -100.0;
fhd = @cec14_func;

for func = 1:28
    optimum = func * 100.0;

    %% Record the best results
    outcome = [];

    fprintf('\n-------------------------------------------------------\n')
    fprintf('Function = %d, Dimension size = %d\n', func, problem_size)

    for problem_size = Dimension_size
      
      %% for each run:
      for run_id = 1:51
          
        lu = decision_range(func, problem_size)';
        max_nfes = 10000 * problem_size;
        
        %% PARAMETER SETTINGS FOR FROFI
        %% PARAMETER SETTINGS FOR EPSILON CONSTRAINTS
        epsilon = 0;
        %% PARAMETER SETTINGS FOR LSHADE
        p_best_rate = 0.11;
        arc_rate = 1.4;
        memory_size = 5;
        popsize = 18 * problem_size;
        
        %% PARAMETER SETTINGS FOR LINEAR POPULATION SIZE REDUCTION (LSPR)
        max_popsize = popsize;
        min_popsize = 4.0;
        
        %% PARAMETER SETTINGS FOR COVARIANCE ADAPTATION MATIRX (CMA)
        xmean = lu / 2;             % mean value of Gaussian distribution (problem_size * 1 vector)
        sigma = 0.3;                % step size of Gaussian distribution
        %TODO CMA-ES原算法中有stop condition 之后测试一下需不需要加上！
        stopfitness = 1e-10;
        % two meanings: 
        %   1) top p% of population based on feasibility or constraint violation  in LSHADE framework
        %   2) used in CMA to update xmean and some other parameters
        mu = p_best_rate * popsize; 
        weights = log(mu + 1/2) - log(1:mu)'; % mu * 1 vector for weighted recombination
        weights = weights / sum(weights); % normalize recombination weights array
        mueff = sum(weights)^2 / sum(weights.^2); % variance-effectiveness of sum w_i x_i
        
        % Strategy parameter setting: Adaptation
        cc = (4 + mueff / problem_size) / (N + 4 + 2 * mueff / problem_size); % time constant for cumulation for C
        cs = (mueff + 2) / (problem_size + mueff + 5); % t-const for cumulation for sigma control
        c1 = 2 / ((problem_size + 1.3)^2 + mueff); % learning rate for rank-one update of C
        cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((problem_size + 2)^2 + mueff)); % and for rank-mu update
        damps = 1 + 2 * max(0, sqrt((mueff - 1) / (problem_size + 1)) - 1) + cs; % damping for sigma  usually close to 1
        
        % Initialize dynamic (internal) strategy parameters and constants
        pc = zeros(problem_size, 1); ps = zeros(problem_size, 1); % evolution paths for C and sigma
        
        % initailize for both two sub-populations pop_fr and pop_ec
        % B defines the coordinate system
        % diagonal D defines the scaling
        % covariance matrix C
        
        B = eye(problem_size, problem_size); 
        D = ones(problem_size, 1); 
        C = B * diag(D.^2) * B';
        invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
        eigeneval = 0; % track update of B and D
        
        chiN = problem_size^0.5 * (1 - 1 / (4 * problem_size) + 1 / (21 * problem_size^2)); % expectation of ||N(0,I)|| == norm(randn(N,1))

        %% Initialize the population
        %% one row in matrix pop indicates one individual ((popsize + 2 )* 1 vector)
%         popold = repmat(lu(1, :), popsize, 1) + rand(popsize, problem_size) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));
        
        % the population is divided into two sub-populations, they have
        % the same (almost the same) parameters but different contraint
        % handling techniques (and parameters related to the techniques)
        popsize_fr = floor(popsize / 2);
        popsize_ec = popsize - popsize_fr;
        max_popsize = floor(max_popsize / 2);
        min_popsize = floor(min_popsize / 2);
        
        pop_fr = zeros(popsize_fr, problem_size);
        pop_ec = zeros(popsize_ec, problem_size);
        % initialize popsize_fr
        k = 1;
        while k <= popsize_fr || k <= popsize_ec
            pop_fr(k,:) =  xmean + sigma * B * (D .* randn(problem_size, 1));
            pop_ec(k,:) =  xmean + sigma * B * (D .* randn(problem_size, 1));
            k = k + 1;
        end
        
        while k <= popsize_fr 
            pop_fr(k, :) = xmean + sigma * B * (D .* randn(problem_size, 1));
            k = k + 1;
        end
        
        while k <= popsize_ec
            pop_ec(k, :) = xmean + sigma * B * (D .* randn(problem_size, 1));
            k = k + 1;
        end
        
        % TODO update epsilon update Covariance matrix and sigma
        % assign members for subpopulation
        pop_fr = assem_pop(pop_fr, popsize_fr, problem_size, xmean,C,D,B,invsqrtC,eigeneval);
        pop_ec = assem_pop(pop_ec, popsize_ec, problem_size, xmean,C,D,B,invsqrtC,eigeneval);
        
        %evaluate both pop_fr and pop_ec
        pop_fr = evalpop(pop_fr, func);
        pop_ec = evalpop(pop_ec, func);

        nfes = 0;
        bsf_fit_var = 1e+30;
        bsf_solution = zeros(1, problem_size);

        %%%%%%%%%%%%%%%%%%%%%%%% for out
        %%% TODO 根据算法修改
%         for i = 1:popsize
% 
%             if fitness(i) < bsf_fit_var
%                 bsf_fit_var = fitness(i);
%                 bsf_solution = pop(i, :);
%             end
% 
%             if nfes > max_nfes; break; end
%         end

        %%%%%%%%%%%%%%%%%%%%%%%% for out

        memory_sf = 0.5 .* ones(memory_size, 1);
        memory_scr = 0.5 .* ones(memory_size, 1);
        memory_pos = 1;

        archive.NP = arc_rate * popsize; % the maximum size of the archive
        archive.pop = zeros(0, problem_size); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions

        %% main loop
        while nfes < max_nfes
            %  TODO test here!
%             pop = popold; % the old population becomes the current population
            
            % generate f and cr for subpopulations respectively
            [f_fr, cr_fr] = gnFCR(pop_fr.popsize);
            [f_ec, cr_ec] = gnFCR(pop_ec.popsize);
            
            ui_fr = gnOffspring(pop_fr, p_best_rate, pop_fr.popsize, f_fr, cr_fr);
            ui_ec = gnOffspring(pop_ec, p_best_rate, pop_ec.popsize, f_ec, cr_ec);

            % evaluate offspring populations of subpopulations
            ui_fr = evalpop(ui_fr, func);
            ui_ec = evalpop(ui_ec, func);
            
            nfes = nfes + popsize;

            %%%%%%%%%%%%%%%%%%%%%%%% for out
            % TODO 选择最佳的个体的逻辑需要改一下！！！
%             for i = 1:popsize
%                 nfes = nfes + 1;
% 
%                 if children_fitness(i) < bsf_fit_var
%                     bsf_fit_var = children_fitness(i);
%                     bsf_solution = ui(i, :);
%                 end
% 
%                 if nfes > max_nfes
%                     break;
%                 end
% 
%             end

            %%%%%%%%%%%%%%%%%%%%%%%% for out
            
            delta_k_fr = []; delta_k_ec = [];
            % TODO CHECK 给元素的值可能没有使用
            suc_f_fr = [];suc_f_ec = [];
            suc_cr_fr = []; suc_cr_ec = [];
            % [pop, archive_fr, archive, suc_f, suc_cr,delta_k] = update_pop_fr(pop,ui,archive,f,cr,delta_k)
            [pop_fr,archive_fr,archive,suc_f_fr,suc_cr_fr,delta_k_fr] = update_pop_fr(pop_fr,ui_fr,archive,f_fr,cr_cr,delta_k_fr);
            [pop_ec,archive_ec,archive,suc_f_ec,suc_cr_ec,delta_k_ec] = update_pop_ec(pop_ec,ui_ec,archive,f_ec_cr_ec,delta_k_ec);
            
            % combining information from subpopulation
            delta_k = [delta_k_fr;delta_k_ec];
            suc_f = [suc_f_fr;suc_f_ec];
            suc_cr = [suc_cr_fr;suc_cr_ec];
            
            
            %% update f and cr memory
            num_success_params = numel(suc_cr);

            if num_success_params > 0
                dif_val = weights_lshade(delta_k);

                %% for updating the memory of scaling factor
                memory_sf(memory_pos) = (dif_val' * (suc_f .^2 )) / (dif_val' * suc_f);

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

            %% resizethe population size of pop_ec and pop_fr
            % TODO 如果同时对连个子种群施加LSPR这样的变化是否太大了？
            pop_ec = resize_pop(max_popsize,min_popsize,pop_ec,max_nfes,nfes);
            pop_fr = resize_pop(max_popsize,min_popsize,pop_fr,max_nfes,nfes);
            
            archive.NP = round(arc_rate * (pop_ec + pop_fr));

            if size(archive.pop, 1) > archive.NP
               rndpos = randperm(size(archive.pop, 1));
               rndpos = rndpos(1:archive.NP);
               archive.pop = archive.pop(rndpos, :);
            end

        end % end of while

        bsf_error_val = bsf_fit_var - optimum;

        if bsf_error_val < val_2_reach
            bsf_error_val = 0;
        end

        fprintf('%d th run, best-so-far error value = %1.8e\n', run_id, bsf_error_val)
        outcome = [outcome bsf_error_val];
        
       end %% end 1 run
    end
   

    fprintf('\n')
    fprintf('mean error value = %1.8e, std = %1.8e\n', mean(outcome), std(outcome))
end %% end 1 function run
