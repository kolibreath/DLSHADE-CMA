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

problem_size = 10;
max_nfes = 10000 * problem_size;

rand('seed', sum(100 * clock));

val_2_reach = 10^(-8);
max_region = 100.0;
min_region = -100.0;
lu = [-100 * ones(1, problem_size); 100 * ones(1, problem_size)];
fhd = @cec14_func;

for func = 1:30
    optimum = func * 100.0;

    %% Record the best results
    outcome = [];

    fprintf('\n-------------------------------------------------------\n')
    fprintf('Function = %d, Dimension size = %d\n', func, problem_size)

    for run_id = 1:51
        %%  parameter settings for L-SHADE
        p_best_rate = 0.11;
        arc_rate = 1.4;
        memory_size = 5;
        popsize = 18 * problem_size;

        max_popsize = popsize;
        min_popsize = 4.0;

        %% Initialize the main population
        popold = repmat(lu(1, :), popsize, 1) + rand(popsize, problem_size) .* (repmat(lu(2, :) - lu(1, :), popsize, 1));
        pop = popold; % the old population becomes the current population

        fitness = feval(fhd, pop', func);
        fitness = fitness';

        nfes = 0;
        bsf_fit_var = 1e+30;
        bsf_solution = zeros(1, problem_size);

        %%%%%%%%%%%%%%%%%%%%%%%% for out
        for i = 1:popsize
            nfes = nfes + 1;

            if fitness(i) < bsf_fit_var
                bsf_fit_var = fitness(i);
                bsf_solution = pop(i, :);
            end

            if nfes > max_nfes; break; end
        end

        %%%%%%%%%%%%%%%%%%%%%%%% for out

        memory_sf = 0.5 .* ones(memory_size, 1);
        memory_scr = 0.5 .* ones(memory_size, 1);
        memory_pos = 1;

        archive.NP = arc_rate * popsize; % the maximum size of the archive
        archive.pop = zeros(0, problem_size); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions

        %% main loop
        while nfes < max_nfes
            pop = popold; % the old population becomes the current population
            [temp_fit, sorted_index] = sort(fitness, 'ascend');

            mem_rand_index = ceil(memory_size * rand(popsize, 1));
            % generate mu_f and mu_cr for Cauchy and Gaussian distribution
            mu_f = memory_sf(mem_rand_index);
            mu_cr = memory_scr(mem_rand_index);

            [f, cr] = gnFCR(popsize);

            r0 = [1:popsize];
            popAll = [pop; archive.pop];
            [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);

            %% DE/current-to-pbest mutation and crossover
            pNP = max(round(p_best_rate * popsize), 2); %% choose at least two best solutions
            randindex = ceil(rand(1, popsize) .* pNP); %% select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions

            % mutation
            vi = pop + f(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
            vi = boundConstraint(vi, pop, lu);

            % crossover
            mask = rand(popsize, problem_size) > cr(:, ones(1, problem_size)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1:popsize)'; cols = floor(rand(popsize, 1) * problem_size) + 1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize problem_size], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop(mask);

            children_fitness = feval(fhd, ui', func);
            children_fitness = children_fitness';

            %%%%%%%%%%%%%%%%%%%%%%%% for out
            for i = 1:popsize
                nfes = nfes + 1;

                if children_fitness(i) < bsf_fit_var
                    bsf_fit_var = children_fitness(i);
                    bsf_solution = ui(i, :);
                end

                if nfes > max_nfes
                    break;
                end

            end

            %%%%%%%%%%%%%%%%%%%%%%%% for out

            dif = abs(fitness - children_fitness);

            %% I == 1: the parent is better; I == 2: the offspring is better
            % TODO wtf???
            I = (fitness > children_fitness);
            goodCR = cr(I == 1);
            goodF = f(I == 1);
            dif_val = dif(I == 1);

            archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));

            [fitness, I] = min([fitness, children_fitness], [], 2);

            popold = pop;
            popold(I == 2, :) = ui(I == 2, :);

            num_success_params = numel(goodCR);

            if num_success_params > 0
                sum_dif = sum(dif_val);
                dif_val = dif_val / sum_dif;

                %% for updating the memory of scaling factor
                memory_sf(memory_pos) = (dif_val' * (goodF.^2)) / (dif_val' * goodF);

                %% for updating the memory of crossover rate
                if max(goodCR) == 0 || memory_scr(memory_pos) == -1
                    memory_scr(memory_pos) = -1;
                else
                    memory_scr(memory_pos) = (dif_val' * (goodCR.^2)) / (dif_val' * goodCR);
                end

                memory_pos = memory_pos + 1;

                if memory_pos > memory_size
                    memory_pos = 1;
                end

            end

            %% for resizing the population size
            plan_popsize = round((((min_popsize - max_popsize) / max_nfes) * nfes) + max_popsize);

            if popsize > plan_popsize
                reduction_ind_num = popsize - plan_popsize;

                if popsize - reduction_ind_num < min_popsize
                    reduction_ind_num = popsize - min_popsize;
                end

                popsize = popsize - reduction_ind_num;

                for r = 1:reduction_ind_num
                    [valBest indBest] = sort(fitness, 'ascend');
                    worst_ind = indBest(end);
                    popold(worst_ind, :) = [];
                    pop(worst_ind, :) = [];
                    fitness(worst_ind, :) = [];
                end

                archive.NP = round(arc_rate * popsize);

                if size(archive.pop, 1) > archive.NP
                    rndpos = randperm(size(archive.pop, 1));
                    rndpos = rndpos(1:archive.NP);
                    archive.pop = archive.pop(rndpos, :);
                end

            end

        end % end of while

        bsf_error_val = bsf_fit_var - optimum;

        if bsf_error_val < val_2_reach
            bsf_error_val = 0;
        end

        fprintf('%d th run, best-so-far error value = %1.8e\n', run_id, bsf_error_val)
        outcome = [outcome bsf_error_val];
    end %% end 1 run

    fprintf('\n')
    fprintf('mean error value = %1.8e, std = %1.8e\n', mean(outcome), std(outcome))
end %% end 1 function run
