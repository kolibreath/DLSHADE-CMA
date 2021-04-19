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

global initial_flag ; % for CEC2017 test function evaluation

rand('seed', sum(100 * clock));

for func = 5:5

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
        bsf_solution(end-1) = inf;
        bsf_solution(end)   = inf;
      
        last_bsf_solution = bsf_solution;
        
        lu = decision_range(func, problem_size)';  % 2 * problem_size matrix
        max_nfes = 10000* problem_size;
        nfes = 0;
        
        %% PARAMETER SETTINGS FOR FROFI
        %% PARAMETER SETTINGS FOR SHADE 
        p_best_rate = 0.11;
        arc_rate = 1.4;         % archive for saving defeated parents
        memory_size = 5;        % memory for successful F and CR
                  
        %% PARAMETER SETTTINGS FOR COVARIANCE MATRIX ADAPTATION 
        % Note: subpopulations update these parameters seperately
        % B defines the coordinate system
        % diagonal D scales the coordinate system
        % covariance matrix C =  B * diag(D.^2) * B'; 
        lambda = 4 + floor(3*log(problem_size)); 
        popsize = floor(lambda / 2);
        B = eye(problem_size, problem_size); 
        D = ones(problem_size, 1); 
        C = B * diag(D.^2) * B';
        invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
        eigeneval = 0; % track update of B and D
        xmean = rand(1,problem_size) .* (lu(2)-lu(1));
        sigma = 0.3;
        
        %% PARAMETER SETTINGS FOR COVARIANCE ADAPTATION MATIRX (CMA)
        cma = assem_cma(problem_size,lambda);
        sigma_lu = [1e-20, min((lu(2)-lu(1))/2)];
      

        %% INTIIALIZE THE POPULATION
        % individual in population will be viewed as (problem_size + 4) * 1 vector
        % the additional 4 'elements' are [g, h, f, conV] 
   
        % the population is divided into two sub-populations, they have
        % the same (almost the same) parameters, while pop_fr using feasibility rule
        % and pop_fo using only fitness to compare between two individuals
        
        popsize_fr = popsize;
        popsize_fo = popsize;
        
        pop_fr = zeros(popsize_fr, problem_size);
        pop_fo = zeros(popsize_fo, problem_size);
        
        %% initialize both subpopulations        
        k = 1;
        while k <= popsize_fr &&  k <= popsize_fo
            pop_fr(k,:) =  (xmean' + sigma * B * (D .* randn(problem_size, 1)))';
            pop_fo(k,:) =  (xmean' + sigma * B * (D .* randn(problem_size, 1)))';
            k = k + 1;
        end
        
        while k <= popsize_fr 
            pop_fr(k, :) = (xmean' + sigma * B * (D .* randn(problem_size, 1)))';
            k = k + 1;
        end
        
        while k <= popsize_fo
            pop_fo(k, :) = (xmean' + sigma * B * (D .* randn(problem_size, 1)))';
            k = k + 1;
        end
        
        clear k; 
        
        % assign members for subpopulation
        pop_fr_struct = assem_pop(pop_fr,popsize_fr,lambda,problem_size,C,D,B,invsqrtC,eigeneval,xmean,sigma);
        pop_fo_struct = assem_pop(pop_fo,popsize_fo,lambda,problem_size,C,D,B,invsqrtC,eigeneval,xmean,sigma);
        
        %% evaluate both pop_fr and pop_fo
        % TODO test! 生成的个体都是可行解！？
        pop_fr = evalpop(pop_fr, func);
        pop_fo = evalpop(pop_fo, func);
        
        pop_fr_struct.pop = pop_fr;
        pop_fo_struct.pop = pop_fo;
        
        nfes = nfes + popsize;
        % TODO 测试种群合适的大小 解决协方差矩阵的问题
        %% INITIALIZATION FOR LSHADE ARCHIVE
        memory_sf = 0.5 .* ones(memory_size, 1);
        memory_scr = 0.5 .* ones(memory_size, 1);
        memory_pos = 1;

        archive.NP = arc_rate * popsize; % the maximum size of the archive
        archive.pop = zeros(0, problem_size); % the solutions stored in te archive
        archive.funvalues = zeros(0, 1); % the function value of the archived solutions

        % there is no need for pop_fr and pop_fo outside of their structs, delete them LOL
        clear pop_fo;
        clear pop_fr;
        
        clear cma;

        %% -------------------------------- main loop -------------------------------
        sigma_gen = 20; 
        gen = 0; % after every sigma_gen generations, output the mean of sigma 
        bsf_index = 0;
        bsf_gen_len = 10 + ceil(30*problem_size/lambda);
        bsf_unchange_counter = 0;
        sigma_record_fr = [];
        sigma_record_fo = [];
        counter = 0;
        while nfes < max_nfes
            memory_w_fr = memory_weights(sigma_record_fr, sigma_gen, memory_size, memory_sf);
            memory_w_fo = memory_weights(sigma_record_fo, sigma_gen, memory_size, memory_sf);
            
            %% generate f and cr for subpopulations respectively
            [f_fr, cr_fr] = gnFCR(pop_fr_struct,memory_size,memory_sf,memory_scr,memory_w_fr);
            [f_fo, cr_fo] = gnFCR(pop_fo_struct,memory_size,memory_sf,memory_scr,memory_w_fo);
            
            % Note: ui_fr and ui_fo are un-evaluated matrix (lambda * problem_size)
            [ui_fr,base_fr] = gnOffspring(pop_fr_struct,lu,archive,nfes,max_nfes,f_fr,cr_fr);
            [ui_fo,base_fo] = gnOffspring(pop_fo_struct,lu,archive,nfes,max_nfes,f_fo,cr_fo);

           
            
            %% evaluate offspring populations of subpopulations
            ui_fr = evalpop(ui_fr, func);
            ui_fo = evalpop(ui_fo, func);
            
            % TODO nfes 是否有重复计算 算多了
            nfes = nfes + pop_fr_struct.lambda + pop_fo_struct.lambda;
            
            % updated subpopulations stored in structs
            [pop_fr_struct,archive_fr,archive,suc_f_fr,suc_cr_fr,delta_k_fr] = update_pop_fr(pop_fr_struct,ui_fr,base_fr,archive,f_fr,cr_fr);
            [pop_fo_struct,archive_fo,archive,suc_f_fo,suc_cr_fo,delta_k_fo] = update_pop_fo(pop_fo_struct,ui_fo,base_fo,archive,f_fo,cr_fo);
            
            % combining information from subpopulation
            delta_k = [delta_k_fr;delta_k_fo];
            suc_f = [suc_f_fr;suc_f_fo];
            suc_cr = [suc_cr_fr;suc_cr_fo];
            
            %TODO remove useless variables
            %% update f and cr memory
            [memory_sf,memory_scr,memory_pos] = update_memory(suc_f,suc_cr,memory_sf,memory_scr,memory_size,memory_pos,delta_k);

            %% resize the population size of pop_fo and pop_fr
           
            % --- Note: population in structs are sorted ---
            [pop_fr_struct,pop_fo_struct,delete_individuald] ...
                = subpop_com(pop_fr_struct,pop_fo_struct, ...
                  archive_fr,archive_fo);

            archive.pop = [archive.pop; delete_individuald];
            archive.NP = round(arc_rate * (pop_fo_struct.popsize + pop_fr_struct.popsize));

            if size(archive.pop, 1) > archive.NP
               rndpos = randperm(size(archive.pop, 1));
               rndpos = rndpos(1:archive.NP);
               archive.pop = archive.pop(rndpos, :);
            end

            % CMA parameters update (populations in pop_fr and pop_fo are sorted)
            [pop_fr_struct] = update_cma(pop_fr_struct,nfes,sigma_lu);
            [pop_fo_struct] = update_cma(pop_fo_struct,nfes,sigma_lu);
            
            if max(pop_fr_struct.D) > 1e7 * min(pop_fr_struct.D) || max(pop_fo_struct.D) > 1e7 * min(pop_fo_struct.D)
               break;
            end
            
           %% update best so far solution
           bsf_solution = find_bsf(pop_fr_struct,pop_fo_struct,bsf_solution);
           
           %% update stop_trigger      
           temp = floor(gen / bsf_gen_len);
           % update last_bsf_solution after every bsf_len generation
           if temp > bsf_index
               last_bsf_solution = bsf_solution;
               if all(last_bsf_solution == bsf_solution)
                   bsf_unchange_counter = bsf_unchange_counter + 1;
               end
               bsf_index = temp;
           end
           % restart only happens at the late stage of the search
           if nfes >= max_nfes * 0.8
                [restart_index_fr, restart_index_fo] = stop_trigger(bsf_unchange_counter,pop_fr_struct,pop_fo_struct);
                if restart_index_fr == 1
                    [pop_fr_struct] = restart_pop(pop_fr_struct,func);
                    nfes = nfes + pop_fr_struct.popsize;
                end
                if restart_index_fo == 1
                    [pop_fo_struct] = restart_pop(pop_fo_struct,func);
                    nfes = nfes + pop_fo_struct.popsize;
                end
            end
           
           %% update sigma record
           sigma_record_fr = [sigma_record_fr; pop_fr_struct.sigma];
           sigma_record_fo = [sigma_record_fo; pop_fo_struct.sigma];
           
           outcome = [outcome;[gen,bsf_solution(end-1)]];
           gen = gen + 1;
           
           outcome = [outcome; [gen, bsf_solution(end-1)] ];
        end % end of while
        fprintf('run= %d, fitness = %d\n, conv = %d\n' ,run_id,bsf_solution(end-1),bsf_solution(end));
        pplot(outcome);
       end %% end 1 run
       
       fprintf("---------------------------------------------------\n");
       fprintf('fitness = %d\n, conv = %d\n' ,bsf_solution(end-1),bsf_solution(end));
       plot(outcome);
    
    end %% end of iterate one problem size
   

    fprintf('\n')
end %% end 1 function run
