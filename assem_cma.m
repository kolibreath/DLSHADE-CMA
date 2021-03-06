function cma = assem_cma(problem_size,lambda)
%ASSEM_CMA_STRUCT construct a struct stores CMA related constants
% input:
    % problem_size  -- problem size
    % lambda        -- size of offspring population
% output:
    % cma           -- construct CMA related information into cma (struct)
% Version 1.5 Author: Shi Zeyuan 734780178@qq.com Date: 2021/4/18

    mu = lambda / 2; % (number of parent) 
    weights = log(mu + 1/2) - log(1:mu); % 1 * mu  vector for weighted recombination
    weights = weights / sum(weights); 
    mu = floor(mu);
    cma.mu = mu;
    cma.weights = weights;
    mueff = sum(weights)^2 / sum(weights.^2); % variance-effectiveness of sum w_i x_i
    cma.mueff = mueff;
    
    % Strategy parameter setting: Adaptation
    cma.cc = (4 + mueff / problem_size) / (problem_size + 4 + 2 * mueff / problem_size); % time constant for cumulation for C
    cs = (mueff + 2) / (problem_size + mueff + 5); % t-const for cumulation for sigma control
    c1 = 2 / ((problem_size + 1.3)^2 + mueff); % learning rate for rank-one update of C
    cma.c1 = c1; cma.cs = cs;
    cma.cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((problem_size + 2)^2 + mueff)); % and for rank-mu update
    cma.damps = 1 + 2 * max(0, sqrt((mueff - 1) / (problem_size + 1)) - 1) + cs; % damping for sigma  usually close to 1
        
    % in the proposed algorithm, Covariance matrix adaptation process give the control 
    % of pc and ps including its update to the evolution of subpopulations
%     % Initialize dynamic (internal) strategy parameters and constants
%     cma.pc = zeros(1, problem_size);  
%     cma.ps = zeros(1, problem_size); % evolution paths for C and sigma
    
    % expectation of ||N(0,I)|| == norm(randn(N,1))
    cma.chiN = problem_size^0.5 * (1 - 1 / (4 * problem_size) + 1 / (21 * problem_size^2)); 
end

