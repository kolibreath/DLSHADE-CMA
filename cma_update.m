function [mean, sigma, C] = cma_update(lambda,mu, mean, sigma, C)
%CMA_UPDATE update mean, sigma (step size) and covariance matrix according
%to pure CMA-ES adaptation rule

    weights = log(mu + 1/2) - log(1:mu)'; % muXone array for weighted recombination
    mu = floor(mu);
    weights = weights / sum(weights); % normalize recombination weights array
    mueff = sum(weights)^2 / sum(weights.^2); % variance-effectiveness of sum w_i x_i

    % Strategy parameter setting: Adaptation
    cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N); % time constant for cumulation for C
    cs = (mueff + 2) / (N + mueff + 5); % t-const for cumulation for sigma control
    c1 = 2 / ((N + 1.3)^2 + mueff); % learning rate for rank-one update of C
    cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((N + 2)^2 + mueff)); % and for rank-mu update
    damps = 1 + 2 * max(0, sqrt((mueff - 1) / (N + 1)) - 1) + cs; % damping for sigma
    % usually close to 1
    % Initialize dynamic (internal) strategy parameters and constants
    pc = zeros(N, 1); ps = zeros(N, 1); % evolution paths for C and sigma
    B = eye(N, N); % B defines the coordinate system
    D = ones(N, 1); % diagonal D defines the scaling
    C = B * diag(D.^2) * B'; % covariance matrix C
    invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
    eigeneval = 0; % track update of B and D
    chiN = N^0.5 * (1 - 1 / (4 * N) + 1 / (21 * N^2)); % expectation of ||N(0,I)|| == norm(randn(N,1))
    
    
end

