function [pop_struct,cma]= update_cma(pop_struct,nfes, sigma_lu)
% UPDATE_CMA update CMA related parameters
% input:
    % pop_struct    -- population struct and population is sorted by
    % sort_ec or sort_fr
    % nfes          -- fitness evaluation 
% output: 
    % pop_struct    -- population struct with its updated memebers
    % cma           -- updated CMA struct 
% Version 1.4 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/24
%UPDATE_CMA update parameters in CMA

    pop = pop_struct.pop;
    problem_size = pop_struct.problem_size;
    pop = pop(:, 1:problem_size);
    xold = pop_struct.xmean;
    popsize = pop_struct.popsize; % mu
    lambda = pop_struct.lambda;
    % in the LSHADE framework, the LPSR strategy is implemented, so the
    % popsize(mu) and parameters in CMA related to it will be updated as well
    % the fastest way is to reinitialized the cma struct with 'new' mu
    cma = assem_cma(problem_size,lambda);

    pop_struct.xmean =  (cma.weights)* pop; 
    % update evolution paths
    mu = popsize;
   
    pop_struct.ps = (1 - cma.cs) * pop_struct.ps ...
             + sqrt(cma.cs * (2 - cma.cs) * cma.mueff) * (pop_struct.xmean - xold) ... 
             * pop_struct.invsqrtC  / pop_struct.sigma;

    hsig = sum(pop_struct.ps.^2) / (1 - (1 - cma.cs)^(2*nfes/lambda)) / pop_struct.problem_size ... 
            < 2 + 4 / (pop_struct.problem_size + 1);
    pop_struct.pc = (1 - cma.cc) * pop_struct.pc ...
            + hsig * sqrt(cma.cc * (2 - cma.cc) * cma.mueff) * (pop_struct.xmean - xold) / pop_struct.sigma;
   
    % Adapt covariance matrix C
    % mu difference vectors (fitness and conv information should be excluded)
    artmp = (1 / pop_struct.sigma) * (pop - repmat(xold, mu, 1));
    artmp = artmp';
    pop_struct.C = (1 - cma.c1 - cma.cmu) * pop_struct.C ... % regard old matrix
           + cma.c1 * (pop_struct.pc' * pop_struct.pc ... % plus rank one update
           + (1 - hsig) * cma.cc * (2 - cma.cc) * pop_struct.C) ... % minor correction if hsig==0
           + cma.cmu * artmp * diag(cma.weights) * artmp'; % plus rank mu update
       
    sigma = pop_struct.sigma * exp((cma.cs / cma.damps) * (norm(pop_struct.ps) / cma.chiN - 1));
    if sigma >= sigma_lu(1) && sigma <= sigma_lu(2)
       pop_struct.sigma = sigma;
    end
    
    
    %TODO check here how to achieve O(N^2)
    if nfes - pop_struct.eigeneval > lambda / (cma.c1 + cma.cmu) / pop_struct.problem_size / 10
        pop_struct.eigeneval  = nfes;
        pop_struct.C = triu(pop_struct.C) + triu(pop_struct.C, 1)';
        [pop_struct.B,pop_struct.D] = eig(pop_struct.C);
        pop_struct.B = pop_struct.B';
        pop_struct.D = sqrt(diag(pop_struct.D));
        pop_struct.invsqrtC = pop_struct.B * diag(pop_struct.D .^ -1) * pop_struct.B';
    end   
    
    
    
end






