function pop = update_cma(pop,cma,nfes)
% UPDATE_CMA update CMA related parameters
% input:
    % pop           -- population struct and population is sorted by
    % sort_ec or sort_fr
    % cma           -- cma struct
% output: 
    % pop           -- updated population and its memebers
% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/20
%UPDATE_CMA update parameters in CMA

    xold = pop.xmean;
    pop.xmean =  weights * pop(1:cma.mu, 1:pop.problem_size); 
    
    % update evolution paths
    pop.ps = (1 - cma.cs) * pop.ps ...
             + sqrt(cma.cs * (2 - cma.cs) * cma.mueff) * (pop.xmean - xold) ... 
             * pop.invsqrtC  / pop.sigma;

    hsig = sum(pop.ps.^2) / (1 - (1 - cma.cs)^(2*nfes/pop.popsize)) / pop.problem_size ... 
            < 2 + 4 / (pop.problem_size + 1);
    pop.pc = (1 - cma.cc) * pop.pc ...
            + hsig * sqrt(cma.cc * (2 - cma.cc) * cma.mueff) * (pop.xmean - xold) / pop.sigma;
   
    % Adapt covariance matrix C
    % mu difference vectors (fitness and conv information should be excluded)
    artmp = (1 / pop.sigma) * (pop(1:cma.mu, 1:pop.problem_size)) - repmat(xold, cma.mu, 1);
    artmp = artmp';
    %TODO 这里检查一下由原来的公式的列向量变成了行向量
    pop.C = (1 - cma.c1 - cma.cmu) * pop.C ... % regard old matrix
           + cma.c1 * (pop.pc * pop.pc' ... % plus rank one update
           + (1 - hsig) * cma.cc * (2 - cma.cc) * pop.C) ... % minor correction if hsig==0
           + cma.cmu * artmp * diag(cma.weights) * artmp'; % plus rank mu update
       
    pop.sigma = pop.sigma * exp((cma.cs / cma.damps) * (norm(pop.ps) / cma.chiN - 1));
    
    %TODO check here how to achieve O(N^2)
%     if nfes - pop.eigenval > pop.popsize / (cma.c1 + cma.cmu) / pop.problem_size / 10
%     end  
    
    pop.eigenval  = nfes;
    pop.C = triu(pop.C) + tri(pop.C, 1)';
    [pop.B,pop.D] = eig(C);
    pop.D = sqrt(diag(pop.D));
    pop.invsqrtC = pop.B * diag(D .^ -1) * pop.B';
end






