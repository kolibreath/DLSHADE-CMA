function pop = assem_pop(pop,popsize, problem_size,xmean,C,D,B,invsqrtC,eigeneval,cma,tag)
%ASSEM_POP assign members for subpopulation
% input:    
    % pop           -- population (pop_fr or pop_ec)
    % popsize       -- population size
    % problem_size  -- problem size
    % xmean         -- mean value of the distribution
    % C,D,B         -- CMA-related parameters
    % invsqrtC      -- C ^ (- 1/2)
    % eigeneval     -- eigene value
    % cma           -- cma struct for Covariance Matrix Adaptation parameters
    % tag           -- tag == 1 -> pop_fr; tag == 2 -> pop_ec
% output:
    % pop           -- struct of pop after assembling

% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

%% 
    pop.popsize = popsize;
    pop.xmean = xmean;
    pop.sigma = cma.sigma;
    pop.C = C;
    pop.D = D;
    pop.B = B;
    pop.invsqrtC = invsqrtC;
    pop.eigeneval = eigeneval;
    pop.problem_size = problem_size;
    pop.tag = tag;
    % each subpopulation has its own evolution path
    pop.pc = cma.pc;
    pop.ps = cma.ps;
    
end

