function pop_struct = assem_pop(pop,popsize, problem_size,xmean,C,D,B,invsqrtC,eigeneval,cma,tag)
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
    % pop_struct           -- struct of pop after assembling

% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

%% 
    pop_struct.pop = pop;
    pop_struct.popsize = popsize;
    pop_struct.xmean = xmean;
    pop_struct.sigma = cma.sigma;
    pop_struct.C = C;
    pop_struct.D = D;
    pop_struct.B = B;
    pop_struct.invsqrtC = invsqrtC;
    pop_struct.eigeneval = eigeneval;
    pop_struct.problem_size = problem_size;
    pop_struct.tag = tag;
    % each subpopulation has its own evolution path
    pop_struct.pc = cma.pc;
    pop_struct.ps = cma.ps;
    
end

