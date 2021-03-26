function pop_struct = assem_pop(pop,popsize, problem_size,C,D,B,invsqrtC,eigeneval,xmean,sigma)
%ASSEM_POP assign members for subpopulation
% input:    
    % pop           -- population (pop_fr or pop_ec)
    % popsize       -- population size
    % problem_size  -- problem size
    % C,D,B         -- CMA-related parameters
    % invsqrtC      -- C ^ (- 1/2)
    % eigeneval     -- eigene value
    % xmean         -- mean value of distribution
    % sigma         -- step  size

% output:
    % pop_struct           -- struct of pop after assembling

% Version 1.4 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

%% 
    pop_struct.pop = pop;
    pop_struct.popsize = popsize;  % refers to mu in ES
    pop_struct.xmean = xmean;
    pop_struct.sigma = sigma;
    pop_struct.C = C;
    pop_struct.D = D;
    pop_struct.B = B;
    pop_struct.invsqrtC = invsqrtC;
    pop_struct.eigeneval = eigeneval;
    pop_struct.problem_size = problem_size;
    % each subpopulation has its own evolution path
    pop_struct.pc = zeros(1,problem_size);
    pop_struct.ps = zeros(1, problem_size);
    
end

