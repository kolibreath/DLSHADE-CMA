function pop = assem_pop(pop,popsize,xmean,C,D,B,invsqrtC,eigeneval)
%ASSEM_POP assign members for subpopulation
% input:    
    % pop           -- population (pop_fr or pop_ec)
    % popsize       -- population size
    % xmean         -- mean value of the distribution
    % C,D,B         -- CMA-related parameters
    % invsqrtC      -- C ^ (- 1/2)
    % eigeneval     -- eigene value

    pop.popsize = popsize;
    pop.xmean = xmean;
    pop.C = C;
    pop.D = D;
    pop.B = B;
    pop.invsqrtC = invsqrtC;
    pop.eigeneval = eigeneval;
    
end
