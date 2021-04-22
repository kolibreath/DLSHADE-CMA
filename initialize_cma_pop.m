function [pop_struct,nfes] = initialize_cma_pop(xmean,sigma,problem_size, nfes, func,lu)
    lambda = 4 + floor(3 * log(problem_size));
    popsize = floor(lambda / 2);
    B = eye(problem_size, problem_size);
    D = ones(problem_size, 1);
    C = B * diag(D.^2) * B';
    invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
    eigeneval = nfes; % track update of B and D

    pop = zeros(popsize, problem_size);
    k = 1;

    while k < popsize
        pop(k, :) = (xmean' + sigma * B * (D .* randn(problem_size, 1)))';
        k = k + 1;
    end

    pop(k, :) = xmean;
    % r0此时就是选择所有个体
    r0 = 1:popsize;
    pop = boundConstraint(pop,pop,r0,lu);
    pop = evalpop(pop, func);
    nfes = nfes + popsize;
    pop_struct = assem_pop(pop, popsize, lambda, problem_size, C, D, B, ...
        invsqrtC, eigeneval, xmean, sigma);
end