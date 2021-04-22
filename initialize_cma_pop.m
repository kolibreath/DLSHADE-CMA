function [pop_struct,nfes] = initialize_cma_pop(xmean,problem_size, nfes, func)
    lambda = 4 + floor(3 * log(problem_size));
    popsize = floor(lambda / 2);
    B = eye(problem_size, problem_size);
    D = ones(problem_size, 1);
    C = B * diag(D.^2) * B';
    invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
    eigeneval = nfes; % track update of B and D
    sigma = 0.3;

    pop = zeros(popsize, problem_size);
    k = 1;

    while k < popsize
        pop(k, :) = (xmean' + sigma * B * (D .* randn(problem_size, 1)))';
        k = k + 1;
    end

    pop(k, :) = xmean;

    pop = evalpop(pop, func);
    nfes = nfes + popsize;
    pop_struct = assem_pop(pop, popsize, lambda, problem_size, C, D, B, ...
        invsqrtC, eigeneval, xmean, sigma);
end