function [f, cr] = gnFCR(popsize, memory_size, memory_sf, memory_scr)
    %GENFCR the size of population is lambda, and gnFCR will generate lambda * 1 vectors of f and crossover rate
    %  input:
    % popsize           -- the size of population
    % memory_size       -- memory size of successful F and CR record
    % memory_sf         -- memory of scale factor f
    % memory_cr         -- memory of scale factor f cr
    %  output:
    % F            -- the vector of scale factor   (vector of lambda * 1)
    % CR           -- the vector of crossover rate (vector of lambda * 1)

    % Version 1.4 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

    %%
    lambda = popsize * 2;
    mem_rand_index = ceil(memory_size * rand(lambda, 1));
    % generate mu_f and mu_cr for Cauchy and Gaussian distribution
    mu_f = memory_sf(mem_rand_index);
    mu_cr = memory_scr(mem_rand_index);

    %% for generating crossover rate
    cr = normrnd(mu_cr, 0.1);
    term_pos = find(mu_cr == -1); % whenever cr is set to terminal value
    cr(term_pos) = 0;
    cr = min(cr, 1);
    cr = max(cr, 0);

    %% for generating scaling factor
    f = mu_f + 0.1 * tan(pi * (rand(lambda, 1) - 0.5));
    pos = find(f <= 0);

    while ~isempty(pos)
        f(pos) = mu_f(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
        pos = find(f <= 0);
    end

    f = min(f, 1);
end
