function [f, cr] = gnFCR(pop_struct,memory_size,memory_sf,memory_scr,weights)
%GENFCR 
    %  input: 
        % pop_struct        -- information of sub populations
        % memory_size       -- memory size of successful F and CR record
        % memory_sf         -- memory of scale factor f
        % memory_cr         -- memory of scale factor f cr
    %  output:
        % F            -- the vector of scale factor   (vector of offpopsize * 1)
        % CR           -- the vector of crossover rate (vector of offpopsize * 1)

% Version 1.4 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18 

%%
%      mem_rand_index = ceil(memory_size * rand(offpopsize, 1));
     offpopsize = pop_struct.lambda;
     temp_rand = rand(offpopsize,1);
     % simple roulette-wheel selection of mean values of distributions 
     mem_rand_index = zeros(offpopsize,1);
     for i = 1 : offpopsize
         for j = 1 : memory_size
             if temp_rand(i) < weights(j)
                 mem_rand_index(i) = j;
                 break;
             end
         end
     end
     % generate mu_f and mu_cr for Cauchy and Gaussian distribution
     mu_f = memory_sf(mem_rand_index);
     mu_cr = memory_scr(mem_rand_index);

     %% for generating crossover rate
     cr = normrnd(mu_cr, 0.1);
     term_pos = find(mu_cr == -1);       % whenever cr is set to terminal value
     cr(term_pos) = 0;
     cr = min(cr, 1);
     cr = max(cr, 0);
      
     %% for generating scaling factor
     f = mu_f + 0.1 * tan(pi * (rand(offpopsize, 1) - 0.5));
     pos = find(f <= 0);

     while ~isempty(pos)
         f(pos) = mu_f(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
         pos = find(f <= 0);
    end

    f = min(f, 1);
end

