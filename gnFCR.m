function [f, cr] = gnFCR(popsize)
%GENFCR generate popsize * 1 F and CR vector
    %  input: 
        % popsize      -- the size of population
    %  output:
        % F            -- the vector of scale factor   (vector of popsize * 1)
        % CR           -- the vector of crossover rate (vector of popsize * 1)

    
     %% for generating crossover rate
     cr = normrnd(mu_cr, 0.1);
     term_pos = find(mu_cr == -1);       % whenever cr is set to terminal value
     cr(term_pos) = 0;
     cr = min(cr, 1);
     cr = max(cr, 0);
      
     %% for generating scaling factor
     % TODO  should be randc in original paper
     f = mu_f + 0.1 * tan(pi * (rand(popsize, 1) - 0.5));
     pos = find(f <= 0);

     while ~isempty(pos)
         f(pos) = mu_f(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
         pos = find(f <= 0);
    end

    f = min(f, 1);
end

