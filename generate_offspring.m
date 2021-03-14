function u = generate_offspring(pop, par_archive, i , N,max_fes,fes,p ,F, CR, xmean, sigma, B, D)
% generate_offspring generating new offspring by DE/current-to-pbetter*/r
% as mutation strategy and implementing binomial crossover as crossover
% strategy

% DE/current-to-pbetter*/r  
% v = xi + F_w(x_pbetter - xi) + F(x_r1* - x_r2*); r1 is selected from
% popluation and r2 is selected from the union of popluation and parent
% archive and then they are ranked by a ranked selection manner
% input parameters
  % pop             -- popluation 
  % par_archive     --  where the index r2 is selected from
  % i               -- index of current individual 
  % N               -- dimension size
  % max_fes         -- maximum of fes
  % p               -- one of top p * popsize individual will be selectedas pbetter
  % fes             -- current fes
  % F               -- scale factor of DE
  % CR              -- crossover factor of DE
  % xmean           -- mean value of population
  % sigma           -- step size of the distribution
  % B, D            -- results of eigendecomposition
  
% output paramters 
  % u  -- the individual (new offspring) after mutation and crossover  procedure
  
  [popsize, rows] = size(pop);   % rows variable indicates the size of dimension
  [par_archive_size, ~] = size(par_archive);
  
  x_pbetter = pop(randi(floor(popsize * p) + 1), :);
  
  % find F_w
  % TODO 这里的F_w 是否需要重置到 [0, 1]之间？?? 这里写了重置！
  if fes <= 0.2 * max_fes
      F_w = 0.7 * F;
  elseif fes > 0.2 * max_fes && fes <= 0.4 * max_fes
      F_w = 0.8 * F;
  else
      F_w = 1.2 * F;
  end 
  
  F_w = min(1, F_w); % place F_w in the upper bound 
  
  % select r1 and r2. 
  % xi , xr1 and xr2 should be different. and xr1 will be better than xr2 
  % which has less constraint violation or has a better fitness
  
  
  % rank selection 
  % TODO 设置比例什么时候选择以fitness 什么时候选择constraint violation
  
  gamma = 0.5;   % sort population by fitness or conV ; minimized both fitness and conv
  if rand < gamma
      % sort  population based on fitness 
      pop = sortrows(pop, rows - 1);
  else
     % sort  population based on conV 
      pop = sortrows(pop, rows);
  end
  
  % TODO 之后可以考虑 rank selection 的问题
  
  r1 = randi(popsize);
  r2 = randi(par_archive_size + popsize);
  
  while r1 == r2 || r1 == i || r2 == i
      r1 = randi(popsize);
      r2 = randi(par_archive_size + popsize);
  end
  
  xr1 = pop(r1, :);
  if r2 > popsize
      xr2 = pop(r2 - popsize, :);
  else
      xr2 = pop(r2, :);
      
  
  % trial vector
  v = xi + F_w * (x_pbetter - xi) + F * (xr1 - xr2);
  
  % crossover 
  cross_rand = rand(1, N);
  cross_rand(randi(N)) = 0;     % at least one decision component to crossover
  cross_ID = cross_rand < CR;   % inherited from v
  v(~cross_ID) = xi(~cross_ID); % inherited from xi
  
  u = v;
  
end

