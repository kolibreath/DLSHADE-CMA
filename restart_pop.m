function pop_struct = restart_pop(pop_struct,func)
%RESTART_POP 此处显示有关此函数的摘要
%   此处显示详细说明
       % pick of the best individuals as a initial point
       problem_size = pop_struct.problem_size;
       popsize = pop_struct.popsize;
       
       initial_index = randi(floor(0.3 * popsize));
       xmean = pop_struct.pop(initial_index,1:problem_size);
       
       popsize = popsize * 2; 
       
       pop = zeros(popsize, problem_size);
       sigma = 0.3;
       
       B = eye(problem_size, problem_size); 
       D = ones(problem_size, 1); 
       C = B * diag(D.^2) * B';
       invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
       
       eigeneval = 0; % track update of B and D
       
       k = 1;
       while k <= popsize
           pop(k,:) =  (xmean' + sigma * B * (D .* randn(problem_size, 1)))';
           k = k + 1;
       end     
       pop = evalpop(pop,func);
     
       lambda = popsize * 2; 
      
       pop_struct = assem_pop(pop,popsize,lambda, ... 
       problem_size,C,D,B,invsqrtC,eigeneval,xmean,sigma);
end

