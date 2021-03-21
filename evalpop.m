function pop = evalpop(pop,func)
%EVALPOP evaluate the population by test function denoted by 'func', after
% evaluation, individual will be a problem_size + 2 (problem_size , fitness , conV)
% input:
    % pop       -- population evaluated
    % func      -- test function no.
    
% output:   
    % pop       -- population after evaluation
% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

    [f, g, h] = CEC2017(pop, func);   % f: fitness g: inequality constraints h:equality constraints
    conV = overall_cv(g, h);
    pop = [pop,g, h, f, conV];
end

function result = overall_cv(g,h)
    h=abs(h)-1e-4;
    cv=[g,h];
    cv(cv < 0) = 0;
    result = sum(cv,2); 
end




