function pop = evalpop(pop,func)
%EVALPOP evaluate the population by test function denoted by 'func', after
% evaluation, individual will be a problem_size + 2 (problem_size , fitness , conV)
% input:
    % pop       -- population evaluated
    % func      -- test function no.
    
% output:   
    % pop       -- population after evaluation
    [f, g, h] = CEC2017(pop, func);   % f: fitness g: inequality constraints h:equality constraints
    conV = overall_cv(g, h);
    pop = [pop_fr,g, h, f, conV];
end
