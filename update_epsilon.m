function epsilon = update_epsilon(epsilon_zero,nfes,fes_control,cp)
%UPDATE_EPSILON update epsilon in epsilon-constraint handling technique
% input:
    % epsilon_zero      -- the theta-th conv of individual of initial population ranked based on conV
    % nfes              -- current fes
    % fes               -- control the value of epsilon
    % cp                -- parameter related to epsilon control
% output:
    % epsilon           -- updated epsilon
% Step
 % epsilon = epsilon_zero * (1 - nfes / fes_control)^ cp  when 0 < nfes < fes_control
 % epsilon = 0     otherwise
 % 1) 
% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/20
    if nfes < fes_control
        epsilon = epsilon_zero * (1 - nfes / fes_control) ^ cp;
    else
        epsilon = 0;
    end
end

