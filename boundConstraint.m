function vi = boundConstraint (vi,pop,r0,lu)

% if the boundary constraint is violated, set the value to be the middle
% of the previous value and the bound
%
% Version: 1.1   Date: 11/20/2007
% Written by Jingqiao Zhang, jingqiao@gmail.com
% Modified by Shi Zeyuan, 734780178@qq.com Date: 03/24/2021

    [~, NP] = size(r0);  % the population size and the problem's dimension
    
    % if violated set a random value in the lower and upper bounds
    %% check the lower bound
    xl = repmat(lu(1, :), NP, 1);
    pos = vi < xl;
    pop = pop(r0,:);
    vi(pos) = (pop(pos) + xl(pos)) / 2;
    
    %% check the upper bound
    xu = repmat(lu(2, :), NP, 1);
    pos = vi > xu;
    pop = pop(r0,:);
    vi(pos) = (pop(pos) + xu(pos)) / 2;

end