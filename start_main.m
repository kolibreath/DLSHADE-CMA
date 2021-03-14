%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        author: Zeyuan Shi
%        University: CCNU
%        email: 734780178@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;
%%  PARAMETERS SETTING
Dimension_set=[10,30,50,100];
xmean = rand(N,1); % mean value of population !! 这里的均值的选取和上下界有关系吗？？？
sigma = 0.5;       % step size of population
N = Dimension_set(1); % current dimension set
ini_popsize = N * 5;
min_popsize = 5;
popsize = ini_popsize; %current popsize
I_fno = 1;  % number of benchmark function 
max_fes = 10000;

mentry_size  = 5;   % size of memory 
M_F = ones(mentry_size) * 0.5;
M_CR = ones(mentry_size) * 0.5;

k = 1;              % the start point of memory update of LSHADE memory

%% init population
B = eye(N , N);
D = ones(N, 1);
C = B * diag(D.^2) * B';
invsqrtC = B * diag(D.^ - 1) * B'; % C^-1/2
chiN = N^0.5 * (1 - 1 / (4 * N) + 1 / (21 * N^2)); 

par_archive = [];   % external archive of parent individuals

% pop = NP * (value of each dimension, fitness, conV )
pop = zeros(popsize, N + 2);
for k = 1 : popsize  % 如果按照CMA-ES的生成方式来生成点，这样的点不是均匀分布在空间中？
    pop(:, k) = xmean + sigma * B * (D .* randn(N, 1));
end

% evaluate popluation
[f,g,h] = CEC2017(pop, I_fno); % f = fitness, g = inequality constraints, h = equality constraint
cv = overall_cv(g, h);
cur_fes = popsize;

% divide popluation into 2 subpopluations
half_index = floor(popsize / 2);

pop_FR = pop(half_index, :);                      % feasibility rule
pop_EC = pop(half_index + 1 : popsize, :);        % epsilon-constraint handling

% MAIN PART

%% 
while cur_fes <= max_fes % TODO 还有CMA-ES的终止条件
    SCR = [];
    SF = [];
    
    [popsize_fr, ~] = size(pop_FR);
    [popsize_EC, ~] = size(pop_EC);
    
    for i = 1 : popsize_fr  % iterate subpopulation using feasibility rule
        r_i = randi(mentry_size);
        
    end
    
    for i = 1 : popsize_EC % iterate subpopulation using epsilon constraint handling
    end
end
%%    









