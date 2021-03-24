function sorted_pop = sort_fo(pop)
%SORT_FO sort population according to fitness
% input:
    % pop           -- population will be sorted
% output:
    % sorted_pop    -- population sorted based on fitness
% Steps:
    % 1) individual has better fitness wins

    % Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18

    [~,fitness_sorted_index] = sort(pop(:,end-1),'ascend');
    sorted_pop = [pop(fitness_sorted_index, :)];
end
