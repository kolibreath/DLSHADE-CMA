function [pop_fr_struct,pop_fo_struct, delete_individual] = subpop_com(pop_fr_struct,pop_fo_struct,archive_fr,archive_fo)
%SUBPOP_COM communication between subpopulations
% inputs:
    % pop_fr_struct                -- population struct using feasibility rule
    % pop_fo_struct                -- population struct using no constraint handling technique 
    % archive_fr                   -- defeated offspring in pop_fr
    % archive_fo                   -- defeated offpspring in pop_fo
% outputs:
    % pop_fr_struct                -- population using feasibility rule    (updated)
    % pop_fo_struct                -- population using no constraint handling  (updated)
    % deleted_individual           -- defeated individuals includes offsprings and parents
%   Steps:
    % 1) find unions of pop_fr and achive_ec, pop_fo and archive_fr
    % 2) sort unions by sort_fr and sort_ec respectively
    % 3) delete overflowing individuals
    
% Version 1.4 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18
    
%%
    pop_fr = pop_fr_struct.pop;
    pop_fo = pop_fo_struct.pop;
    
    union_fr = [pop_fr; archive_fo];
    union_fo = [pop_fo; archive_fr];
    
    sorted_union_fr = sort_fr(union_fr);
    sorted_union_fo = sort_fo(union_fo);
    
    % delete 
    pop_fr = sorted_union_fr(1:pop_fr_struct.popsize,:);
    pop_fo = sorted_union_fo(1:pop_fo_struct.popsize,:);
    
    % deleted 
    delete_individual = sorted_union_fr(pop_fr_struct.popsize+1:end,:);
    delete_individual = [delete_individual;sorted_union_fo(pop_fo_struct.popsize+1:end,:)];
    
    %TODO 这里两个种群是否交换了？ 是否删除了不要的解
    pop_fr_struct.pop = pop_fr;
    pop_fo_struct.pop = pop_fo;
end

