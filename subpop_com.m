function [pop_fr,pop_ec, delete_individual] = subpop_com(pop_fr,pop_ec,archive_fr,archive_ec)
%SUBPOP_COM communication between subpopulations
% inputs:
    % pop_fr                -- population using feasibility rule
    % pop_ec                -- population using epsilon constraint
    % archive_fr            -- defeated offspring in pop_fr
    % archive_ec            -- defeated offpspring in pop_ec
% outputs:
    % pop_fr                -- population using feasibility rule    (updated)
    % pop_ec                -- population using epsilon constraint  (updated)
    % deleted_individual    -- defeated individuals includes offsprings and parents
%   Steps:
    % 1) find unions of pop_fr and achive_ec, pop_ec and archive_fr
    % 2) sort unions by sort_fr and sort_ec respectively
    % 3) delete overflowing individuals
    
% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18
    
%%
    union_fr = [pop_fr; archive_ec];
    union_ec = [pop_ec; archive_fr];
    
    [len_union_fr, ~] = size(union_fr);
    [len_union_ec, ~] = size(union_ec);
    
    sorted_union_fr = sort_fr(union_fr);
    sorted_union_ec = sort_ec(union_ec);
    
    % delete 
    pop_fr = sorted_union_fr(1:pop_fr.popsize,:);
    pop_ec = sorted_union_ec(1:pop_ec.popsize,:);
    
    % deleted 
    delete_individual = [sorted_union_fr(pop_fr.popsize+1:len_union_fr,:)];
    delete_individual = [delete_individual;sorted_union_ec(pop_ec.popsize+1:len_union_ec,:)];
end

