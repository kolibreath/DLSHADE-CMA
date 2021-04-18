function bsf_solution = find_bsf(pop_fr_struct,pop_fo_struct,bsf_solution)
%FIND_BSF 此处显示有关此函数的摘要
%   此处显示详细说明
 k = 1;
 while k < pop_fr_struct.popsize && k < pop_fo_struct.popsize
    off_fr = pop_fr_struct.pop(k,:);
    off_fo = pop_fo_struct.pop(k,:);
    % better in conv
    if off_fr(end) < bsf_solution(end)
       bsf_solution = off_fr;
       % equal conv better in fitness
       elseif off_fr(end) == bsf_solution(end) && off_fr(end-1) < bsf_solution(end-1)
         bsf_solution = off_fr;
    end
                
    % better in conv
    if off_fo(end) < bsf_solution(end)
       bsf_solution = off_fo;
       % equal conv better in fitness
    elseif off_fo(end) == bsf_solution(end) && off_fo(end-1) < bsf_solution(end-1)
       bsf_solution = off_fo;
    end
    k =k +1;
 end
 
end
       


