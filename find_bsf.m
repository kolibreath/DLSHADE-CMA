function bsf_solution = find_bsf(pop_struct, bsf_solution)
%FIND_BSF 此处显示有关此函数的摘要
%   此处显示详细说明
   k = 1;
   while k < pop_struct.popsize
      off = pop_struct.pop(k,:);
      % better in conv
      if off(end) < bsf_solution(end)
         bsf_solution = off;
         % equal conv better in fitness
      elseif off(end) == bsf_solution(end) && off(end-1) < bsf_solution(end-1)
         bsf_solution = off;
      end
                  
      k =k +1;
   end
 
end
       


