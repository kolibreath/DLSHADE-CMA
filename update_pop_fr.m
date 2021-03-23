function [pop_struct, archive_fr, archive, suc_f, suc_cr,delta_k] = update_pop_fr(pop_struct,ui,archive,f,cr)
%COMPARE_FR update population using FROFI
% input:
    % pop_struct   -- struct of population matrix, problem_size etc.
    % ui           -- offspring population of pop
    % archive      -- archive stores defeated parents
    % f            -- scale factor used by ui
    % cr           -- crossover rate used by ui
% output:
    % pop_struct   -- struct saves updated population 
    % archive_fr   -- defeated offspring in ui_fr
    % archive      -- updated archive
    % suc_f        -- successful scale factor by ui
    % suc_cr       -- successful crossover rate by ui
    % both suc_f and suc_cr are column vector (suc_num * 1)
    % delta_k      -- improvement of successful offspring compared to
    % parent considering both fitness and conv
    
% Steps: 
%  1) if x1 and x2 are both feasible and the one with better fitness will
%  be selected. 
%  2) if x1 is feasible and x2 is not, x1 will be selected.
%  3) if x1 and x2 are both infeasible, the one with less conV will be
%  seleted
%  4) According to FROFI, ones from ui not selected by feasibilty rule but
%  better in fitness will be store in an archive which is used to replace
%  some individuals of pop.

%  NOTE: If the selected one is from pop, it will go to archive (shared by 
%  subpopulations), otherwise the one from ui will be stored at archive_fr

% Version 1.2 Author: Shi Zeyuan 734780178@qq.com Date: 2021/3/18
    
%% 
    suc_cr = [];
    suc_f = [];
    archive_fr = [];
    
    archive_frofi = [];
    
    delta_k = [];
    popsize = pop_struct.popsize;
    pop = pop_struct.pop;

    % compare between offspring individual and parent
    for k = 1 : popsize
       cur_par = pop(k,:);
       cur_off = ui(k,:);
       
       par_fit = cur_par(end - 1);
       off_fit = cur_off(end - 1);
       
       par_conv = cur_par(end);
       off_conv = cur_off(end);
       
       %% between two infeasible solutions, the one with smaller conV wins
       if par_conv ~=0 && off_conv ~= 0
          if par_conv > off_conv
             %parent is defeated
             [pop,archive,suc_f,suc_cr, delta_k] = replace_record(pop,k,cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k);
          else
              % ui is defeated 
             [archive_fr, archive_frofi] = save_archive(archive_fr, archive_frofi, cur_par, cur_off);
          end
       %% if both feasible, better fitness is preferred
       elseif par_conv == 0 && off_conv == 0
           if par_fit > off_fit 
              %parent is defeated
              [pop,archive,suc_f,suc_cr, delta_k] = replace_record(pop,k,cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k);
           else
               % ui is defeated 
              [archive_fr, archive_frofi] = save_archive(archive_fr, archive_frofi, cur_par, cur_off);
           end
       %% feasible one is better than its infeasible counterpart
       else 
           if off_conv == 0 
             %parent is defeated
             [pop,archive,suc_f,suc_cr, delta_k] = replace_record(pop,k,cur_off,archive,cur_par,suc_f,suc_cr,f,cr,delta_k);
           else 
                % ui is defeated 
              [archive_fr, archive_frofi] = save_archive(archive_fr, archive_frofi, cur_par, cur_off);
           end
       end
    end
   
    %% applying FROFI replacement strategy
    [length,~] = size(archive_frofi);
    if length ~= 0
        pop = replacement(pop_struct, archive_frofi);
    end
    
    pop_struct.pop = pop;  % rewrite updated population into pop_struct
end
    

function [archive_fr, archive_frofi] = save_archive(archive_fr, archive_frofi, cur_par, cur_off)
% SAVE_ARCHIVE when ui is defeated, check if its fitness better than its
% parent, if it is better than parent, save it into archive_frofi, beforing
% applying replacement strategy in FROFI; otherwise save it into
% archive_fr, which will be later used for replacing individuals in pop_ec.
% input: 
    % archive_fr            -- save defeated offspring from pop_fr
    % arhive_frofi          -- save defeated offspring from pop_fr, if its fitness better than its parent
    % cur_par               -- current parent
    % cur_off               -- current offspring
% output    
    % archive_fr            -- save defeated offspring from pop_fr
    % arhive_frofi          -- save defeated offspring from pop_fr, if its fitness better than its parent
    
    if cur_par(end-1) > cur_off(end-1)
        archive_frofi = [archive_frofi; cur_off];
    else
        archive_fr = [archive_fr; cur_off];
    end
end

function  pop = replacement(pop_struct,archive_frofi)
% FRORI replacement strategy
% input:
    % pop_struct      -- population struct
    % archive_frofi   -- defeated individual but have better fitness than parent
% output:
    % pop            -- updated population 

%%
 % calculate the size of the population p(popsize) and the number of dimensions(n) of
 % the tested function
 pop = pop_struct.pop;
 popsize = pop_struct.popsize;
 n = pop_struct.problem_size;
 
 % the maximum number of vectors to be replaced
 N=round(max(5,n/2)); 
 
 % the number of parts to be divided
 MRN=round(popsize/N);
 
 objF = pop(:, end-1);
 conV = pop(:, end);
 
 recordobjF = archive_frofi(:, end-1);
 recordconV = archive_frofi(:, end);
 % sort the objective function value in descendant order
 [~,sortindex]=sort(-objF); 
 
 % divide the population into Nf parts according to their objective function values 
 % in descendant order and execute the replacement operation
 for i=1:floor(popsize/MRN)
  
   % calculate the current number of the recorded vectors  
   len=length(archive_frofi);
   
   % when the recored set is not empty, excuted the replacement operation
   if len~=0 
      
      % calculate index of the vector which has maximum degree of
      % constraint violation in the ith part
      subConV=conV(sortindex((i-1)*MRN+1:(i-1)*MRN+MRN));
      [~,maxIndex]=max(subConV);     
      maxIndex=(i-1)*MRN+maxIndex;
      maxSubConVIndex=sortindex(maxIndex);  
      
      % calculate index of the vector which has minimum degree of
      % constraint violation in the recorded population
      [~,minRecordconVIndex]=min(recordconV);
      
      % replacement according to the objective function value
      if recordobjF(minRecordconVIndex) < objF(maxSubConVIndex)
          
          % replacement  
          pop(maxSubConVIndex,:) = archive_frofi(minRecordconVIndex,:);
          
          % delete the corresponding vector from the recorded population
          recordconV(minRecordconVIndex)=[];
          recordobjF(minRecordconVIndex)=[];
          archive_frofi(minRecordconVIndex,:)=[];     
     end   
  end
 end

end