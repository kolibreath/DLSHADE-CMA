function  pop = replacement(pop,archive)
% FRORI replacement strategy
% input:
    % pop             -- population
    % archive_frofi   -- defeated individual but have better fitness than parent
% output:
    % p               -- updated population

%%
 % calculate the size of the population p(popsize) and the number of dimensions(n) of
 % the tested function
 popsize = pop.popsize;
 n = pop.problem_size;
 
 % the maximum number of vectors to be replaced
 N=round(max(5,n/2)); 
 
 % the number of parts to be divided
 MRN=round(popsize/N);
 
 objF = pop(:, end-1);
 conV = pop(:, end);
 
 recordobjF = archive(:, end-1);
 recordconV = archive(:, end);
 % sort the objective function value in descendant order
 [~,sortindex]=sort(-objF); 
 
 % divide the population into Nf parts according to their objective function values 
 % in descendant order and execute the replacement operation
 for i=1:floor(popsize/MRN)
  
   % calculate the current number of the recorded vectors  
   len=length(archive);
   
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
          pop(maxSubConVIndex,:) = archive(minRecordconVIndex,:);
          
          % delete the corresponding vector from the recorded population
          recordconV(minRecordconVIndex)=[];
          recordobjF(minRecordconVIndex)=[];
          archive(minRecordconVIndex,:)=[];     
     end   
  end
 end

end