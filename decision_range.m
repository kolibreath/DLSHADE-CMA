
function range = decision_range(I_fno,D)

range = ones(D,2);      %行代表变量数，列代表上下限

switch I_fno
    case {1,2,3,8,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25,26,27}
        range(:,1)      = -100;
        range(:,2)      = 100;
    case {4,5,9}
        range(:,1)      = -10;
        range(:,2)      = 10;
    case {6}
        range(:,1)      = -20;
        range(:,2)      = 20;  
    case {7,19,28}
        range(:,1)      = -50;
        range(:,2)      = 50;         
end
end