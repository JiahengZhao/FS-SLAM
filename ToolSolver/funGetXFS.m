function outX = funGetXFS(Xstate,Zstate)
% This function gets the module of state vector, which is used in LMSBA.
% Author: Jiaheng Zhao

outX = Xstate(:,1)' * Xstate(:,1);
outX = sqrt(outX);

end