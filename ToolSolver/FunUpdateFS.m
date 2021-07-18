function Xstate = FunUpdateFS(Xstate, Delta)
% This fucntion update the delt in 3 cases.
% Author: Jiaheng Zhao

Xstate(:,1) = Xstate(:,1) + Delta;

Xstate(Xstate(:,2)==1,1) = WrapX(Xstate(Xstate(:,2)==1,1),Xstate(Xstate(:,2)==1,:)); %wrap

end