% nonlinear discrete time dynamic system motion model
function [xstate1] = Motionmodel(xstate0,input,noise)
if isempty(noise)
    noise = [0 0 0];
end
dx = input(1)+noise(1);
dy = input(2)+noise(2);
dphi = input(3)+noise(3);
R = theta2R(xstate0(3));
xstate1(1:2)=xstate0(1:2)+R*[dx;dy];
xstate1(3)=xstate0(3)+dphi;
xstate1(3)=wrapToPi(xstate1(3)); % limit to [-pi pi]
xstate1 = xstate1';
end
