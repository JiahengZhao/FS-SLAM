function [ptx, pty] = beamFun(x0,y0,ang,range)
% y = y0+k*(x-x0)
if ang ~= pi/2 && ang~= -pi/2
    k = tan(ang);   % tangent of beam function
    sqt=sqrt((range^2)/(1+k^2)); 
    ptx_tmp = [x0+sqt;  x0-sqt ];
    pty_tmp = [y0+k*sqt; y0-k*sqt];
    temp =  atan2(pty_tmp-y0,ptx_tmp-x0);
    
    ptx = ptx_tmp(abs(temp-ang)<1e-15);
    pty = pty_tmp(abs(temp-ang)<1e-15);
elseif ang ==  -pi/2
    ptx=x0;
    pty = y0-range;
elseif ang ==  pi/2
    ptx=x0;
    pty = y0+range;
end
end