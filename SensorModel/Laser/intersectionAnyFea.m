function [fval, info] = intersectionAnyFea(F,ang,x0,y0)

v = num2cell(F);
% [a, b, c, d, e, f, g, h, ii, jj, k, l, m] = deal(v{:}); % assume the feature is fixed.
if ang ~= pi/2 && ang~= -pi/2
    n = tan(ang);
    
    pa0 = 90*n^4 - 4*n + 10;
    pa1 = 4*y0 - 24*n - 4*n*x0 + 360*n^4*x0 - 360*n^3*y0 + 10*n^2 + 720*n^3 + 72;
    pa2 = 540*n^4*x0^2 - 1080*n^3*x0*y0 + 2160*n^3*x0 + 20*n^2*x0 + 540*n^2*y0^2 - 2160*n^2*y0 + 2180*n^2 - 24*n*x0 - 20*n*y0 - 10*n + 24*y0 + 192;
    pa3 = 90*n^4*x0^4 - 360*n^3*x0^3*y0 + 720*n^3*x0^3 + 540*n^2*x0^2*y0^2 - 2160*n^2*x0^2*y0 + 2180*n^2*x0^2 - 360*n*x0*y0^3 + 2160*n*x0*y0^2 - 4360*n*x0*y0 + 2924*n*x0 + 90*y0^4 - 720*y0^3 + 2180*y0^2 - 2924*y0 + 32159/20;
    pa4 = 360*n^4*x0^3 - 1080*n^3*x0^2*y0 + 2160*n^3*x0^2 + 10*n^2*x0^2 + 1080*n^2*x0*y0^2 - 4320*n^2*x0*y0 + 4360*n^2*x0 - 20*n*x0*y0 - 10*n*x0 - 360*n*y0^3 + 2160*n*y0^2 - 4360*n*y0 + 2924*n + 10*y0^2 + 10*y0 + 260;
    
    Pa0 = (3*pa1^2)/(8*pa0^2) - pa2/pa0;
    Pa1 = pa3/pa0 - (3*pa1^4)/(256*pa0^4) + (pa1^2*pa2)/(16*pa0^3) - (pa1*pa4)/(4*pa0^2);
    Pa2 = pa1^3/(8*pa0^3) + pa4/pa0 - (pa1*pa2)/(2*pa0^2);
    
    PA0 = - 256*Pa1^3 - 4*Pa0^3*Pa2^2 - 16*Pa0^4*Pa1 + 128*Pa0^2*Pa1^2 + 27*Pa2^4 + 144*Pa0*Pa2^2*Pa1;
    PA1 = 128*Pa0^2*Pa1^2 - 4*Pa0^3*Pa2^2 - 16*Pa0^4*Pa1 - 256*Pa1^3 + 27*Pa2^4 + 144*Pa0*Pa2^2*Pa1;
    PA2 = - 2*Pa0^3 + 3*3^(1/2)*PA0^(1/2) + 27*Pa2^2 + 72*Pa0*Pa1;
    
    PAA0 = - Pa0^3/27 + (3^(1/2)*PA0^(1/2))/18 + Pa2^2/2 + (4*Pa0*Pa1)/3;
    PAA1 = Pa0^2 + (12*pa3)/pa0 - (9*pa1^4)/(64*pa0^4) + 9*PAA0^(2/3) + 6*Pa0*PAA0^(1/3) + (3*pa1^2*pa2)/(4*pa0^3) - (3*pa1*pa4)/pa0^2;
    PAA2 = 6*((3^(1/2)*PA1^(1/2))/18 - Pa0^3/27 + Pa2^2/2 + (4*Pa0*Pa1)/3)^(1/6);
    PAA3 = - 12*Pa1*PAA1^(1/2) - 9*PAA0^(2/3)*PAA1^(1/2) - Pa0^2*PAA1^(1/2) + 12*Pa0*PAA0^(1/3)*PAA1^(1/2);
    PAA4 = 6*((3^(1/2)*PA1^(1/2))/18 - Pa0^3/27 + Pa2^2/2 + (4*Pa0*Pa1)/3)^(1/6)*(Pa0^2 + (12*pa3)/pa0 - (9*pa1^4)/(64*pa0^4) + 9*((3^(1/2)*PA1^(1/2))/18 - Pa0^3/27 + Pa2^2/2 + (4*Pa0*Pa1)/3)^(2/3) + 6*Pa0*((3^(1/2)*PA1^(1/2))/18 - Pa0^3/27 + Pa2^2/2 + (4*Pa0*Pa1)/3)^(1/3) + (3*pa1^2*pa2)/(4*pa0^3) - (3*pa1*pa4)/pa0^2)^(1/4);
    
    x1 = pa1/(4*pa0) - PAA1^(1/2)/PAA2 - (PAA3 - 3*6^(1/2)*Pa2*PA2^(1/2))^(1/2)/PAA4;
    x2 = pa1/(4*pa0) - PAA1^(1/2)/PAA2 + (PAA3 - 3*6^(1/2)*Pa2*PA2^(1/2))^(1/2)/PAA4;
    x3 = pa1/(4*pa0) + PAA1^(1/2)/PAA2 - (PAA3 + 3*6^(1/2)*Pa2*PA2^(1/2))^(1/2)/PAA4;
    x4 = pa1/(4*pa0) + PAA1^(1/2)/PAA2 + (PAA3 + 3*6^(1/2)*Pa2*PA2^(1/2))^(1/2)/PAA4;
    
    realX = [isreal(x1) isreal(x2) isreal(x3) isreal(x4)];
    if ~any(realX)%(~isreal(x1)) && (~isreal(x2)) && (~isreal(x3)) && (~isreal(x4))
        fval = [];
        info = 0; % 0 means no intersection
        return
    end
    X = [x1 x2 x3 x4]';
    X = X(realX);
    Y = y0 + n.*(X - x0);
    dist = (Y-y0).^2 + (X-x0).^2;
    temp1 = atan2(Y-y0,X-x0);
    [~,id]=min(dist);
    if abs(temp1(id) - ang) < 1e-15 % <
        info = 1;
        fval = [X(id); Y(id)];
    else
        info = 0;
        fval = [];
        return
    end
else
    
    pb0 = (2*x0^3)/45 - (4*x0^2)/15 + (5*x0)/9 - 2/5;
    pb1 = x0^4/9 - (8*x0^3)/9 + (8*x0^2)/3 - (32*x0)/9 + 3199/1800;
    
    Pb0 = 27*pb0^4 + 128*(x0/9 - 2/9)^2*pb1^2 - 256*pb1^3 - 16*(x0/9 - 2/9)^4*pb1 - 4*(x0/9 - 2/9)^3*pb0^2 + 144*(x0/9 - 2/9)*pb0^2*pb1;
    Pb1 = pb0^2/2 + (4*(x0/9 - 2/9)*pb1)/3 - (x0/9 - 2/9)^3/27 + (3^(1/2)*Pb0^(1/2))/18;
    Pb2 = 27*pb0^2 + 72*(x0/9 - 2/9)*pb1 - 2*(x0/9 - 2/9)^3 + 3*3^(1/2)*Pb0^(1/2);
    Pb3 = 3*6^(1/2)*Pb2^(1/2)*pb0;
    
    PB0 = - (128*x0)/3 + 6*(x0/9 - 2/9)*Pb1^(1/3) + 9*Pb1^(2/3) + (x0/9 - 2/9)^2 + 32*x0^2 - (32*x0^3)/3 + (4*x0^4)/3 + 3199/150;
    PB1 =    6*(x0/9 - 2/9)*Pb1^(1/3) - (128*x0)/3 + 9*Pb1^(2/3) + (x0/9 - 2/9)^2 + 32*x0^2 - (32*x0^3)/3 + (4*x0^4)/3 + 3199/150;
    if PB0 == PB1
        disp('PB0=PB1')
    end
    PB2 = - 9*Pb1^(2/3)*PB0^(1/2) - (x0/9 - 2/9)^2*PB0^(1/2) - 12*pb1*PB0^(1/2) + 12*(x0/9 - 2/9)*Pb1^(1/3)*PB0^(1/2);
    PB4 = 6*Pb1^(1/6);
    PB3 = PB4*PB1^(1/4);
    
    y1 = - PB0^(1/2)/PB4 - (PB2 - Pb3)^(1/2)/PB3 + 2;
    y2 = - PB0^(1/2)/PB4 + (PB2 - Pb3)^(1/2)/PB3 + 2;
    y3 =   PB0^(1/2)/PB4 - (PB2 + Pb3)^(1/2)/PB3 + 2;
    y4 =   PB0^(1/2)/PB4 + (PB2 + Pb3)^(1/2)/PB3 + 2;
    
    realY = [isreal(y1) isreal(y2) isreal(y3) isreal(y4)];
    if ~any(realY)%(~isreal(x1)) && (~isreal(x2)) && (~isreal(x3)) && (~isreal(x4))
        fval = [];
        info = 0; % 0 means no intersection
        return
    end
    Y = [y1 y2 y3 y4]';
    Y = Y(realY);
    dist = (Y-y0).^2;
    [~,id]=min(dist);
    fval = [x0;Y(id)];
    info = 1;
    
    
end
end