function [fval, info] = intersection(F,ang,x0,y0)
phiF = F(3);
cF = cos(phiF);
sF = sin(phiF);
xF = F(1);
yF = F(2);
r1 = F(4);
r2 = F(5);
if ang ~= pi/2 && ang~= -pi/2
    k = tan(ang);
    m1 = cF+k*sF;
    m2 = sF-k*cF;
    n1 = (y0-yF)*sF - xF*cF - k*x0*sF;
    n2 = (y0-yF)*cF + xF*sF - k*x0*cF;
    a = m1^2*r2^2 + m2^2*r1^2;
    b = 2*m1*n1*r2^2 - 2*m2*n2*r1^2;
    c = n1^2*r2^2 + n2^2*r1^2 - r1^2*r2^2;
    delt = b^2-4*a*c;
    if delt < 0
        fval = [];
        info = 0; % 0 means no intersection
    else
        %         info = 1;
        x1 = (-b + sqrt(delt))/(2*a) ;
        x2 = (-b - sqrt(delt))/(2*a) ;
        y1 = y0+k*(x1-x0);
        y2 = y0+k*(x2-x0);
        dist1 = (y1-y0)^2+(x1-x0)^2;
        temp1 =  atan2(y1-y0,x1-x0); % because x1y1 and x2y2 lie on the same side of the beam
        if abs(temp1 - ang) < 1e-15
            info = 1;
        else
            info = 0;
            fval = [];
            return
        end
        dist2 = (y2-y0)^2+(x2-x0)^2;
        if dist1 < dist2
            fval = [x1;y1];
        else
            fval = [x2;y2];
        end
    end
else
    p1=(x0-xF)*cF;
    p2=(x0-xF)*sF;
    q1=yF*sF;
    q2=yF*cF;
    a=r2^2*sF^2+r1^2*cF^2;
    b=2*(p1-q1)*sF*r2^2-2*(p2-q2)*cF*r1^2;
    c=(p1-q1)^2*r2^2+(p2-q2)*r1^2-r1^2*r2^2;
    delt = b^2-4*a*c;
    if delt < 0
        fval = [];
        info = 0; % 0 means no intersection
    else
        info = 1;
        y1 = (-b + sqrt(delt))/(2*a);
        y2 = (-b - sqrt(delt))/(2*a);
        x1 = x0;
        if ang == pi/2
            fval = [x1;y2];
        elseif ang == -pi/2
            fval = [x1;y1];
        end
    end
end
end