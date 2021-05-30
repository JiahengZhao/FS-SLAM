function [ alpha ] = LineSearchMy(feaOccurredID, Xstate,Zstate,Ellipse,Line,invCovOdom,invCovEllipse,invCovLine, d)
%LINESEARCH Summary of this function goes here
%   Using line search to determine the exact step for 
%   descending. alpha * d is the actual descending step
%   x:current point
%   d:descending direction
%   F:f(x+0*h)
%   P:covariance matrix

alpha = 0.5;   %initial alpha
step = 0.25;   %increase step
%0 < gamma_1 < gamma_2 < 1
gamma_1 = 0.1;  
gamma_2 = 0.9;

[fval_odom, fval_Ellipse, fval_line] = FuncfAllPoints(feaOccurredID,Xstate,Zstate,Ellipse);
F = [fval_odom;fval_Ellipse;fval_line];
J = FunJacAllPoints(feaOccurredID,Xstate,Ellipse);

% invCov
num1 = length(fval_odom);
num2 = length(fval_Ellipse);
num3 = length(fval_line);
num0 = size(F,1);
if (~isempty(invCovEllipse)) && (isempty(invCovLine)) % only Ellipse
    b = [repmat(invCovOdom,1, num1/3) repmat(invCovEllipse,1, num2)]';
elseif (isempty(invCovEllipse)) && (~isempty(invCovLine)) % only Line
    b = [repmat(invCovOdom,1, num1/3) repmat(invCovLine,1, num3)]';
elseif (~isempty(invCovEllipse)) && (~isempty(invCovLine)) % Ellipse and Line
    b = [repmat(invCovOdom,1, num1/3) repmat(invCovEllipse,1, num2) repmat(invCovLine,1, num3)]';
end
P = spdiags(b,0,num0,num0);

phi_0 = F' * P * F;
phi_0_derivative = d' * J' * P * F;

[Xstate_alpha,Ellipse_alpha] = FunUpdateAllPoints(Xstate,Ellipse,[],alpha*d);
[fval_odom, fval_Ellipse, fval_line] = FuncfAllPoints(feaOccurredID,Xstate_alpha,Zstate,Ellipse_alpha);
F_alpha = [fval_odom;fval_Ellipse;fval_line];
J_alpha = FunJacAllPoints(feaOccurredID,Xstate_alpha,Ellipse_alpha);

phi_alpha = F_alpha' * P * F_alpha;
phi_alpha_derivative = d' * J_alpha' * P * F;

k = 0;
if((phi_alpha>phi_0+gamma_1*phi_0_derivative*alpha)&&(phi_alpha_derivative<gamma_2*phi_0_derivative))
    alpha = 0;
    return;
end
while(~((phi_alpha<=phi_0+gamma_1*phi_0_derivative*alpha)&&(phi_alpha_derivative>=gamma_2*phi_0_derivative))&&k<10)
    if(phi_alpha>phi_0+gamma_1*phi_0_derivative*alpha)
        alpha = alpha - step;
        step = step / 2;
        [Xstate_alpha,Ellipse_alpha] = FunUpdateAllPoints(Xstate,Ellipse,[],alpha*d);
        [fval_odom, fval_Ellipse, fval_line] = FuncfAllPoints(feaOccurredID,Xstate_alpha,Zstate,Ellipse_alpha);
        F_alpha = [fval_odom;fval_Ellipse;fval_line];
        phi_alpha = F_alpha' * P * F_alpha;
        phi_alpha_derivative = d' * J_alpha' * P * F;
    elseif(phi_alpha_derivative<gamma_2*phi_0_derivative)
        alpha = alpha + step;
        step = step / 2;
        [Xstate_alpha,Ellipse_alpha] = FunUpdateAllPoints(Xstate,Ellipse,[],alpha*d);
        [fval_odom, fval_Ellipse, fval_line] = FuncfAllPoints(feaOccurredID,Xstate_alpha,Zstate,Ellipse_alpha);
        F_alpha = [fval_odom;fval_Ellipse;fval_line];
        phi_alpha = F_alpha' * P * F_alpha;
        phi_alpha_derivative = d' * J_alpha' * P * F;
    end
    k = k + 1;
end

end

