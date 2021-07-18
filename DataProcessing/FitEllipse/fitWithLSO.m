  function [fitP,invCovFea,isill] = fitWithLSO( Xi, Yi, CovScan,Theta_odom)
%FITWITHLSO Ellipse fitting with Least-Squares based on orthogonal distance
%
% [1] Sung Joon Ahn, W. Rauh, and M. Recknagel, "Ellipse fitting and
% parameter assessment of circular object targets for robot vision",
% Intelligent Robots and Systems, 1999.
%
% AUTHOR Sebastian Dingler <s.dingler@gmail.com>
%        Karlsruhe Institute of Technology (KIT), Germany
%
% DATE   22.12.2014

% Calculate covariance of feature
% Upgraded by Jiaheng Zhao
% 2019.03.22

% initial guess
p = fitWithLSA(Xi,Yi);
[Xc,Yc,a,b,alpha] = Matrix2AngleForm(p(1),p(2),1-p(1),p(3),p(4),p(5));
fitP=[Xc,Yc,a,b,alpha]';
if ~isreal(fitP)
    invCovFea=[];
    isill = 1;
    return
end
% ratio = a/b;
% if ratio <= 1.1


% step size
lambda=0.1;
isill = 0;
for k=1:50
    J = [];
    X_new2 = [];
    for i=1:length(Xi)
        r = XY2xy([Xi(i);Yi(i)],alpha,[Xc;Yc]);
        x_new = getOrthoPoint(r(1),r(2),a,b); % orthogonal point in local frame
        X_new = xy2XY2(x_new, alpha, [Xc;Yc]);  % orthogonal point in global frame
        xi = XY2xy([Xi(i);Yi(i)],alpha,[Xc;Yc]);
        J = [ J;calcJacobianMatrix(a,b,x_new(1),x_new(2),alpha,xi(1),xi(2))];
        X_new2 = [X_new2; [Xi(i);Yi(i)]-X_new]; % L1 norm of input points and orthogonal points(in global)
    end
    r=-pinv(J) * X_new2;
    % update
    Xc=Xc-lambda*r(1);
    Yc=Yc-lambda*r(2);
    a=a-lambda*r(3);
    b=b-lambda*r(4);
    alpha=alpha-lambda*r(5);
    
    ratio = a/b;
    
    if k == 50
        CovScan_Full = diag(repmat(CovScan,length(Xi),1)); % CovScan is 2 by 1 vector. CovScan_Full is 2*n by 2*n covariance matrix
        InV_CovScan_Full = CovScan_Full \ eye(size(CovScan_Full,1));
        rcon = rcond(J' * InV_CovScan_Full * J);
        
        if rcon<1e-10
%             disp(['                Illed conditioned, rcond is ', num2str(rcon),'.']);
            isill=2;
        elseif ratio>3  % assume we don't have weird ellipse.
%             disp(['                Fr1/Fr2 is to large: ', num2str(ratio),'.']);
            isill=2;
        end
        %         CovFea_tmp = (J' * InV_CovScan_Full * J) \ eye(5);
        invCovFea_tmp = (J' * InV_CovScan_Full * J);
        %         if abs(ratio-1) <= 0.1
        %            CovFea_tmp = pinv(J' * InV_CovScan_Full * J);
        %            CovFea_tmp(:,5) = CovFea_tmp(:,5)*10000;
        %            CovFea_tmp(5,:) = CovFea_tmp(5,:)*10000;
        %         end
        
        %         CovFea_tmp = pinv(J' * InV_CovScan_Full * J);
        
    end
end
alpha = wrapToPi(alpha); % Wrap alpha ---- something wrong may occur.

% limit feature's angle in 0 to pi.
if ~isempty(Theta_odom)
    FeaAng_global = Theta_odom + alpha;
    if ~isreal(FeaAng_global)
        FeaAng_global
    end
    FeaAng_global = wrapToPi(FeaAng_global);
    % if feature angle is 0 or pi or -pi, normalize it to 0.
    if abs(FeaAng_global) < 1e-10 || abs(abs(FeaAng_global./pi)-1) < 1e-10
        FeaAng_global = 0;
        alpha = FeaAng_global - Theta_odom;
    end
    if FeaAng_global < 0
%         alpha = alpha + pi;
    elseif FeaAng_global > pi
%         alpha = alpha - pi;
    end
end

invCovFea=orderCov(invCovFea_tmp);
fitP = [Xc,Yc,alpha,a,b]';
end

% Jacobian matrix at the orthogonal contacting point on ellipse
function r = calcJacobianMatrix(a, b, x, y, alpha, xi, yi)
C = cos(alpha);
S = sin(alpha);
R = [C S;-S C];
B1 = [b^2 * x * C-a^2 * y * S; b^2 * (yi-y) * C+a^2 * (xi-x) * S];
B2 = [b^2 * x * S+a^2 * y * C;b^2 * (yi-y) * S-a^2 * (xi-x) * C];
B3 = [a * (b^2-y^2); 2 * a * y * (xi-x)];
B4 = [b * (a^2-x^2);-2 * b * x * (yi-y)];
B5 = [(a^2-b^2) * x * y; (a^2-b^2) * (x^2 - y^2 - x * xi + y * yi)];
B = [B1 B2 B3 B4 B5];
Qk = [b^2 * x a^2*y; (a^2-b^2)*y+b^2*yi (a^2-b^2)*x-a^2*xi];
% r = pinv(R)*pinv(Qk)*B;
% r = R'*pinv(Qk)*B;
r = R'*Qk\B;
end

% Covariance matrix in order: x y phi r1 r2
function cov = orderCov(Cov)
cov = Cov;
cov(:,[3 4 5])= cov(:,[5 3 4]);
cov([3 4 5],:) = cov([5 3 4],:);
end
