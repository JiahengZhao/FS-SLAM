function [V, center, newpt] = fitWithFS( N, Xi, Yi, center_in, isAdd)
% FITWITHFS fit the feature with fourier series
% Author: Jiaheng Zhao, Tiancheng Li
% Define: click point counter-clockwizely.

if isempty(center_in)
    % Fit the feature with ellipse to obtain the center
    %         [fitP,~,~] = fitWithLSO( Xi, Yi, [1;1],0);
    %         center = fitP(1:2,1);
    % Fit the feature with circle to obtain the center.
    [center,~,initialc]= fitWithCircle( Xi, Yi);
    center = initialc;
else
    center = center_in;
end

M = length(Xi);
if size(Xi,1) > size(Xi,2)
    pt = [Xi';Yi'];
else
    pt = [Xi;Yi];
end

% Decentralization
pt_noCenter = pt  - center;
% Decentralized d(theta)
% theta
theta = atan2(pt_noCenter(2,:),pt_noCenter(1,:));
theta = theta .* (theta >= 0) + (theta + 2 * pi) .* (theta < 0);
theta_start = theta(1);
theta_end=theta(end);
r = sqrt(pt_noCenter(2,:).^2 + pt_noCenter(1,:).^2);
r_end = r(end);
r_start = r(1);
[theta,od] = sort(theta);
r=r(od);
% theta_start = theta(1);
% theta_end=theta(end);
% r = sqrt(pt_noCenter(2,:).^2 + pt_noCenter(1,:).^2);
% r=r(od);
% r_end = r(end);
% r_start = r(1);
% find empty area
if isAdd
    %% There are some corner cases.
    if (theta_end - theta_start) > 0
        if max(abs(diff(theta))) > pi/3
            theta_i = wrapTo2Pi(theta_start:deg2rad(5):theta_end);
            if r_end > r_start
                r_i=r_end : -abs((r_end - r_start)/length(theta_i)): r_start;
            elseif r_end < r_start
                r_i=r_end : abs((r_end - r_start)/length(theta_i)) : r_start;
            else
                r_i=linspace(r_end, r_start, length(theta_i)+1);
            end
            r_i = flip(r_i);
            r_i(1) = [];
        else
            theta_i=[wrapTo2Pi(theta_end:deg2rad(5):2*pi),wrapTo2Pi(0:deg2rad(5):theta_start)];
            if r_end > r_start
                r_i=r_end : -abs((r_end - r_start)/length(theta_i)): r_start;
            elseif r_end < r_start
                r_i=r_end : abs((r_end - r_start)/length(theta_i)) : r_start;
            else
                r_i=linspace(r_end, r_start, length(theta_i)+1);
            end
            r_i(1) = [];
        end
    else
        if max(abs(diff(theta))) > pi/3
            theta_i = wrapTo2Pi(theta_end:deg2rad(5):theta_start);
            if r_end > r_start
                r_i=r_end : -abs((r_end - r_start)/length(theta_i)): r_start;
            elseif r_end < r_start
                r_i=r_end : abs((r_end - r_start)/length(theta_i)) : r_start;
            else
                r_i=linspace(r_end, r_start, length(theta_i)+1);
            end
            r_i(1) = [];
        else
            theta_i=[wrapTo2Pi(0:deg2rad(5):theta_end),wrapTo2Pi(theta_start:deg2rad(5):2*pi)];
            sn = length(0:deg2rad(5):theta_end);
            sn1 = length(theta_i);
            if r_end > r_start
                r_i=r_end : -abs((r_end - r_start)/length(theta_i)): r_start;
            elseif r_end < r_start
                r_i=r_end : abs((r_end - r_start)/length(theta_i)) : r_start;
            else
                r_i=linspace(r_end, r_start, length(theta_i)+1);
            end
            r_i = flip(r_i);
            r_i = [r_i(end-sn+1:end), r_i(1:sn1-sn+1)];
            r_i(1) = [];
        end
    end
    
   
    pt_Comple=[r_i.*cos(theta_i);r_i.*sin(theta_i)];
    pt_noCenter_sup=[pt_noCenter,pt_Comple];
    theta_sup = [theta, theta_i];
    [theta_sup,od1] = sort(theta_sup);
    r_sup = [r, r_i];
    r_sup = r_sup(od1);
    M_sup = length(r_sup);
else
    r_sup=r;
    theta_sup = theta;
    M_sup = M;
end

% Equations in the paper.
C0 = sum(r_sup);
Ck = [];
Ck2 = [];
P = zeros(2*N+1);

for k = 0:N
    if k > 0
        Ck = [Ck; sum(r_sup.*cos(k.*theta_sup))];
        Ck2 = [Ck2; sum(r_sup.*sin(k.*theta_sup))];
    end
    for n = 1:N
        if k==0
            P(1,n+1) = sum(cos(n.*theta_sup));
            P(1,N+n+1) = sum(sin(n.*theta_sup));
        else
            P(k+1,n+1) = sum(cos(k.*theta_sup) .*cos(n.*theta_sup));
            P(k+1,N+n+1) = sum(cos(k.*theta_sup) .* sin(n.*theta_sup));
            P(N+k+1,n+1) = sum(sin(k.*theta_sup) .* cos(n.*theta_sup));
            P(N+k+1,N+n+1) = sum(sin(k.*theta_sup) .* sin(n.*theta_sup));
        end
    end
end
P(1,1) = length(theta_sup);
P(2:N+1,1) = P(1,2:N+1)';
P(N+2:end,1) = P(1,N+2:end);

C =[C0; Ck; Ck2];
V = P\C;
if rcond(P) < 1e-14
    V=[];
    center=[];
    newpt=[];
    return;
end
    

an = V(1:N+1,1);
bn = V(N+2:end,1);
bn = [0;bn];
dtheta = zeros(1,M_sup);
% theta = atan2(pt_noCenter(2,:),pt_noCenter(1,:));
for i = 1:N+1
    dtheta = dtheta + an(i).*cos((i-1).*theta_sup) + bn(i).*sin((i-1).*theta_sup);
end
% to point
newpt = [dtheta.*cos(theta_sup); dtheta.*sin(theta_sup)] + center;
end