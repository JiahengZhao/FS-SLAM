function [V, center, newpt] = fitWithFS_update( N, Xi, Yi, center_in, isAdd)
% FITWITHFS fit the feature with fourier series
% Author: Jiaheng Zhao, Tiancheng Li
% Define: click point counter-clockwizely.

if isempty(center_in)
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
% theta_start = theta(1);
% theta_end=theta(end);
r = sqrt(pt_noCenter(2,:).^2 + pt_noCenter(1,:).^2);
% r_end = r(end);
% r_start = r(1);
[theta,od] = sort(theta);
r=r(od);

% find empty area
step = 0;
if isAdd
    %%
    bk = find(abs(diff(theta).*180/pi)>15); % break id. may be more than 1.
    res_deg = trimmean(diff(theta),60);
    theta_1 = theta(bk);
    r_1 = r(bk);
    theta_2 = theta(bk+1); % insert angle between theta_1 and theta_2
    r_2 = r(bk+1);
%     newtheta = theta;
    theta_i = cell(1,length(theta_1));
    r_i = theta_i;
    for index = 1 : length(theta_1)
        if theta_1(index) < theta_2(index)
            step = res_deg;
        else 
            step = -res_deg;
        end
        theta_i{1,index} = wrapTo2Pi(theta_1(index) : step : theta_2(index));
        
%         if r_1 < r_2
            rstep = (r_2 - r_1)/(length(theta_i{1,index})-1);
%         else
%             rstep = -(r_1 - r_2)/length(theta_i{1,index});
%         end
        r_i{1,index} = r_1 : rstep: r_2;
    end
    theta_b = wrapTo2Pi(theta(end) : step : (theta(1)+2*pi));
    rstep =  (r(1) - r(end))/(length(theta_b)-1);

    r_b =  r(end) : rstep: r(1);
    pt_b = [r_b.*cos(theta_b);r_b.*sin(theta_b)];

    r_i = cat(2,r_i{:});
    theta_i = cat(2,theta_i{:});
    pt_Comple=[r_i.*cos(theta_i);r_i.*sin(theta_i)];
    pt_noCenter_sup=[pt_noCenter, pt_Comple, pt_b];
    theta_sup = [theta, theta_i, theta_b];
    [theta_sup,od1] = sort(theta_sup);
    r_sup = [r, r_i, r_b];
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
% for k = 0:N
%     if k ~= 0
%         Ck = [Ck; sum(r_sup.*cos(k.*theta_sup))];
%         Ck2 = [Ck2; sum(r_sup.*sin(k.*theta_sup))];
%     end
%     for n = 0:N
%         P(k+1,n+1) = sum(cos(k.*theta_sup) .*cos(n.*theta_sup));
%         P(k+1,N+n+1) = sum(cos(k.*theta_sup) .* sin(n.*theta_sup));
%         P(N+k+1,n+1) = sum(sin(k.*theta_sup) .* cos(n.*theta_sup));
%         P(N+k+1,N+n+1) = sum(sin(k.*theta_sup) .* sin(n.*theta_sup));
%     end
% end
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