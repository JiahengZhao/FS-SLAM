function [p_up, cov_up] = fitWithEKF(x, y,init)

if nargin < 3
    init = [];
    p = ones(5,1);
end
M = [];
M_x = [];
f=[];
[xcent,ycent]=pol2cart(init(2),init(1));
dsize = init(3)/2;
% p = (dsize^2*xcent^2+idsize^2*ycent^2-dsize^4).*[dsize^2 0 dsize^2 -xcent*dsize -ycent*dsize]';
denominater = (xcent^2+ ycent^2)/(dsize^2) -1;
p=[1./(dsize^2*denominater); 0; 1./(dsize^2*denominater); -xcent/(dsize^2*denominater); -ycent/(dsize^2*denominater)];

% predict
p_pre = p;
cov_pre = eye(5);


% update

for i = 1 : length(x)
    %-------subject to A+C=1
%     M = [M;x(i)^2-y(i)^2 2*x(i)*y(i) 2*x(i) 2*y(i) 1]; % n by 5
%     M_x = [M_x; 2*(x(i)*p_pre(1)+y(i)*p_pre(2)+p_pre(3)) 2*(-y(i)*p_pre(1)+x(i)*p_pre(2)+p_pre(4)+y(i))];
%     f = [f; p_pre(1)*(x(i)^2-y(i)^2) + 2*x(i)*y(i) * p_pre(2) + 2*x(i)*p_pre(3) + 2*y(i)*p_pre(4) + p_pre(5) + y(i)^2]; % n by 1
    %---- subject to F = 1
    M = [M;x(i)^2 2*x(i)*y(i) y(i)^2  2*x(i) 2*y(i)]; % n by 5
    f = [f; p_pre(1)*x(i)^2 + 2*x(i)*y(i) * p_pre(2) + p_pre(3)*y(i)^2 + 2*x(i)*p_pre(4) + 2*y(i)*p_pre(5) + 1]; % n by 1
end
Mcov_pre = (M * cov_pre)'; % 5 by n
s = M * Mcov_pre; % n by n
K = Mcov_pre * pinv(s); % 5 by n
p_up = p_pre - K * f;
cov_up = (eye(5) - K*M)*cov_pre;
cov_up(3,:)=100*cov_up(3,:);
cov_up(:,3)=100*cov_up(:,3);
end

