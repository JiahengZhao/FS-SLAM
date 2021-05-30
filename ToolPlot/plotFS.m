function p = plotFS(center,V,color)
% PLOTFS plot fouries series using center.
% center = [center1 center2 ... centern] 2 by n;
% V = [v1 v2 ... vn] 2N+1 by n
if nargin < 3
    color = 'r';
end
nF = size(center,2);
N = (size(V,1)-1)/2;
theta = -pi:deg2rad(1):pi;
dtheta = zeros(nF,length(theta));

for nid = 1:nF
    an = V(1:N+1,nid);
    bn = V(N+2:end,nid); bn = [0;bn];
    for i = 1:N+1
        dtheta(nid,:) = dtheta(nid,:) + an(i) .* cos((i-1).*theta) + bn(i) .* sin((i-1).*theta);
    end
    pts = [dtheta(nid,:).*cos(theta); dtheta(nid,:) .* sin(theta)];
    pts = pts + center(:,nid);
    p(nid) = plot(pts(1,:),pts(2,:),'Color',color,'LineWidth',4);
end
end