function [center,r,rc]= fitWithCircle( Xi, Yi)
% FITWITHCIRCLE fit the feature with circle.
% Author: Jiaheng Zhao.
% center = [c1; c2], radius = r
n= length(Xi);
% initialc = [sum(Xi)/n; sum(Yi)/n];

minx = min(Xi);
maxx = max(Xi);
centx = (minx + maxx) / 2;
miny = min(Yi);
maxy = max(Yi);
centy = (miny + maxy) / 2;
dist2 = (Xi - centx).^2 + (Yi - centy).^2;
rc = [centx;centy];
[~, idx] = min(dist2);
bestx = Xi(idx);
besty = Yi(idx);
initialc = [bestx; besty];

initialr = sqrt(sum((Xi - initialc(1)).^2 + (Yi - initialc(2)).^2)/n);

% F = @(x) sum(((Xi-x(1)).^2 + (Yi-x(2)).^2 - x(3).^2).^2);
F = @(x) sum(abs((Xi-x(1)).^2 + (Yi-x(2)).^2 - x(3).^2));

options = optimoptions(@fminunc,'Display','none');
x = fminunc(F,[initialc;initialr],options);
center = x(1:2,1);
r = x(3);

end