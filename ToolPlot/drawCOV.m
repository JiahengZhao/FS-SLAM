function h = drawCOV(mean,cov,color,sdwidth)
if ~exist('sdwidth', 'var');sdwidth = 1; end
if ~exist('color','var');color = 'y';end
if ~exist('npts', 'var'), npts = []; end
if ~exist('axh', 'var'), axh = gca; end

CV = GetCov(cov,mean(1),mean(2),sdwidth);

switch numel(mean)
    case 2
%         plot(CV(1,:),CV(2,:), '-','Color',color); % draw the ellispse
        h{1}=plot(mean(1),mean(2), ':','Color',color); % draw the estimated value 
%         h{2}=fill(CV(1,:), CV(2,:),color,'FaceAlpha',0.1,'LineStyle','none');
        h{2}=fill(CV(1,:), CV(2,:),color,'FaceAlpha',0,'EdgeColor',color,'LineStyle','-.','LineWidth',0.2);
%        h= show2d(mean,cov,color,sdwidth,npts,axh);
        
    case 3
       h= show3d(mean,cov,color,sdwidth,npts,axh);

end


%-----------------------------
function h = show2d(means, C,color, sdwidth, npts, axh)
if isempty(npts), npts=50; end
% plot the gaussian fits
tt=linspace(0,2*pi,npts)';
x = cos(tt); y=sin(tt);
ap = [x(:) y(:)]';
[v,d]=eig(C); 
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means, 1, size(ap,2)); 
h = plot(bp(1,:), bp(2,:), '-', 'parent', axh);
fill(bp(1,:), bp(2,:),color,'FaceAlpha',0.3,'LineStyle','none');

%-----------------------------
function h = show3d(means, C, color,sdwidth, npts, axh)
if isempty(npts), npts=20; end
[x,y,z] = sphere(npts);
ap = [x(:) y(:) z(:)]';
[v,d]=eig(C); 
% edit by Jiaheng Zhao
v = real(v);
d = real(d);
if any(d(:) < 0)
   fprintf('warning: negative eigenvalues\n');
   d = max(d,0);
end
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means, 1, size(ap,2)); 
xp = reshape(bp(1,:), size(x));
yp = reshape(bp(2,:), size(y));
zp = reshape(bp(3,:), size(z));
h = surf(axh, xp,yp,zp);
