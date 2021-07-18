function p=plotFeature(F,color,flag,FeaID)
if nargin < 3
    flag = 'default';
    FeaID = [];
end
if isempty(F)
    p=[];
    return 
end
[~,num] = size(F);
% Fpoint = F(1:2,:); % 2 by n;
% p.center = plot(Fpoint(1,:),Fpoint(2,:),'o','color',color); 
hold on;
for i = 1:num
    % 95% confidence level corresponds to s^2 = 5.991
%     s = sqrt(5.991);
    s=1;
    mu = F(1:2,i)';
    phi = F(3,i);
    R = theta2R(phi);
    r1 = F(4,i);
    r2 = F(5,i);
%     cov_iden = [r1^2 0; 0 r2^2];
%     cov_real = R*cov_iden*R'; 
    
    % change the plot part
    if ~exist('sdwidth', 'var'), sdwidth = 1; end
    if ~exist('npts', 'var'), npts = 50; end
    if ~exist('axh', 'var'), axh = gca; end
    % plot the gaussian fits
    tt=linspace(0,2*pi,npts)';
    x = cos(tt); y=sin(tt);
    ap = [x(:) y(:)]';
%     [v,d1]=eig(cov_real);
    d = sdwidth *[r1 0;0 r2 ];
%     d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
    bp = (R*d*ap) + repmat(mu', 1, size(ap,2));
    p(i) = plot(bp(1,:), bp(2,:), '-', 'parent', axh);
%     p(i)=Plot_gaussian_ellipsoid(mu,cov_real,s);
    %     set(p.ellipse{i},'color',color,'LineWidth',2);
    set(p(i),'color',color,'LineWidth',3);
    switch flag
        case 'default'
            text(mu(1)+0.3,mu(2)+0.1,num2str(i));
        case 'seeFeature'    % seen feature in each step
            text(mu(1)+0.3,mu(2)-0.5,FeaID);
        case 'notext'
           
        case 'gt'
             set(p(i),'LineStyle','--');
    end
    
    arrowLength=0;
%     dx=arrowLength*cos(phi)+mu(1);
%     dy=arrowLength*sin(phi)+mu(2);
%     arrow fixlimits;   
% Instead of arrow-----draw groundtruth
    
    dx = mu(1)+ arrowLength * cos(phi);
    dy = mu(2)+ arrowLength * sin(phi);
%     p.arrow{i}=plot([mu(1);dx],[mu(2);dy],'Color',color,'LineWidth',2);
%     arrow(mu,[dx, dy],'EdgeColor',color,'FaceColor',color,'TipAngle',10,'Length',10,'BaseAngle',20);
%     annotation('arrow',[mu(1);mu(1)+dx],[mu(2);mu(2)+dy],'Color',color);
%     quiver(mu(1),mu(2),mu(1)+dx,mu(2)+dy,0.1,'Color',color);
    
end
hold off;
end