function [linePloter,borderROI] = drawLine(LineXstate,isDrawRange,lineID,color,drift)
if nargin < 4
    color = 'r';
    drift = [];
end
if isDrawRange
    borderROI = getrect;
    borderROI = [borderROI(1) borderROI(1)+borderROI(3) borderROI(2) borderROI(2)+borderROI(4)];
else
%     borderROI = [-0.6 14.6 -4.1 4.1];
    borderROI = [-20 20 -10 105];
%     lineRange = [-0.5 -0.5 -4 4;
%         -0.5 14.5 4 4;
%         -0.5 14.5 -4 -4;
%         14.5 14.5 -4 4;
%         8 8 1 4;
%         8 8 -4 -1;
%         -0.5 2.5 1 1;
%         4 8 1.5 1.5];% In old order. xmin xmax ymin ymax norm_angle;
%     borderROI = [-20 25 -8 18];
end
% yy = borderROI(3):0.01:borderROI(4);
% xx = borderROI(1):0.01:borderROI(2);
% ax = gca;
axis(borderROI);
hold on;
for i = 1:max(LineXstate(:,3))
    thisline = LineXstate(LineXstate(:,3)==i,1);
    linePloter{i} =  refline(-cot(thisline(1)),thisline(2)/sin(thisline(1)));
    set(linePloter{i},'Color',color);
    set(linePloter{i},'LineWidth',3);
%     refline(-cot(x(1)),x(2)/sin(x(1)));hold on; set(l1,'Color','r');
    
    %%
%     thisLine = [cos(thisline(1)); sin(thisline(1)); -thisline(2)];
%     if ~isempty(drift)
%         T = [1 0 0; 0 1 0; drift' 1];
%         thisLine = T * thisLine;
%     end
%     if (abs(thisline(1)) < 1*pi/180) || (abs(thisline(1)-pi) < 1*pi/180) ||  (abs(thisline(1)+pi) < 1*pi/180)% x
%         x = -(thisLine(2).*yy+thisLine(3))./thisLine(1); 
%         valid = (x<=borderROI(2)+0.1) & (x >= borderROI(1)-0.1) ;
%         if ~isDrawRange
%             ori_lineID = lineID(lineID(:,1) == i,2);
%             lineRange(ori_lineID,:);
% %             avoid = 
%         end
%         
%         if ~isempty(find(valid,2))
%             linePloter(i) = plot(x(valid), yy(valid), 'Color',color,'LineWidth',2);
%         else
%             fprintf('Line is over sized. ID: %s.\n',i);
%         end
%     else
%         y =  -(thisLine(1).*xx+thisLine(3))./thisLine(2);
%         valid = (y<=borderROI(4)+0.1) & (y >= borderROI(3)-0.1) ;
%         if ~isempty(find(valid,2))
%             linePloter(i) = plot(xx(valid), y(valid), 'Color',color,'LineWidth',2);
%         else
%             fprintf('Line is over sized. ID: %s.\n',i);
%         end
%     end
end
axis equal 
axis(borderROI);
end