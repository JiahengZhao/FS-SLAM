function [ang, center, isB] = calVarCur(thispt)
% Author: Jiaheng Zhao.
% Calculate the angle betwen the normal line in both end
%                        _________
%                       /              \  <----|          For example, this two arrow (there should be an angle, but I 
%        ---->      |                                       don't know how to draw it in comment.       ^_^

ptFitted = thispt;
N = LineNormals2D(ptFitted');
t1 = N(1,:) +N(2,:)+N(3,:);
t2 = N(end,:)+N(end-1,:)+N(end-2,:);

v1 = [t1/norm(t1) 0];
v2 = [t2/norm(t2) 0];

ang = atan2(norm(cross(v1,v2)), dot(v1,v2));

isB = false;
% divariaty
sp = ptFitted(:,1:3); sp = sum(sp,2)/3;
ep = ptFitted(:,end-2:end); ep = sum(ep,2)/3;

cn = v1(2)/v1(1) - v2(2)/v2(1);
if abs(cn) > 1e-6
    nd = ep(2) - sp(2) + (v1(2)/v1(1)*sp(1) - v2(2)/v2(1) * ep(1));
    x = nd / cn;
    y = (x - sp(1)) *  v1(2)/v1(1)  + sp(2);
    center = [x y];
    if norm(max(ptFitted,[],2) - min(ptFitted,[],2)) > 10 % if the size of points is larger than 10, then it's border
        isB = true;
    end
else
    if size(ptFitted,2) >40
        isB = true;
        center = [];
    else
        center = [];
    end
end

end
