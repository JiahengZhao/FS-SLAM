function data = transScan(scan,pose)
t = pose(1:2);
theta = pose(3);
R = theta2R(theta);
% data = repmat(t,1,size(scan,2)) + R * scan;
data = t + R * scan;

end