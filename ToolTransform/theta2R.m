function R = theta2R(theta)
% For coordinate transformation, x2 = R*x1 + t set the contour clockwise as the positive direction.
R = [cos(theta) -sin(theta);
    sin(theta) cos(theta)];
end