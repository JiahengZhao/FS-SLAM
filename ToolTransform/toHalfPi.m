function angle = toHalfPi(angle)
% This function limits angle into [-pi/2  pi/2)
%       First, wrap angle into [-pi pi]
%       Then, limit to half pi. Noted angle doesn't equal pi/2

q = (angle < -pi) | (pi < angle);
angle(q) = To2Pi(angle(q) + pi) - pi;

p1 = (angle >= pi/2);
p2 =(angle < -pi/2);
angle(p1) = angle(p1) - pi;
angle(p2) = angle(p2) +pi;

end

function lambda = To2Pi(lambda)
positiveInput = (lambda > 0);
lambda = mod(lambda, 2*pi);
lambda((lambda == 0) & positiveInput) = 2*pi;
end