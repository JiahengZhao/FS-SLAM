function [x0, y0, a, b, angle] = Matrix2AngleForm(A, B, C, D, E, F)
%MATRIX2ANGLEFORM Computes the angle form of matrix form of an ellipse
%  See: http://en.wikipedia.org/wiki/Matrix_representation_of_conic_sections#Reduced_equation
%
% AUTHOR Sebastian Dingler <s.dingler@gmail.com>
%        Karlsruhe Institute of Technology (KIT), Germany
%
% DATE   22.12.2014

B = 2*B;            %----------------> Comment
D = 2*D;            %----------------> Comment
E = 2*E;            %----------------> Comment
AQ = [A B/2 D/2;B/2 C E/2;D/2 E/2 F];           %----------------> Comment
A33 = [A B/2;B/2 C];                                %----------------> Comment
% Center
center = pinv(A33)*[-D/2;-E/2];             %----------------> Comment
x0 = center(1);                     %----------------> Comment
y0 = center(2);                     %----------------> Comment
% Major and minor axis
e = eig(A33);                       %----------------> Comment
[V,D] = eig(A33);                  %----------------> Comment
% [u,s,v] = svd(A33)
% if det(u) < 0 
%     u(2,;) = -u(2,:);
% end
% 
% if det(V) < 0 
%     V(:,1) = -V(:,1);
% end
a = sqrt(-det(AQ)/(det(A33)*e(1)));                 %----------------> Comment
b = sqrt(-det(AQ)/(det(A33)*e(2)));                     %----------------> Comment
% angle
angle = atan2(-V(2,1),-V(1,1)); % This is right! Corrected by Jiaheng Zhao, Shoudong Huang, Liang Zhao   %----------------> Comment
% angle = acos(-dot(V(:,1),[1;0])/norm(V(:,1)));

end