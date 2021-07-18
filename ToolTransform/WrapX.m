function X = WrapX(X,Xstate)
id_pose = find(Xstate(:,2)==1); % pose id.
id_feature = find(Xstate(:,2)==2);
num_pose = size(id_pose,1)/3;
num_fea = size(id_feature,1)/5;
id_theta=[repmat([0;0;1],num_pose,1);repmat([0;0;1;0;0],num_fea,1)];
id_theta=logical(id_theta);
theta = X(id_theta);
X(id_theta)=wrapToPi(theta);
end