function [Xstate,Zstate,lfeaID] = lmapIni(CenterState,ScanState,control,fsN)
% LMAPINI initiate local map
% Author: Jiaheng Zhao

% step number.
step = size(control,1)/3;
global scanPara
initial = [0;0;0]; % all start at [0;0;0;]
ZstateOdom=control;

for i=1:step
    if i == 1         % [0 0 0] is not included
        Xstate_pose(1:3,1) = Motionmodel(initial,control(3*i-2:3*i,:),[]);
    else
        Xstate_pose(3*i-2:3*i,1) = Motionmodel(Xstate_pose(3*i-5:3*i-3,1),control(3*i-2:3*i,:),[]);
    end
        Xstate_pose(3*i-2:3*i,2) = 1; % pose state is 1
        Xstate_pose(3*i-2:3*i,3) = i; % pose id
end
Xstate_pose(:,1)=wrapX(Xstate_pose(:,1),Xstate_pose);% Wrap theta of Xstate.
CenterState(:,4) = CenterState(:,4) - control(1,4);
ZstateOdom(:,3) = control(:,3) - control(1,4);
ZstateOdom(:,4) = control(:,4) - control(1,4);
[Xstate_feature,lfeaID] = invFuncfFS(Xstate_pose,CenterState,ScanState,scanPara,[]); % all feature are in the initial frame.
Xstate = [Xstate_pose; Xstate_feature];
Zstate.odom = ZstateOdom;
Zstate.center = CenterState;
Zstate.pts = ScanState;
Zstate.fsN = fsN;

end

function X = wrapX(X,Xstate)
id_pose = find(Xstate(:,2)==1); % pose id.
num_pose = size(id_pose,1)/3;
id_theta=repmat([0;0;1],num_pose,1);
id_theta=logical(id_theta);
theta = X(id_theta);
X(id_theta)=wrapToPi(theta);
end