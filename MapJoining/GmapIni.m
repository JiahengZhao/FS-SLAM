function [Xstate,Zstate,lfeaID] = GmapIni(lmap,gfsN)
% GMAPINI initiate global map
% Author: Jiaheng Zhao

 centerState = [];
 scanState = [];
 anbnState = [];
% step number.
step = length(lmap);
control = cat(2,lmap.endPose);
initial = [0;0;0]; % all start at [0;0;0;]
ZstateOdom=zeros(3*step,4);

for i=1:step
    if i == 1         % [0 0 0] is not included
        Xstate_pose(1:3,1) = Motionmodel(initial,control(:,i),[]);
    else
        Xstate_pose(3*i-2:3*i,1) = Motionmodel(Xstate_pose(3*i-5:3*i-3,1),control(:,i),[]);
    end
        Xstate_pose(3*i-2:3*i,2) = 1; % pose state is 1
        Xstate_pose(3*i-2:3*i,3) = i; % pose id
        ZstateOdom(3*i-2:3*i,:) = [control(:,i) ones(3,1)  i*ones(3,1) (i-1)*ones(3,1)];
        
        centerState = [centerState;lmap(i).Zstate.center];
        scanState = [scanState;lmap(i).Zstate.pts];
        anbnState = [anbnState;lmap(i).Zstate.anbn];
end
Xstate_pose(:,1)=wrapX(Xstate_pose(:,1),Xstate_pose);% Wrap theta of Xstate.

[Xstate_feature,lfeaID] = invFuncfFS(Xstate_pose,centerState,scanState,gfsN,[]); % all feature are in the initial frame.
Xstate = [Xstate_pose; Xstate_feature];
Zstate.odom = ZstateOdom;
Zstate.center = centerState;
Zstate.pts = scanState;
Zstate.fsN = gfsN;

end

function X = wrapX(X,Xstate)
id_pose = find(Xstate(:,2)==1); % pose id.
num_pose = size(id_pose,1)/3;
id_theta=repmat([0;0;1],num_pose,1);
id_theta=logical(id_theta);
theta = X(id_theta);
X(id_theta)=wrapToPi(theta);
end
