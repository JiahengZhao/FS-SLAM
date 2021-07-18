function [Xstate,Zstate,invCovMatrix,feaOccurredID] = ...
    InitializeFSProb(CenterState,ScanState,control,invCovFea)
% Author: Jiaheng Zhao
% INITIALIZEPROB Initialization Xstate and Zstate.
%                           1                    2                     3
%       Xstate = [ value, pose->1 feature->2, id_this]
%                           1                    2                     3            4
%       Zstate = [ value, pose->1 feature->2, id_this, id_relativeto]
%
global scanPara
step = size(control,1);
Cov_pose = scanPara.covodom; % in vector
invCov_pose = inv(Cov_pose);
fsN = scanPara.fsN_local;
fsN_rect = scanPara.fsN_rect;
fsN_border = scanPara.fsN_border;
if nargin < 5 % For global calculating
    [invCovMatrix] =  mergCov(invCov_pose, step,invCovFea);
else
    invCovMatrix = [];
end

initial = [0;0;0]; % all start at [0;0;0;]
ZstateOdom=zeros(3*step,4);

for i=1:step
    if i == 1         % [0 0 0] is not included
        Xstate_pose(1:3,1) = Motionmodel(initial,control(i,:),[]);
    else
        Xstate_pose(3*i-2:3*i,1) = Motionmodel(Xstate_pose(3*i-5:3*i-3,1),control(i,:),[]);
    end
        Xstate_pose(3*i-2:3*i,2) = 1; % pose state is 1
        Xstate_pose(3*i-2:3*i,3) = i; % pose id
        %4 columns: value, 1 for odom 2 for feature, id, reference pose id
        ZstateOdom(3*i-2:3*i,:) = [control(i,:)' ones(3,1)  i*ones(3,1) (i-1)*ones(3,1)];
end
[Xstate_feature,feaOccurredID] = invFuncfFS(Xstate_pose,CenterState,ScanState,scanPara,[]);
Xstate_pose(:,1)=wrapX(Xstate_pose(:,1),Xstate_pose);% Wrap theta of Xstate.
Xstate = [Xstate_pose; Xstate_feature];
Zstate.odom = ZstateOdom;
Zstate.center = CenterState;
Zstate.pts = ScanState;
Zstate.fsN = fsN;
Zstate.fsN_rect = fsN_rect;
Zstate.fsN_border = fsN_border;

end

function X = wrapX(X,Xstate)
id_pose = find(Xstate(:,2)==1); % pose id.
num_pose = size(id_pose,1)/3;
id_theta=repmat([0;0;1],num_pose,1);
id_theta=logical(id_theta);
theta = X(id_theta);
X(id_theta)=wrapToPi(theta);
end


function [invcov] = mergCov(invCovPose,num_pose,invCovFea)
% MERGCOV align CovPose and CovFea into diagonal covariance matrix
% invCovPose is 3 by 3 diagonal matrix.
% invCovFea is M by 7.
nRowFea = size(invCovFea,1);
nRowPose = 3*num_pose;
% n = nRowPose +nRowFea;
% m= 3*num_pose +length(CovFea)*0.2*3;
% cov = sparse(n,n);
pcov = diag(invCovPose);
pcov = repmat(pcov,num_pose,1); % repmat cov pose part to vector
% cov(1:3*num_pose,1:3*num_pose) = diag(pcov);
% Sparse matrix. ID of matrix: from up to down
poseid1 = 1:nRowPose;
poseid2 = poseid1;
% cov = sparse(id1,id2,pcov); % pose
tmp = 1 : nRowFea; %(3*num_pose+1) : n;
feaid2 = repmat((tmp+nRowPose),5,1);
feaid2 = feaid2(:);
feaid1 = zeros(5*nRowFea,1);% preallocated
tmpidSmall = repmat(1:5,1,5);
tmpidSmall = tmpidSmall(:)+nRowPose;
feaval = feaid1;
for j = 1: (nRowFea*0.2)
    %     cov((3*num_pose-4+5*j):(3*num_pose+5*j),(3*num_pose-4+5*j):(3*num_pose+5*j)) = invCovFea(5*j-4:5*j,1:5);
    tmpCov = invCovFea(5*j-4:5*j,1:5);
    feaval(25*j-24:25*j,1) = tmpCov(:);
    feaid1(25*j-24:25*j,1) = tmpidSmall + 5*(j-1);
end
id1 = [poseid1';feaid1];
id2 = [poseid2';feaid2];
val = [pcov;feaval];
invcov = sparse(id1,id2,val);

end
