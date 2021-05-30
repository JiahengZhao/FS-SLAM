function [pose_end,CenterState,ScanState,anbnState] = transFrameFS(Xstate,ScanState_pre,fsN)
% TRANSFRAMEFS transfer zstate to the last pose frame
% Author: Jiaheng Zhao
pose = Xstate(Xstate(:,2)==1,:);
step = max(pose(:,3));
pose_end=pose(end-2:end,1);
R_end = theta2R(pose_end(3));
feaState = Xstate(Xstate(:,2)==2,:);
nfea = size(feaState,1)/(2*fsN+3);
% scan pts
ScanState_tmp = ScanState_pre;
ScanState = cell(1, size(ScanState_pre,2));
for j = 1:step
    PoseState_ThisStep = pose(3*j-2:3*j,1);
    phi = PoseState_ThisStep(3);
    R = theta2R(phi);
    transFun = @(x) R_end'*(R*x + PoseState_ThisStep(1:2)) - R_end'*pose_end(1:2);
    ScanState_tmp(j,:) = cellfun(transFun,ScanState_pre(j,:),'UniformOutput',false);
end
for k = 1:size(ScanState_pre,2)
    % pts
    ScanState{1,k} = cat(2,ScanState_tmp{:,k});
end
% center, anbn
m = 2*fsN+3;
CenterState = zeros(2*nfea,4);
CenterState(:,2) = 2;
anbn_pre = zeros(2*fsN+1,nfea);
anbn = anbn_pre;
anbnState = zeros(length(anbn_pre(:)),4);
anbnState(:,2)=4;
for i = 1: nfea
    % center
    CenterState(2*i-1:2*i,1) = R_end' * feaState(m*i - (m-1):m*i - (m-2),1) - R_end' * pose_end(1:2);
    CenterState(2*i-1:2*i,3) = feaState(m*i - (m-1):m*i - (m-2),3);
    %anbn
    anbn_pre(:,i) = feaState(m*i - (m-3):m*i,1);
    % test transrform anbn
    for k = 1:fsN
        Rk = theta2R(k*pose_end(3));
        anbn([k,k+fsN+1],i)  =  Rk' * anbn_pre([k,k+fsN+1],i);
    end
    anbnState((2*fsN+1)*i-(2*fsN):(2*fsN+1)*i,1) = anbn(:,i);
    anbnState((2*fsN+1)*i-(2*fsN):(2*fsN+1)*i,3) = feaState(m*i - (m-1),3);
end
end