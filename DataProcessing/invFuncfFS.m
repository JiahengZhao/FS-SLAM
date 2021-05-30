function  [Xstate_Feature,newid_seeFea] = invFuncfFS(Xstate_pose,obsState,scanState,scanPara,feaOccurredID)
% INVFUNCF solve initial Xstate_Feature via given ObsThisStep.
% OUTPUT:
%             Xstate_Feature returns xstate for feature.
%             newid_seeFea returns all of the seen featurns id with new order (only for global)
% Noted newid_seeFea should have different column dimension against feaOccurredID
if nargin ~= 6
    CovScan = [];
end
if isempty(obsState)
    Xstate_Feature = [];
    newid_seeFea = 0;
    return
end
if isstruct(scanPara)
    fsN_f = scanPara.fsN_local;
    fsN_border = scanPara.fsN_border;
    fsN_rect = scanPara.fsN_rect;
else
    fsN_f = scanPara;
    fsN_border = 0;
    fsN_rect = 0;
end

if isempty(feaOccurredID)
    %------------------For Global solving------------------%
    step = size(Xstate_pose,1)/3;
    id_seeFea = unique(obsState(:,3));
    newid_seeFea = [(1:size(id_seeFea))', id_seeFea];
    %------rectF and borderF has different fsN -----%
    nF = size(find(id_seeFea<10),1);
    nB = size(find(id_seeFea>9 & id_seeFea < 100),1);
    nR = size(find(id_seeFea>99),1);
    nxf = (3+2*fsN_f)*nF + (3+2*fsN_border)*nB + (3+2*fsN_rect)*nR;
    Xstate_Feature = zeros(nxf,3); % cetner + an bn
    Xstate_Feature(:,2)=2; % feature is 2;
    tmpZstate = cell(step,size(scanState,2));
    for j=1:step
        obs_ThisStep = obsState(obsState(:,4)==j,:);
        id_seeFea_ThisStep = unique(obs_ThisStep(:,3));
        pose_ThisStep = Xstate_pose(3*j-2:3*j,1);
        phi = pose_ThisStep(3);
        R = theta2R(phi);
        t = pose_ThisStep(1:2);
        transFun = @(x) R*x + t;
        tmpZstate(j,:) = cellfun(transFun,scanState(j,:),'UniformOutput',false);
        for i = 1:size(id_seeFea_ThisStep,1)
            Feeindex_temp = newid_seeFea(newid_seeFea(:,2)==id_seeFea_ThisStep(i),1);
            Feeindex_old = newid_seeFea(newid_seeFea(:,2)==id_seeFea_ThisStep(i),2);
            if Feeindex_old > 9 && Feeindex_old < 100
                fsN = fsN_border;
                Feeindex_old_1 = Feeindex_old/10;
                Feeindex_old = Feeindex_old_1 + nF;
                xid = (3+2*fsN_f)*nF + (3+2*fsN)*Feeindex_old_1-(2+2*fsN) : (3+2*fsN_f)*nF ...
                    + (3+2*fsN)*Feeindex_old_1;
            elseif Feeindex_old > 99
                fsN = fsN_rect;
                Feeindex_old_1 = Feeindex_old/100;
                Feeindex_old = Feeindex_old_1 + nF + nB;
                xid = (3+2*fsN_f)*nF + (3+2*fsN_border)*nB  + (3+2*fsN)*Feeindex_old_1-(2+2*fsN) : (3+2*fsN_f)*nF ...
                    + (3+2*fsN_border)*nB + (3+2*fsN)*Feeindex_old_1;
            else
                fsN = fsN_f;
                xid = (3+2*fsN)*Feeindex_temp-(2+2*fsN):(3+2*fsN)*Feeindex_temp;
            end
            if  Xstate_Feature(xid,1) == 0 
                scan_temp = transScan(scanState{j,Feeindex_old},pose_ThisStep);
                if isempty(scan_temp)
                    continue;
                end
                [V, ~, ~] = fitWithFS( fsN, scan_temp(1,:),scan_temp(2,:), [],true); % isAdd sometimes cause trouble. 
                if isempty(V)
                    continue;
                end
                Xstate_Feature(xid(1:2),1) = R*obs_ThisStep(2*i-1:2*i,1)+t;
                Xstate_Feature(xid(3:end),1) = V;
                Xstate_Feature(xid,3) = Feeindex_temp;
            end
        end
    end
end



