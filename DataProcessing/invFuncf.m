function  [Xstate_Feature,newid_seeFea] = invFuncf(Xstate_pose,obsState,feaOccurredID,CovScan)
% INVFUNCF solve initial Xstate_Feature via given ObsThisStep.
% OUTPUT:
%             Xstate_Feature returns xstate for feature.
%             newid_seeFea returns all of the seen featurns id with new order (only for global)
% Noted newid_seeFea should have different column dimension against feaOccurredID
if nargin ~= 4
    CovScan = [];
end
if isempty(obsState)
    Xstate_Feature = [];
    newid_seeFea = 0;
    return
end

if ~iscell(obsState)
    if isempty(feaOccurredID)
        %------------------For Global solving------------------%
        step = size(Xstate_pose,1)/3;
        id_seeFea = unique(obsState(:,3));
        newid_seeFea = [(1:size(id_seeFea))', id_seeFea];
        if step == 1
            phi = Xstate_pose(3);
            R = theta2R(phi);
            t = Xstate_pose(1:2);
            for i = 1:size(id_seeFea,1)
                Xstate_Feature(5*i-4:5*i-3,1) = R*obsState(5*i-4:5*i-3,1)+t;
                Xstate_Feature(5*i-2,1) =obsState(5*i-2,1) +phi;
                Xstate_Feature(5*i-2,1) = wrap(Xstate_Feature(5*i-2,1)); % wrap
                Xstate_Feature(5*i-1:5*i,1)=obsState(5*i-1:5*i,1);
            end
        else
            Xstate_Feature = zeros(5*size(id_seeFea,1),3);
            Xstate_Feature(:,2)=2; % feature is 2;
            for j=1:step
                obs_ThisStep = obsState(obsState(:,4)==j,:);
                id_seeFea_ThisStep = unique(obs_ThisStep(:,3));
                pose_ThisStep = Xstate_pose(3*j-2:3*j,1);
                phi = pose_ThisStep(3);
                R = theta2R(phi);
                t = pose_ThisStep(1:2);
                for i = 1:size(id_seeFea_ThisStep,1)
                    Feeindex_temp = newid_seeFea(newid_seeFea(:,2)==id_seeFea_ThisStep(i),1);
                    if  Xstate_Feature(5*Feeindex_temp-4:5*Feeindex_temp,1) == 0
                        Xstate_Feature(5*Feeindex_temp-4:5*Feeindex_temp-3,1) = R*obs_ThisStep(5*i-4:5*i-3,1)+t;
                        Xstate_Feature(5*Feeindex_temp-2,1) = obs_ThisStep(5*i-2,1) +phi;
                        Xstate_Feature(5*Feeindex_temp-2,1) = wrapToPi(Xstate_Feature(5*Feeindex_temp-2,1)); % wrap
                        Xstate_Feature(5*Feeindex_temp-1:5*Feeindex_temp,1)=obs_ThisStep(5*i-1:5*i,1);
                        Xstate_Feature(5*Feeindex_temp-4:5*Feeindex_temp,3)=Feeindex_temp;
                    end
                end
            end
        end
        %------------------for Incremental solving------------------%
    else
        nSeeFea = size(feaOccurredID,1); % number of feature
        Xstate_Feature = zeros(5*nSeeFea,3);
        Xstate_Feature(:,2)=2;
        for i = 1:nSeeFea
            poseIDSeeFea = feaOccurredID(i,3);
            obsThispose = obsState(obsState(:,4)==poseIDSeeFea,:); % All observations at this step.
            thisFeaObs =  obsThispose(obsThispose(:,3) == feaOccurredID(i,2),1);% The first time seeing this feature
            pose_ThisStep = Xstate_pose(3*poseIDSeeFea-2:3*poseIDSeeFea,1); % pose id of this feature
            phi = pose_ThisStep(3);
            R = theta2R(phi);
            t = pose_ThisStep(1:2);
            Xstate_Feature(5*i-4:5*i-3,1) = R * thisFeaObs(1:2,1) + t;
            Xstate_Feature(5*i-2,1) = phi + thisFeaObs(3,1);
            Xstate_Feature(5*i-1:5*i,1) =  thisFeaObs(4:5,1);
            Xstate_Feature(5*i-4:5*i,3) = i;
        end
        
    end
else
    if isempty(feaOccurredID)
        %------------------For Global solving------------------%
        step = size(Xstate_pose,1);
        id_seeFea = unique(cell2mat(obsState(cell2mat(obsState(:,3)),3)));
        newid_seeFea = [(1:size(id_seeFea))', id_seeFea];
        if step == 1
            phi = Xstate_pose(3);
            R = theta2R(phi);
            t = Xstate_pose(1:2);
            for i = 1:size(id_seeFea,1)
                Xstate_Feature(5*i-4:5*i-3,1) = R*obsState(5*i-4:5*i-3,1)+t;
                Xstate_Feature(5*i-2,1) =obsState(5*i-2,1) +phi;
                Xstate_Feature(5*i-2,1) = wrap(Xstate_Feature(5*i-2,1)); % wrap
                Xstate_Feature(5*i-1:5*i,1)=obsState(5*i-1:5*i,1);
            end
        else
            Xstate_Feature = zeros(5*size(id_seeFea,1),3);
            Xstate_Feature(:,2)=2; % feature is 2;
            for j=1:step
                obs_ThisStep = obsState(obsState(:,4)==j,:);
                id_seeFea_ThisStep = unique(obs_ThisStep(:,3));
                pose_ThisStep = Xstate_pose(3*j-2:3*j,1);
                phi = pose_ThisStep(3);
                R = theta2R(phi);
                t = pose_ThisStep(1:2);
                for i = 1:size(id_seeFea_ThisStep,1)
                    Feeindex_temp = newid_seeFea(newid_seeFea(:,2)==id_seeFea_ThisStep(i),1);
                    if  Xstate_Feature(5*Feeindex_temp-4:5*Feeindex_temp,1) == 0
                        Xstate_Feature(5*Feeindex_temp-4:5*Feeindex_temp-3,1) = R*obs_ThisStep(5*i-4:5*i-3,1)+t;
                        Xstate_Feature(5*Feeindex_temp-2,1) = obs_ThisStep(5*i-2,1) +phi;
                        Xstate_Feature(5*Feeindex_temp-2,1) = wrap(Xstate_Feature(5*Feeindex_temp-2,1)); % wrap
                        Xstate_Feature(5*Feeindex_temp-1:5*Feeindex_temp,1)=obs_ThisStep(5*i-1:5*i,1);
                        Xstate_Feature(5*Feeindex_temp-4:5*Feeindex_temp,3)=Feeindex_temp;
                    end
                end
            end
        end
        %------------------for Incremental solving------------------%
    else
        step = max(Xstate_pose(:,3));
        nSeeFea = size(feaOccurredID,1); % number of feature
        Xstate_Feature = zeros(5*nSeeFea,3);
        w = 4;
        feaCount = zeros(nSeeFea,1);
        for i = 1:nSeeFea
%             fmean=[];
            % Use the first time and next few poses this feature occured pose to initialize
            poseIDSeeFea = feaOccurredID(i,3);
            
            % Get a rough fitted parameter.
            isill = 1;
            while isill
                obsThispose = obsState(cell2mat(obsState(:,4))==poseIDSeeFea,:); % All observations at this step. Original ID
                thisFeaObs =  cell2mat(obsThispose(cell2mat(obsThispose(:,3)) == ...
                    feaOccurredID(i,2),1));% The first time seeing this feature
                if poseIDSeeFea > max(Xstate_pose(:,3))
                    error(" Can't find an initial fitted Ellipse (Feature original id is %d). Please check the data.",...
                        feaOccurredID(i,2));
                end
                if isempty(thisFeaObs)
                    break
                end
                pose_ThisStep = Xstate_pose(Xstate_pose(:,3)==poseIDSeeFea,1); % pose id of this feature
                phi = pose_ThisStep(3);
                R = theta2R(phi);
                t = pose_ThisStep(1:2);
                if ~isempty(thisFeaObs) && size(thisFeaObs,2) > 10
                    [ang] = calVarCur(thisFeaObs);
                    if  ang>=0.4  % in order to get over complex fitted result.
                        [f,invcovf,isill] = fitWithLSO(thisFeaObs(1,:),thisFeaObs(2,:), CovScan,phi);
                        if ~isreal(invcovf)
                            disp(['Find Complex Fitted Ellipse via Feature ',...
                                num2str(i),'at Pose ',num2str(poseIDSeeFea),'. Move to next pose...']);
                            poseIDSeeFea = poseIDSeeFea+1;
                            continue;
                        elseif isill == 1
                            disp(['Find Complex Fitted Ellipse via Feature ',...
                                num2str(i),'at Pose ',num2str(poseIDSeeFea),'. Move to next pose...']);
                            poseIDSeeFea = poseIDSeeFea+1;
                            continue;
                        elseif isill == 2
                            disp(['Find ill conditioned Ellipse via Feature ',...
                                num2str(i),'at Pose ',num2str(poseIDSeeFea),'. Move to next pose...']);
                            poseIDSeeFea = poseIDSeeFea+1;
                            continue;
                        end
                        % Can't take average directly. Error.
%                         feaCount(i) = feaCount(i) + 1;
%                         if feaCount(i) <= w && poseIDSeeFea < step
%                             isill = 1;
%                             poseIDSeeFea = poseIDSeeFea+1;
%                             fmean = [fmean f];
%                         end
                        
                    end
                end
            end
%             f = mean(fmean,2);
            Xstate_Feature(5*i-4:5*i,1) = [R * f(1:2,1) + t;phi + f(3,1); f(4:5,1)];
            Xstate_Feature(5*i-4:5*i,2) = 2;
            Xstate_Feature(5*i-4:5*i,3) = i;
        end
    end
end


