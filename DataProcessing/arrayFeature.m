function  [feaOccurredID,invCovFeaOutput,newObsState] = arrayFeature(obsState,invCovFea)
% ARRAYFEATURE  determine the number of seen feature and return to feaOccurredID.
% Also return reordered inverse covariance of feature part

% feaOccurredID with format: [newID   preID   occuredStepID]
% Example: if at the first step see feature 1 3 4, then feaOccurredID is
%                                        1   1   1
%                                        2   3   1
%                                        3   4   1
% The second step see feature 1 4 5, feaOccurredID is
%                                        1   1   1
%                                        2   3   1
%                                        3   4   1
%                                        4   5   2
%
if isempty(obsState)
    feaOccurredID = [0 0 0];
    invCovFeaOutput = invCovFea;
    newObsState = [];
    return
end
if isempty(invCovFea) % for test new algorithm
    newObsState = obsState;
    nStep = max(cell2mat(obsState(:,4)));

    for i = 1: nStep
        thisObs = obsState(cell2mat(obsState(:,4)) == i,:); % Observation at this step
        thisOriFeaID = unique(cell2mat(thisObs(:,3))); % Original feature id at this step
        nThisOriFeaID = size(thisOriFeaID,1); % number of this feature
      
        if i == 1
            feaOccurredID = [ (1:nThisOriFeaID)',  thisOriFeaID, i.*ones(size(thisOriFeaID,1),1)];
            for j = 1:nThisOriFeaID
                newObsState{(cell2mat(obsState(:,4))==i) ... 
                    & (cell2mat(obsState(:,3))==thisOriFeaID(j)),3} =...
                    feaOccurredID(feaOccurredID(:,2) == thisOriFeaID(j),1);
            end
        else
            [liA,~]=ismember(thisOriFeaID,feaOccurredID(:,2)); % Judge whether occurs this observation
            idNewFea = find(liA == 0);
            
            if ~isempty(idNewFea) %  If new feature occurs
                preID = size(feaOccurredID,1);
                tmpID = 1:size(idNewFea,1); % update new feature's id
                tmpID = preID + tmpID; % Join
                feaOccurredID = [ feaOccurredID;
                    tmpID,  thisOriFeaID(idNewFea), i.*ones(size(idNewFea,1),1)];% New feature
            end

            for j = 1:nThisOriFeaID
                newFeaID = feaOccurredID(feaOccurredID(:,2) == thisOriFeaID(j),1); % find new id
                newObsState{(cell2mat(obsState(:,4))==i) & ...
                    (cell2mat(obsState(:,3))==thisOriFeaID(j)),3}=newFeaID;
            end

        end
        
    end

    invCovFeaOutput = invCovFea;
    
else
    %%
    newObsState = obsState;
    nStep = max(obsState(:,4));
    % tmpThisInvCov = [];
    for i = 1: nStep
        thisObs = obsState(obsState(:,4) == i,:); % Observation at this step
        thisOriFeaID = unique(thisObs(:,3)); % Original feature id at this step
        nThisOriFeaID = size(thisOriFeaID,1); % number of this feature
        tmpInvCov = invCovFea(invCovFea(:,7)==i,:); % temp invCovFea of this step
        
        
        if i == 1
            feaOccurredID = [ (1:nThisOriFeaID)',  thisOriFeaID, i.*ones(size(thisOriFeaID,1),1)];
            for j = 1:nThisOriFeaID
                %             tmpThisInvCov = [tmpThisInvCov;
                %                 tmpInvCov(tmpInvCov(:,6)==thisOriFeaID(j),:)];
                %             tmpThisInvCov(5*j-4:5*j,6) = feaOccurredID(feaOccurredID(:,2) == thisOriFeaID(j),1);
                newObsState((newObsState(:,4)==i) == (newObsState(:,3)==thisOriFeaID(j)),3)=feaOccurredID(feaOccurredID(:,2) == thisOriFeaID(j),1);
            end
        else
            [liA,locB]=ismember(thisOriFeaID,feaOccurredID(:,2)); % Judge whether occurs this observation
            idNewFea = find(liA == 0);
            
            if ~isempty(idNewFea) %  If new feature occurs
                preID = size(feaOccurredID,1);
                tmpID = 1:size(idNewFea,1); % update new feature's id
                tmpID = preID + tmpID; % Join
                feaOccurredID = [ feaOccurredID;
                    tmpID',  thisOriFeaID(idNewFea), i.*ones(size(idNewFea,1),1)];% New feature
            end
            %         tmptmpThisInvCov = zeros(5*size(feaOccurredID,1),7);
            for j = 1:nThisOriFeaID
                newFeaID = feaOccurredID(feaOccurredID(:,2) == thisOriFeaID(j),1); % find new id
                %             tmptmpThisInvCov(5*newFeaID-4:5*newFeaID,:) = tmpInvCov(tmpInvCov(:,6)==thisOriFeaID(j),:);
                %             tmptmpThisInvCov(5*newFeaID-4:5*newFeaID,6) = newFeaID;
                newObsState((newObsState(:,4)==i) == (newObsState(:,3)==thisOriFeaID(j)),3)=newFeaID;
            end
            %         tmptmpThisInvCov(all(tmptmpThisInvCov==0,2),:)=[];
            %         tmpThisInvCov = [tmpThisInvCov; tmptmpThisInvCov];
        end
        
    end
    % invCovFeaOutput = tmpThisInvCov;
    invCovFeaOutput = invCovFea;
end
end