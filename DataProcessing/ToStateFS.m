function fs= ToStateFS(scan,scanPara)

% main axis should always be biger than minor axis.
% Setting all features' angle in global frame varies from 0 to Pi.
if isfield(scan,'rectF')
    nR = size(scan.rectF,2);
else
    nR = 0;
end
if isfield(scan, 'borderF')
    nB = size(scan.borderF,2);
else 
    nB = 0;
end
nF = size(scan.feature,2);
numFea = (nF+nR+nB);
step = size(scan.scan,2);
% step=10;
obs=[];
% CovScan = scanPara.covscan; % Defined by scan points x, y;
invCovFea = [];

%-------------Compare center to filt feature--------%
rCheck = zeros(2,numFea);
rUncRange=rCheck;
% This is a test for sliding window: use previous w fitted r1 r2 as a criteria
w=15; % sliding window size
swCenter = cell(numFea,1);% save r1 r2

%---------------------------------------------%
feaCount= zeros(numFea,1); % count feature occured times.
ptObs = cell(step,numFea);
disp('Start processing data...');
for i = 1:step % step id
    Thisscan = scan.scan{i};
    if isempty(Thisscan)
        continue;
    end

    %%
    for j=1:numFea
        if j>nF && j <= nF+nB
            fID = 10*(j-nF);
        elseif j > nF + nB
            fID = 100*(j-nF-nB);
        else
            fID = j;
        end
            FifFeaPtc = Thisscan(2:3,Thisscan(1,:)==fID);
            ptObs{i,j} = FifFeaPtc;
            % when points number > 7, processing
            if ~isempty(FifFeaPtc) && size(FifFeaPtc,2) > 7
                [ang,check_center,isB] = calVarCur(FifFeaPtc);
                if  ang>=0.4  % in order to get over complex fitted result.
                    [center,r]= fitWithCircle(FifFeaPtc(1,:),FifFeaPtc(2,:));
                    if ~isB && norm(check_center'-center) >1.5 
                        continue;
                    end
                    if isB && r > 40
                        continue;
                    end
                    feaCount(j) = feaCount(j) + 1;
                    if feaCount(j) <=w
                        swCenter{j} = [swCenter{j} r];
                    else
                        rCheck(:,j) = trimmean(swCenter{j}(:,1:w),40,2);
                        rUncRange(:,j) = abs((r-rCheck(:,j)))./r;
                        if rUncRange(1,j) > 1 || rUncRange(2,j) > 1
                            feaCount(j) = feaCount(j) - 1;
                            disp(['                Step ', num2str(i), ': Feature ',num2str(j),': radius beyond range. Average radius is ',...
                                num2str(rCheck(1,j)), '. Current radius is :',num2str(r), '. Ignored...']);
                            continue;
                        end
                        swCenter{j}(:,1:w-1) = swCenter{j}(:,2:w);
                        swCenter{j}(:,w)=r;
                    end
                    obs = [ obs; center 2*ones(2,1) fID*ones(2,1) i*ones(2,1) ];
                end % end if ang > 
            end % end if isempty
    end % end all feature
    
    if i >= step/6 && i < step/6 + 1
        disp('_________ 10% data pocessed.');
    elseif i >= step/3 && i < step/3 + 1
        disp('_________ 30% data pocessed.');
    elseif i >= 4*step/6 && i < 4*step/6 + 1
        disp('_________ 50% data pocessed.');
    elseif i >= 5*step/6 && i < 5*step/6 + 1
        disp('_________ 80% data pocessed.');
    end
    
end

[Xstate,Zstate,invCovMatrix,feaOccurredID] = ...
    InitializeFSProb(obs,ptObs,scan.odom,invCovFea);

% store all previous variables without line.
fs.Xstate = Xstate;
fs.Zstate = Zstate;
fs.invCovMatrix = invCovMatrix;
fs.feaOccurredID = feaOccurredID;
fs.gfsN = scanPara.fsN_global;

disp('_________ 100% data pocessed.');
end

