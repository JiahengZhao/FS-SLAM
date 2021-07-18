function [Jac,CovFS] = FunJacFS(feaOccurredID,Xstate,Zstate, refsID)
%FUNJAC solve Jacobian of f(x).
%       Jac is odometry + center + fs
% Author: Jiaheng Zhao
if nargin < 4
    refsID = [];
end
Pose_state = Xstate(Xstate(:,2)==1,:);
step = max(Pose_state(:,3));
% global ws;
% using sparse matrix. Odompart should be with the same size of Featurepart
id1Pose = [];
id2Pose = [];
valPose = []; 
global noise
fsN_f = Zstate.fsN;
if isfield(Zstate, 'fsN_rect')
    fsN_rect = Zstate.fsN_rect;
else
    fsN_rect = 0;
end
if isfield(Zstate, 'fsN_border')
    fsN_border = Zstate.fsN_border;
else
    fsN_border = 0;
end
nF = size(find(feaOccurredID(:,2)<10),1);
nB = size(find(feaOccurredID(:,2)>9 & feaOccurredID(:,2) < 100),1);
nR = size(find(feaOccurredID(:,2)>99),1);

ZstateCenter = Zstate.center;
ZstatePts = Zstate.pts;
ZstateOdom = Zstate.odom;
id1Center = [];
id2Center = [];
valCenter = [];
numAllPt = cellfun('size',ZstatePts,2); % in old id order, need change to new id order
cid1 = feaOccurredID(:,2) > 9 & feaOccurredID(:,2) < 100;
cid2 = feaOccurredID(:,2) > 99;
tmp_feaOccurredID = feaOccurredID;
tmp_feaOccurredID(cid1,2) = feaOccurredID(cid1,2)/10+nF;
tmp_feaOccurredID(cid2,2) = feaOccurredID(cid2,2)/100 + nF + nB;
numAllPt = numAllPt(:, tmp_feaOccurredID(:,2)'); % Change to new id order
numAllPt_save = zeros(size(numAllPt));

nFea = size(feaOccurredID,1);
id1fs = cell(1,nFea);
id2fs = id1fs;
valfs = id1fs;
CovFS = id1fs; % covariance of FS

nallSeeFea = 0;
for i=1:step % solve value by pose number.
    %% ----------------------------------------       Pose      -----------------------------------------------%%
    PoseState_ThisStep = Pose_state(3*i-2:3*i,1);
    phi = PoseState_ThisStep(3);
    R = theta2R(phi);
    dR = difftheta2R(phi);
    % Odometry part of Jacobbian
    if i == 1
        id1Pose = [1;2;3];
        id2Pose = [1;2;3];
        valPose = [-1;-1;-1];
    else % here should array sparse(i,j,v).
        PoseState_LastStep = Pose_state(3*i-5:3*i-3,1);
        phi_last = PoseState_LastStep(3);
        R_last = theta2R(phi_last);
        dR_last = difftheta2R(phi_last);
        
        id1PoseTmp = 3*i-2:3*i;
        id1PoseTmp = repmat(id1PoseTmp',6,1);
        id1Pose =  [id1Pose ; id1PoseTmp];
        
        id2PoseTmp = 3*i-5:3*i;
        id2PoseTmp = repmat(id2PoseTmp,3,1);
        id2PoseTmp = id2PoseTmp(:);
        id2Pose = [id2Pose;id2PoseTmp];
        tmpJPosepart = [ R_last', -dR_last'*(PoseState_ThisStep(1:2)-PoseState_LastStep(1:2)), -R_last',zeros(2,1);
            0 0 1 0 0 -1];
        tmpJPosepart=tmpJPosepart(:);
        valPose = [valPose;tmpJPosepart];
    end
    %% ---------------------------------  Center Jac  & FS Jac -------------------------------------%%
    % Feature Points part of Jacobbian
    thisStepZstateCenter = ZstateCenter(ZstateCenter(:,4)==i,:);
    thisSeeFea = unique(thisStepZstateCenter(:,3));
    nthisSeeFea = length(thisSeeFea);
    tmpCount = 1;
    for sid=1:nFea % current seen feature
        newsid = sid;
        oldsid = feaOccurredID(feaOccurredID(:,1)==newsid,2);
        [a,~]=ismember(oldsid,thisSeeFea);

        % for erct. The order is always " General feature -> Border feature -> rectangle feature
        if oldsid > 9 && oldsid < 100
            Feeindex_old_1 = oldsid/10;
            Feeindex_old = Feeindex_old_1 + nF;
            fsN = fsN_border;
            xid2_center = (3+2*fsN_f)*nF + (3+2*fsN)*Feeindex_old_1 - (2+2*fsN) : ...
                (3+2*fsN_f)*nF + (3+2*fsN)*Feeindex_old_1 - (1+2*fsN);
            xid2_fs = (3+2*fsN_f)*nF + (3+2*fsN)*Feeindex_old_1 - (2+2*fsN) : ...
                (3+2*fsN_f)*nF + (3+2*fsN)*Feeindex_old_1;
        elseif oldsid > 99
            fsN = fsN_rect;
            Feeindex_old_1 = oldsid/100;
            Feeindex_old = Feeindex_old_1 + nF + nB;
            xid2_center = (3+2*fsN_f)*nF + (3+2*fsN_border)*nB + (3+2*fsN)*Feeindex_old_1 - (2+2*fsN) : ...
                (3+2*fsN_f)*nF + (3+2*fsN_border)*nB + (3+2*fsN)*Feeindex_old_1 - (1+2*fsN);
            xid2_fs = (3+2*fsN_f)*nF + (3+2*fsN_border)*nB + (3+2*fsN)*Feeindex_old_1 - (2+2*fsN) : ...
                (3+2*fsN_f)*nF + (3+2*fsN_border)*nB + (3+2*fsN)*Feeindex_old_1;
        else
            fsN = fsN_f;
            Feeindex_old = oldsid;
            xid2_center = (3+2*fsN)*sid - (2+2*fsN) : ...
                                 (3+2*fsN)*sid - (1+2*fsN);
            xid2_fs = (3+2*fsN)*sid - (2+2*fsN) : ...
                                 (3+2*fsN)*sid;
        end
        
        if a
            if tmpCount > nthisSeeFea
                error('tmpCount is larger than the number of seen feature');
            end
            % below is center jac
            id1Center_tmp = 3*step+2*nallSeeFea+2*tmpCount-1:3*step+2*nallSeeFea+2*tmpCount; % this is related to the number of cost function.
            id1Center_tmp = repmat(id1Center_tmp',5,1  );
            id1Center = [id1Center; id1Center_tmp(:)];
            
            id2Center_tmp = [3*i-2:3*i, 3*step+xid2_center];
%             id2Center_tmp = [3*i-2:3*i 3*step+(2*fsN+3)*newsid-(2*fsN+2):3*step+(2*fsN+3)*newsid-(2*fsN+1)];
            id2Center_tmp = repmat(id2Center_tmp,2,1);
            id2Center = [id2Center; id2Center_tmp(:)];
            jC = thisStepZstateCenter(thisStepZstateCenter(:,3)==oldsid,1);
            valCenter_tmp = [eye(2,2) dR*jC -eye(2,2)];
            valCenter = [valCenter; valCenter_tmp(:)];
            tmpCount = tmpCount+1;
        end
        % below is fs jac
        thiFeaXstate = Xstate((Xstate(:,2)==2) & (Xstate(:,3) == newsid),1);
        an = thiFeaXstate(3:fsN+3); bn =  thiFeaXstate(fsN+4:end); bn = [0;bn];
        numAllPt_save(i,newsid) = numAllPt(i,newsid); % sometimes one scan doesn't have fitted center
        thisNumPt = numAllPt_save(i,newsid); %here is new id order
        if i==1
            thisStartid1=1;
        else
            thisStartid1 = sum( numAllPt_save(1:i-1,newsid))+1; % this is the start id at ith step for sth ellipse
        end
        thisRangeid2 = [3*i-2:3*i 3*step+xid2_fs]';
%         thisRangeid2 = [3*i-2:3*i 3*step+(2*fsN+3)*newsid-2*fsN-2:3*step+(2*fsN+3)*newsid]';
        
        jhatxy = R * ZstatePts{i,Feeindex_old};
        hatxy = jhatxy + PoseState_ThisStep(1:2) - thiFeaXstate(1:2);
        thetak = atan2(hatxy(2,:),hatxy(1,:)); %thetak=thetak';
        rk = sqrt(hatxy(1,:).^2+hatxy(2,:).^2);  % rk = rk';
        tmpxy1 = 2.*hatxy(2,:).^2./(hatxy(1,:).*rk.^2);
        tmpxy1x_test = hatxy(2,:)./(rk.^2);
        tmpxy1y_test = hatxy(1,:)./(rk.^2);


        drdtjx = hatxy(1,:)./rk; drdtjy = hatxy(2,:)./rk;
        
        tmpxy2 = 2 .* hatxy(2,:)./(hatxy(1,:)+hatxy(2,:).^2./hatxy(1,:));
        tmpxy2_test = 1./(1+(hatxy(2,:)./hatxy(1,:)).^2);
         
        tmpxy3 = jhatxy(1,:)./hatxy(1,:) + jhatxy(2,:).*hatxy(2,:)./(hatxy(1,:).^2);
        drdphi = (jhatxy(1,:).* hatxy(2,:) - jhatxy(2,:).*hatxy(1,:))./rk;
        
        dBdtx = zeros(size(thetak)); dBdty= dBdtx; dBdphi = dBdtx; dBdCx = dBdtx; dBdCy = dBdtx; 
        dBdpx = dBdtx; dBdpy = dBdtx;
        dBdan = zeros(fsN+1, length(thetak)); dBdbn = dBdan; % later we need to transpose it.
        for nfs = 0:fsN
            ns = nfs.*sin(nfs.*thetak);
            nc = nfs.*cos(nfs.*thetak);
            dBdan(nfs+1,:) = cos(nfs.*thetak);
            dBdbn(nfs+1,:) = sin(nfs.*thetak);
            
            dcosdtjx = ns.*tmpxy1x_test;             
            dsindtjx =-nc.*tmpxy1x_test;
            dcosdtjy = - ns.*tmpxy1y_test;
            dsindtjy =nc.*tmpxy1y_test;            
            
            dBdtx = dBdtx + an(nfs+1).*dcosdtjx + bn(nfs+1).*dsindtjx;
            dBdty = dBdty +  an(nfs+1).*dcosdtjy + bn(nfs+1).*dsindtjy;
            
            dcosdphi = -ns .* tmpxy2_test .* tmpxy3;
            dsindphi = nc .* tmpxy2_test .* tmpxy3;
            dBdphi = dBdphi +  an(nfs+1).*dcosdphi + bn(nfs+1).*dsindphi;
            
            %--- covariance --%
            dBdpx = dBdpx + an(nfs+1).*dcosdtjx + bn(nfs+1).*dsindtjx;
            dBdpy = dBdpy + an(nfs+1).*dcosdtjy + bn(nfs+1).*dsindtjy;
        end
        dBdCx = -dBdtx + drdtjx;    dBdCx = dBdCx';
        dBdCy = -dBdty + drdtjy;      dBdCy = dBdCy';
        dBdtx = dBdtx - drdtjx;        dBdtx = dBdtx';
        dBdty = dBdty - drdtjy;        dBdty = dBdty';
        dBdphi = dBdphi - drdphi;   dBdphi = dBdphi';
        dBdbn(1,:) = [];    dBdbn = dBdbn';
        dBdan = dBdan';

        tmpJacthis =  [dBdtx dBdty dBdphi dBdCx dBdCy dBdan dBdbn] ;
        tmpJacthis = tmpJacthis';
        
        tmpreid1 = refsID{i,Feeindex_old}; tmpreid1 = tmpreid1';
        
        thisid1 = repmat(tmpreid1,(2*fsN+6),1);
        thisid2 =  repmat(thisRangeid2 ,1,thisNumPt);
        thisval = tmpJacthis(:);

        id1fs{1,newsid} = [ id1fs{1,newsid}; thisid1(:) ];
        id2fs{1,newsid} = [ id2fs{1,newsid}; thisid2(:) ];
        valfs{1,newsid} = [ valfs{1,newsid}; thisval ];
        
        % covariance
        Jdp = [ dBdtx dBdty ] * R; % This is the Jacobian of p
        sig_z = noise.^2.* eye(2,2);
        CovFS{1,sid} = [CovFS{1,sid}; sum(Jdp * sig_z .* Jdp,2)];  % lemma 1
    end
    nallSeeFea = nallSeeFea + nthisSeeFea;
    
end % end step
%% All Jac
startIDAllpt_fs = sum(numAllPt_save,1); %
numOdomCenterRow = max(id1Center);
id1fsFinal = [];
for i=1:nFea
    if i==1
        id1fsFinal = [id1fsFinal; id1fs{i} + numOdomCenterRow];
    else
        id1fsFinal = [id1fsFinal; id1fs{i} + sum(startIDAllpt_fs(1:i-1)) + numOdomCenterRow];
    end
end
id2fsFinal = cat(1,id2fs{:});
valfsFinal = cat(1,valfs{:});

id1 = [id1Pose;id1Center; id1fsFinal];
id2 = [id2Pose;id2Center; id2fsFinal];
val = [valPose; valCenter; valfsFinal];
Jac = sparse(id1,id2,val);

end

