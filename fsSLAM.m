% Fourier Series SLAM. Assume data association is done.
% Author: Jiaheng Zhao

clear all;

global scan scanPara

experiment = 'demo_carpark';  % demo_simu - simulation | demo_techlab - techlab data | demo_carpark

switch experiment
    case 'demo_simu'
        finaldata = fullfile(pwd,filesep,'Data',filesep,'ICRA',filesep,'demo_simu.mat'); % This is the longer one
        scan = load(finaldata);
        scanPara = setLidarParameters(); % Set lidar parameters
        isRect = true;
        isGoodInitial = false;
    case 'demo_techlab'
        finaldata = fullfile(pwd,filesep,'Data',filesep,'ICRA', filesep,'demo_techlab.mat'); %
        scan = load(finaldata);
        scan = filtscan(scan,7); % filt some data
        scanPara = setLidarParameters(); % Set lidar parameters
        scan.borderF = [];
        scan.rectF = [];
        isRect = false;
        isGoodInitial = false;
    case 'demo_carpark'
        finaldata = fullfile(pwd,filesep,'Data',filesep,'ICRA', filesep,'carPark.mat'); %
        scan = load(finaldata);
        scan = filtscan(scan,7); % filt some data
        scanPara = setLidarParameters2(); % this is for the robot in car park
        isRect = true;
        isGoodInitial = true; % sometimes we still need to find a good initial guess.
end


t0 = tic;
fs = ToStateFS(scan,scanPara);

problem.fs = fs;

%initial guess
if isGoodInitial
    goodInitFile = fullfile(pwd,filesep,'Data',filesep,'ICRA',filesep,['goodInitFile_',experiment,'.mat']); %
    if isfile(goodInitFile)
        gIXstate = load(goodInitFile);
        problem.fs.Xstate = gIXstate.tsXstate;
    end
else
    goodInitFile = '';
end

%%
problem.option.solver = 'lmsba';
problem.option.weight.odom = [2500,2500,5e6];
problem.option.weight.center = 1;
problem.option.weight.fs = 10; % this parameter is not used.
% problem.option.weight.odom = [50,50,500];
% problem.option.weight.center = 1;
% problem.option.weight.fs = 100;

t1 = tic;
problem = LSsolverFS(problem); % This is a solver with LSBA
t2 = toc(t1);
t3 = toc(t0);
disp(['Total runtime: ', num2str(t3),' sec.']);

XstateFS = problem.solution.Xstate;

% ready to plot?
PlotXstatewithOdom1 = XstateFS(XstateFS(:,2)==1,:);
PlotXstatewithFS1 = XstateFS(XstateFS(:,2)==2,:);
ec = (PlotXstatewithOdom1(1:2,1)-scan.groundtruth(1:2,1));

if isRect
    PlotXstatewithFS1 = XstateFS(XstateFS(:,2)==2 & XstateFS(:,3) <size(scan.feature,2)+1 ,:);
    PlotXstatewithFS2 = XstateFS(XstateFS(:,2)==2 & XstateFS(:,3) >size(scan.feature,2) & XstateFS(:,3) <size(scan.feature,2)+size(scan.borderF,2)+1 ,:);
    PlotXstatewithFS3 = XstateFS(XstateFS(:,2)==2 & XstateFS(:,3) >size(scan.feature,2)+size(scan.borderF,2) ,:);
    tmp1 = reshape(PlotXstatewithOdom1(:,1),3,[]);tmp1(1:2,:) = tmp1(1:2,:) - ec;
    PlotXstatewithOdom1(:,1) = tmp1(:);
else
    ec = (PlotXstatewithOdom1(1:2,1)-scan.groundtruth(1:2,1));
    tmp1 = reshape(PlotXstatewithOdom1(:,1),3,[]);tmp1(1:2,:) = tmp1(1:2,:) - ec;
    PlotXstatewithOdom1(:,1) = tmp1(:);
    tmp2 =  reshape(PlotXstatewithFS1(:,1),2*fs.Zstate.fsN+3,[]);tmp2(1:2,:) = tmp2(1:2,:) -ec;
    PlotXstatewithFS1(:,1) = tmp2(:);
end

% Go to plot
f1 = figure;
titleStr = strcat('\fontsize{24}','Solver: ',problem.option.solver, '.');
gt = scan.groundtruth;

colorFS = [0.8500 0.3250 0.0980];
colorIm = [0 0.4470 0.7410];

increFig=increPloter(PlotXstatewithOdom1,gt,[],[],colorFS,4,colorFS,2); % plot feature
delete(increFig(2).figGt);delete(increFig(1).figGt); delete(increFig(1).figPs);delete(increFig(2).figPs);
axis equal
backcolor = [ 0.9020    0.9020    0.9020];
box on;

% test an bn plot
fsgt = plotFS(scan.feature(1:2,:),scan.feature(3:end,:),'k');
if ~isRect
    FS1 = reshape(PlotXstatewithFS1(:,1),2*fs.Zstate.fsN+3,[]);
    fs1 = plotFS(FS1(1:2,:),FS1(3:end,:),colorFS);
else
    for s3 = 1:size(scan.borderF,2)
        tmp = scan.borderF(:,s3);
        plot ([tmp(1) tmp(1)+tmp(3) tmp(1) tmp(1) tmp(1)+tmp(3) tmp(1)+tmp(3)  tmp(1)+tmp(3)  tmp(1)]', ...
            [tmp(2) tmp(2) tmp(2) tmp(2)+tmp(4) tmp(2)+tmp(4) tmp(2)  tmp(2)+tmp(4)  tmp(2)+tmp(4)]','k-', 'LineWidth', 3);
    end
    for s2 = 1:size(scan.rectF,2)
        tmp = scan.rectF(:,s2);
        plot ([tmp(1) tmp(1)+tmp(3) tmp(1) tmp(1) tmp(1)+tmp(3) tmp(1)+tmp(3)  tmp(1)+tmp(3)  tmp(1)]', ...
            [tmp(2) tmp(2) tmp(2) tmp(2)+tmp(4) tmp(2)+tmp(4) tmp(2)  tmp(2)+tmp(4)  tmp(2)+tmp(4)]','k-', 'LineWidth', 3);
    end
    FS1 = reshape(PlotXstatewithFS1(:,1),2*fs.Zstate.fsN+3,[]);
    FS1(1:2,:) = FS1(1:2,:) - ec;
    fs1 = plotFS(FS1(1:2,:),FS1(3:end,:),colorFS);
    FS2= reshape(PlotXstatewithFS2(:,1),2*fs.Zstate.fsN_border+3,[]);
    FS2(1:2,:) = FS2(1:2,:) - ec;
    fs2 = plotFS(FS2(1:2,:),FS2(3:end,:),colorFS);
    FS3 = reshape(PlotXstatewithFS3(:,1),2*fs.Zstate.fsN_rect+3,[]);
    FS3(1:2,:) = FS3(1:2,:) - ec;
    fs3 = plotFS(FS3(1:2,:),FS3(3:end,:),colorFS);
end
Pose_state1 = PlotXstatewithOdom1(PlotXstatewithOdom1(:,2)==1,:);
step = max(Pose_state1(:,3));
tmpZstate = cell(step,size(scan.feature,2)+size(scan.rectF,2)+size(scan.borderF,2));
for i = 1:step
    ri = i;
    PoseState_ThisStep = Pose_state1(3*ri-2:3*ri,1);
    phi = PoseState_ThisStep(3);
    R = theta2R(phi);
    transFun = @(x) R*x + PoseState_ThisStep(1:2);
    tmpZstate(i,:) = cellfun(transFun,fs.Zstate.pts(ri,:),'UniformOutput',false);
    testtmpZstate{1,i} = transFun(scan.scan{1,ri}(2:3,:));
end

% plot points
% testAllPoints = cat(2,testtmpZstate{:});
% testAllPoints_fea = cat(2,tmpZstate{:});
title('Trajectory comparison')
set(gca,'FontSize',38);

%%
% Below is to adjust the initial guess and save as new initial value.
if ~isfile(goodInitFile) && isGoodInitial
    figure;
    hold on;
    axis equal;
    tsXstate = XstateFS;
    thisfeept = cell(1, size(fs.feaOccurredID,1));
    for i = 1 : size(fs.feaOccurredID,1)
        oid = fs.feaOccurredID(i,2);
        if oid > 9 && oid <100
            fsN = scanPara.fsN_border;
            comp = true;
        elseif oid > 99
            fsN = scanPara.fsN_rect;
        else
            fsN = scanPara.fsN_local;
        end
        thisfeept{1,i} = cat(2,tmpZstate{:,i});
        [V, center, ~] = fitWithFS_update( fsN, thisfeept{1,i}(1,:),thisfeept{1,i}(2,:), [],true); % isAdd sometimes cause trouble.
        tsXstate(tsXstate(:,2)==2 & tsXstate(:,3)==i,1) = [center;V];
        p1 = plot(thisfeept{1,i}(1,:),thisfeept{1,i}(2,:),'b.');
        p2 = plotFS(center,V,'g');
    end
    save(goodInitFile,'tsXstate');
end

function lidar = setLidarParameters2()
% SETLIDARPARAMETERS Set Lidar parameters.
% Lidar Parameters on the robot.
lidar.covscan = ([1 ; 1].*1e-2).^2;
lidar.covodom =  diag([0.002 0.002 6.2500e-04]); %diag([1 1 1]);%
lidar.angle_min = -4.71238899231;
lidar.angle_max =  pi/2;
lidar.angle_increment = 0.010471976;
lidar.npoints   = 600;
lidar.range_min = 1;
lidar.range_max = 50;
lidar.scan_time = 0;
lidar.time_increment  = 0;
angles = lidar.angle_min + [ 0 :lidar.npoints-1]' .* lidar.angle_increment ;
lidar.angles = wrapToPi(angles);
lidar.fsN_local = 4;
lidar.fsN_global = 9;
lidar.fsN_rect = 15;
lidar.fsN_border = 33;
end
