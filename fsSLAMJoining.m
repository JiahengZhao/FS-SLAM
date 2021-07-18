% Test Fourier Series SLAM map joining. Assume data association is finished.
% Author: Jiaheng Zhao

clear all;

% load data
finaldata = fullfile(pwd,filesep,'Data',filesep,'ICRA', filesep,'demo_techlab.mat'); %
scan = load(finaldata);
scan = filtscan(scan,5); % filter some data

global scanPara
% load lidar parameters.
scanPara = setLidarParametersMapJoining(); 

t1 = tic;
% convert to my state
fs = ToStateFS(scan,scanPara);

% each local map includes (nLmap+1) poses.
nLmap = 5; 

% start submap joining
[gmap,lmap]  = mapjoining(fs, nLmap);

t2 = toc(t1);
disp(['Total runtime: ', num2str(t2),' sec.']);

%% Plot
Xstate = gmap.solution.Xstate;
PlotXstatewithOdom = Xstate(Xstate(:,2)==1,1);
PlotXstatewithFS = Xstate(Xstate(:,2)==2,1);

f1 = figure;
gt = scan.groundtruth; 
colorFS = [0.8500 0.3250 0.0980];
colorIm = [0 0.4470 0.7410];
hold on;

gtp = plot(scan.groundtruth(1,:),scan.groundtruth(2,:),'k','LineWidth',1);
fsgt = plotFS(scan.feature(1:2,:),scan.feature(3:end,:),'k');
FS = reshape(PlotXstatewithFS(:,1),2*fs.gfsN+3,[]);
fs1 = plotFS(FS(1:2,:),FS(3:end,:),colorFS); % plot feature
sbmappose = reshape(PlotXstatewithOdom,3,[]);
pfs1 = plot(sbmappose(1,:),sbmappose(2,:),'Color',colorFS,'LineStyle',':','LineWidth',1.5,...
    'Marker','o','MarkerFaceColor',colorFS,'MarkerEdgeColor','k','MarkerSize',5);
pfs2 = plot([sbmappose(1,:); sbmappose(1,:)+ 0.5.*cos(sbmappose(3,:))],...
    [sbmappose(2,:); sbmappose(2,:)+ 0.5.*sin(sbmappose(3,:))],'Color',colorFS,'LineStyle','-','LineWidth',2);
axis equal
backcolor = [ 0.9020    0.9020    0.9020];
box on;
title('Submap Joining');
set(gca,'FontSize',36);
xlabel('x/m');
ylabel('y/m')


function lidar = setLidarParametersMapJoining()
lidar.covscan = ([1 ; 1].*1e-2).^2;
lidar.covodom =  diag([0.002 0.002 6.2500e-04]); 
lidar.angle_min = -1.9199064;
lidar.angle_max =  1.9252524;
lidar.angle_increment = 0.005816;
lidar.npoints   = 662;
lidar.range_min = 0.05;
lidar.range_max = 25;
lidar.scan_time = 0;
lidar.time_increment  = 0;
angles = lidar.angle_min + [ 0 :lidar.npoints-1]' .* lidar.angle_increment ;
lidar.angles = wrapToPi(angles);
lidar.fsN_local = 7; % recommend < 7
lidar.fsN_global = 3; % recommend < 5
lidar.fsN_rect = 7; %useless
lidar.fsN_border = 15; %useless
end
