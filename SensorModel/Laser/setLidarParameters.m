function lidar = setLidarParameters()
% SETLIDARPARAMETERS Set Lidar parameters.

lidar.useAllpoints = 0; % 0 for the original version.
%% Real Fetch Parameters
lidar.covscan = ([1 ; 1].*1e-2).^2;
lidar.covodom =  diag([0.002 0.002 6.2500e-04]); %diag([1 1 1]);%
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

lidar.fsN_local = 4; % 4 for the traj, 7 for the submap joing
lidar.fsN_global = 9;
lidar.fsN_rect = 7;
lidar.fsN_border = 15;