% Add the folders of the repository to the MATLAB search path.
matlabdir = fileparts(mfilename('fullpath'));
addpath(fullfile(matlabdir, '.'))
addpath(genpath(fullfile(matlabdir, filesep, 'Data')))
addpath(genpath(fullfile(matlabdir, filesep, 'DataProcessing')))
addpath(genpath(fullfile(matlabdir, filesep, 'SensorModel')))
addpath(genpath(fullfile(matlabdir, filesep, 'ToolPlot')))
addpath(genpath(fullfile(matlabdir, filesep, 'ToolSolver')))
addpath(genpath(fullfile(matlabdir, filesep, 'ToolTransform')))
addpath(genpath(fullfile(matlabdir, filesep, 'MapJoining')))
% outdir = fullfile(matlabdir, 'output');
% if exist(outdir,'dir')
%     addpath(outdir)
% end
% addpath(fullfile(matlabdir, 'script'))

% Clear all created variables.
clear 'matlabdir'
