% remove search folder.

matlabdir = fileparts(mfilename('fullpath'));
rmpath(fullfile(matlabdir, '.'))
rmpath(genpath(fullfile(matlabdir, filesep, 'Data')))
rmpath(genpath(fullfile(matlabdir, filesep, 'DataProcessing')))
rmpath(genpath(fullfile(matlabdir, filesep, 'SensorModel')))
rmpath(genpath(fullfile(matlabdir, filesep, 'ToolPlot')))
rmpath(genpath(fullfile(matlabdir, filesep, 'ToolSolver')))
rmpath(genpath(fullfile(matlabdir, filesep, 'ToolTransform')))
rmpath(genpath(fullfile(matlabdir, filesep, 'MapJoining')));
clear 'matlabdir'
