function [gmap,lmap] = mapjoining(fs, nLmap)
% MAPJOINING map joining method
% Input: fs - processed data
%           nLmap - number of steps at each local map.
% Author: Jiaheng Zhao

% generate all local maps;
disp('Start local map generating...');
lmap = generateLmap(fs, nLmap);

% joining all local maps;
disp('Start map joining...');
t1 = tic;
gmap = joining(fs,lmap);
t2 = toc(t1);
disp(['Optimizing runtime: ', num2str(t2),' sec.']);

end