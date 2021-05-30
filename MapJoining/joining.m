function gmap = joining(fs,lmap)
% JOINING join sub maps together.
% Author: Jiaheng Zhao

% Initialize parameter. all feature are in the initial frame.
[Xstate,Zstate,lfeaID] = GmapIni(lmap,fs.gfsN);

fs.Xstate = Xstate; fs.Zstate = Zstate; fs.invCovMatrix = 1; fs.feaOccurredID = lfeaID;

problem.noLine = []; problem.Ellipse = []; problem.Line = []; problem.AnyFea = [];
problem.fs = fs; problem.option.solver = 'lmsba';  problem.option.display = true; 
problem.option.weight.odom = [10,10,10];
problem.option.weight.center = 1;
problem.option.weight.fs = 150;
gmap = LSsolverFS(problem);

end