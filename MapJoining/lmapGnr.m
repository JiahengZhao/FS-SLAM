function  lmap = lmapGnr(CenterState,ScanState_pre,control,fsN,gfeaID)
% LMAPGNR generate local map.
% Author: Jiaheng Zhao

% Initialize parameter. all feature are in the initial frame.
[Xstate,Zstate,lfeaID] = lmapIni(CenterState,ScanState_pre,control,fsN);

fs.Xstate = Xstate; fs.Zstate = Zstate; fs.invCovMatrix = 1; fs.feaOccurredID = lfeaID;
problem.noLine = []; problem.Ellipse = []; problem.Line = []; problem.AnyFea = [];
problem.fs = fs; problem.option.solver = 'lmsba';  problem.option.display = false; 
problem.option.weight.odom = [50,50,50];
problem.option.weight.center =1;
problem.option.weight.fs = 10;
problem = LSsolverFS(problem);

%%
% transfer feature to the last frame.
% compare local feature id with global feature id.
problem.solution.Xstate = rearrangeID(problem.solution.Xstate,lfeaID,gfeaID);
% transfer points anbn pts to the last pose frame.
 [pose_end,CenterState,ScanState,anbn] = transFrameFS(problem.solution.Xstate,ScanState_pre,fsN);

lmap.id = [];
lmap.Xstate = problem.solution.Xstate;
lmap.Zstate.center = CenterState; % need to transfer to the end pose
lmap.Zstate.pts = ScanState; % need to transfer to the end pose
lmap.Zstate.anbn = anbn;
lmap.localFeaID = lfeaID;
lmap.fsN = fsN;
lmap.endPose = pose_end;
end

function XstateFea = rearrangeID(Xstate,lfeaID,gfeaID)
    XstateFea = Xstate;
    for i = 1: size(lfeaID,1)
        thisInGfeaID = gfeaID(gfeaID(:,2) == lfeaID(i,2),1); % this id in global fea id
        XstateFea((Xstate(:,3) == lfeaID(i,1)) & (Xstate(:,2)==2),3) = thisInGfeaID;
    end
end