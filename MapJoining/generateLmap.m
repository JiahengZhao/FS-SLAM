function lmap = generateLmap(fs, nLmap)
% GENERATELMAP generate all local maps.
% Author: Jiaheng Zhao

nstep = max(fs.Xstate(fs.Xstate(:,2)==1,3));
nmap = ceil(nstep/nLmap);
lmap = struct([]);
lmap(1).id = 1;
lmap(1).Xstate = [];
lmap(1).Zstate.center =[]; lmap(1).Zstate.pts =[]; lmap(1).Zstate.anbn =[];
lmap(1).localFeaID = [];
lmap(1).fsN = [];
lmap(1).endPose = [];
% build local map.
for ni = 1: nmap
    obs = fs.Zstate.center((fs.Zstate.center(:,4)<= (ni*nLmap)) &  (fs.Zstate.center(:,4)> ((ni-1)*nLmap)),:);
    if ni == nmap
        ptObs = fs.Zstate.pts(((ni-1)*nLmap+1) : end,:);
    else
        ptObs = fs.Zstate.pts(((ni-1)*nLmap+1) : ni*nLmap,:);
    end
    control = fs.Zstate.odom((fs.Zstate.odom(:,3)<=(ni*nLmap)) & (fs.Zstate.odom(:,3)>((ni-1)*nLmap)),:);
    lmap(ni) = lmapGnr(obs,ptObs,control,fs.Zstate.fsN,fs.feaOccurredID);
    lmap(ni).id = ni;
    lmap(ni).Zstate.center(:,4) = ni;
    lmap(ni).Zstate.anbn(:,4) = ni;
    disp(['Generate local map: ', num2str(ni), '/',num2str(nmap)]);
end
 disp('                   local map generation finished. ');
end
