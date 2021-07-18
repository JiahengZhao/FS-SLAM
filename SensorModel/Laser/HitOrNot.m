function [fval,info,range] = HitOrNot(f,border,ang,x0,y0,maxrange,AnyFea)
% find hitted point
% AnyFea: m by n, m is the number of features, n is the number of parameters.
if nargin<7
    AnyFea = []; 
end
fval = [];
info = [];
numFea = size(f,2);
% numBorder = size(border,2);

F = [f border]; % new feature including border.
nF = size(F,2);
for idf = 1:nF
    [tmpf,tmpinfo] = intersection(F(:,idf),ang,x0,y0); % Hit or not with specific feature
    if idf > numFea
        tmpinfo = tmpinfo * 2;
    end
    fval = [fval tmpf];
    info = [info; tmpinfo]; 
end
if ~isempty(AnyFea)
    
    for ida = 1:size(AnyFea,1)
        [tmpfa,tmpinfoa] = intersectionAnyFea(AnyFea(ida,:),ang,x0,y0);
    end
    if ~isempty(tmpfa)
        fval = [fval tmpfa];
        info = [info; tmpinfoa];
    end
    
end

[a,id]=max(info); % judge how many intersections are.
if a>0
    if size(fval,2)~=1 % if more than one intersection, choose the nearer one.
        tmp_id = find(info>0); % find intersection of feature and border
        dist=vecnorm( fval-[x0;y0],2,1);
        [~,id_2]=min(dist); % find the nearest point
        if min(dist) > maxrange
            info = 0;
            fval = nan(3,1);
            range = 60;
            return
        end
        %----------For RAL Simulation, assign line segments into the same id--------%
        % line 5 = line 6, line 7 = line 8
        if tmp_id(id_2) > numFea && tmp_id(id_2) <= nF% border id
            if (tmp_id(id_2)-numFea) == 6
                flagid = 50;
            elseif (tmp_id(id_2)-numFea) == 7 || (tmp_id(id_2)-numFea) == 8
                flagid = 60;
            else
                flagid = 10 * (tmp_id(id_2)-numFea); %  Set 0 as suffix;
            end
        elseif tmp_id(id_2) > nF
            flagid = 100*(tmp_id(id_2)-nF);
        else
            flagid = tmp_id(id_2);
        end
             fval = [flagid; fval(:,id_2)];
             info = 1;
             range = dist(id_2);
    else
        dist=norm( fval-[x0;y0]);
        if dist > maxrange
            info = 0;
            fval = nan(3,1);
            range = 60;
            return
        end
        if id > numFea && id <= nF % boader
            id = 10 * (id-numFea);
        elseif id > nF
            id = 100 * (id - nF);
        end
        fval = [id ; fval];
        info = 1;
        range = dist;
    end
else
    fval = nan(3,1);
    info = 0;
    range = 60;
end
end