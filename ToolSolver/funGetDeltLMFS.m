function [delt,delt_sum,Info,solveInfo] = funGetDeltLMFS(Jac,fval_odom, fval_center, fval_fs,Lambda, CovFS)
% This function get delt for post_count method for Levenberg-Marquart
global ws;
sigma = sparse(1:(length(fval_odom)+length(fval_center)+length(fval_fs)),...
    1:(length(fval_odom)+length(fval_center)+length(fval_fs)),...
    [repmat(ws.odom, length(fval_odom)/3,1);...
    ws.center .* ones(length(fval_center),1);...
    ws.fs .*ones(length(fval_fs),1)]);
tmp = cat(1,CovFS{:}); tmp = 1./tmp;
sigma = sparse(1:(length(fval_odom)+length(fval_center)+length(fval_fs)),...
    1:(length(fval_odom)+length(fval_center)+length(fval_fs)),...
    [repmat(ws.odom, length(fval_odom)/3,1);...
    ws.center .* ones(length(fval_center),1);...
    tmp]);
fval = [fval_odom;fval_center;fval_fs];

[~,nColumNum] = size(Jac);

C = sparse(1:nColumNum,1:nColumNum,diag(Jac' * Jac))*Lambda; % general selection

%%
solveInfo = Jac' * sigma * Jac;
Info = solveInfo + C;
E = Jac' * sigma * fval;
delt = -Info\E;
if ~isreal(delt)
    delt;
end
% delt = -pinv(Info)*E;
delt_sum = sqrt(delt'*delt);
clear Jac fval Lambda C
end