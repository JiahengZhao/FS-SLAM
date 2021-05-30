function [error,invCov] = FunSolveFS(fval_odom, fval_center, fval_fs, CovFS)
% This function solve fval for 3 cases.
global ws;
tmp = cat(1,CovFS{:}); tmp = 1./tmp;
sigma = sparse(1:(length(fval_odom)+length(fval_center)+length(fval_fs)),...
    1:(length(fval_odom)+length(fval_center)+length(fval_fs)),...
    [repmat(ws.odom, length(fval_odom)/3,1);...
    ws.center .* ones(length(fval_center),1);...
    tmp]);
fval = [fval_odom; fval_center; fval_fs];
% num1 = length(fval_odom);
% num2 = length(fval_center);
% num3 = length(fval_fs);
% num0 = size(fval,1);
invCov=[];
error = fval'*sigma* fval;
error = 0.5 * sqrt(error);
end