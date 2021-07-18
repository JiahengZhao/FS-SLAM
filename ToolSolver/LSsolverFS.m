function problem = LSsolverFS(problem)

% LSSOLVER solves problem by least square.
% INPUT:
%               PROBLEM
% OUTPUT:
%               XSTATE                Estimated state vector.
%               ERROR_SUM         Error of state vector.
%               COV                     Estimated covariance of xstate
% Author: Jiaheng Zhao

if isfield(problem.option, 'solver')
    method = problem.option.solver;
end

global ws noise scan
noise = 0.05; % Observation error. unit: m. 
if isfield(problem.option, 'weight')
    ws.odom = reshape(problem.option.weight.odom,3,[]);
    ws.center =problem.option.weight.center;
    ws.fs = problem.option.weight.fs;
else
    ws.odom = [1 1 10];
    ws.center = 1;
    ws.fs = 1;
end

if isfield(problem.option, 'display')
    isdisplay = problem.option.display;
else
    isdisplay = true;
end

if isfield(problem,'fs')
    feaOccurredID = problem.fs.feaOccurredID;
    Xstate = problem.fs.Xstate;
    Zstate = problem.fs.Zstate;
%     if isfield(problem.option,'isfixcov')
%         isfixcov = problem.option.isfixcov;
%     else
%         isfixcov = true;
%     end
end

% Default Argument input
iter_max = 500; % max iteration
error_sum = 9999999999999; % initial error
iteration = [];

switch method
    case 'lmsba' % this L-M method uses sba.
        iter = 0;
        Factor=2;
        t = 1e-3;
        e1 = 1e-10;
        e2 = 1e-10;
        e3 = 1e-10;
        e4 = 0;
        Stop = 0;
        
        %% post-count
        % all points, ellipse and line.
        % pre check
        [fval_odom, fval_center, fval_fs, refsID] = ...
            FuncfFS(feaOccurredID,Xstate,Zstate);
        
        [Jac, covFS] = FunJacFS(feaOccurredID,Xstate,Zstate, refsID);
        
        [error_sum,~] = FunSolveFS(fval_odom, fval_center, fval_fs, covFS);
        
        if isdisplay
            fprintf('Iteration: %d. Initial Error is %.8f\n', iter, error_sum);
        end
        error_pre = error_sum;
        
        fval = [fval_odom; fval_center; fval_fs];
        % add weight
        sigma = 1;
        %-  --%
        iteration = 0.5*sqrt(fval'*sigma*fval);
        A = Jac'*sigma*Jac;
        Lambda = t*max(diag(A));
        G = Jac'*sigma*fval;
        g = max(abs(G));
        if g<=e1
            Stop = 1;
            if isdisplay
                fprintf('Iteration: %d. Reason: 1\n',iter);
            end
        end
        while iter <= iter_max && Stop ~= 1
            iter = iter + 1;
            P = -1;
            innerIter=1;
            while Stop~=1 && P<=0 && innerIter <= 50 % this part is to find a reasonable region.
                if innerIter == 50
                    Stop = 1;
                    if isdisplay
                        fprintf('Inner Iteration Max!!!\n');
                    end
                end
                innerIter = innerIter + 1;
                
                [delt,delt_sum,~,~] = funGetDeltLMFS(Jac,fval_odom, fval_center, fval_fs,Lambda, covFS);
                outX = funGetXFS(Xstate,Zstate);
                
                if delt_sum <= e2*(outX + e2) % If delt is too small, quit
                    Stop = 1;
                    if isdisplay
                        fprintf('Iteration: %d. Reason: 2. Delt_sum: %d.\n',iter,delt_sum);
                    end
                else
                    
                    Xstate= FunUpdateFS(Xstate,delt); % update to Xstate

                    [fval_odom, fval_center, fval_fs, refsID] = ...
                        FuncfFS(feaOccurredID,Xstate,Zstate); % solve fval
                    fval = [fval_odom; fval_center; fval_fs];
                    % add weight
                    [~, covFS] = FunJacFS(feaOccurredID,Xstate,Zstate, refsID);
                    [error_sum,~] = FunSolveFS(fval_odom, fval_center, fval_fs, covFS);
                    
                    P = (error_pre - error_sum)/(delt' * (Lambda * delt - G)); %  Liang here use + G. But mine is - G
                    if P > 0
                        if sqrt(error_pre) - sqrt(error_sum) < e4*sqrt(error_pre)
                            Stop = 1;
                            if isdisplay
                                fprintf('Iteration: %d. Reason: 3\n',iter);
                            end
                        end
%                         fval = [fval_odom; fval_center; fval_fs];
                        Jac = FunJacFS(feaOccurredID,Xstate,Zstate, refsID);
                        % add weight--%
                        G = Jac'*sigma*fval;
                        g = max(abs(G));
                        if Stop ==1 || g<=e1
                            Stop = 1;
                            if isdisplay
                                fprintf('Iteration: %d. Reason: 1\n',iter);
                            end
                        end
                        Lambda = Lambda * max(1/3, 1-(2*P-1)^3);
                        Factor = 2;
                    else % error increases...
                        delt = -delt;
                        Xstate= FunUpdateFS(Xstate,delt); % update to Xstate
                        Lambda = Factor * Lambda;
                        Factor = Factor * 2;
                    end % end if P
                end
            end % end inner while
            if sqrt(error_sum) <= e3
                Stop = 1;
                if isdisplay
                    fprintf('Iteration: %d. Reason: 4\n',iter);
                end
            end
            if P>0
                if isdisplay
                    fprintf('Iteration: %d. Error is %.8f\n', iter, error_sum);
                end
                iteration = [iteration; 0.5*sqrt(fval'*fval)];
                error_pre = error_sum;
            end
        end
        
                [~,~,Info] = funGetDeltLMFS(Jac,fval_odom, fval_center, fval_fs,Lambda, covFS);

        
end
cov = Info \ eye(length(Info));
% cov=[];

problem.solution.Xstate = Xstate;
problem.solution.error_sum = error_sum;
problem.solution.cov = cov;
problem.solution.iteration = iteration;
end

