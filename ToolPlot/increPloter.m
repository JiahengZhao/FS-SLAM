classdef increPloter
    properties
        figLaser,figGt,figPs,figLink,figFea,figGraph,Xstate,figCov;
    end
    methods
        function obj = increPloter(Xstate,Groundtruth,cov,scan,color,linewidth,colorCOV,sdwidth,filterSize)
            %             warning off;
            if nargin ~= 0
                obj(3) = obj;
                obj(1).figGraph = [];
                obj(1).figCov = [];
                step=max(Xstate(Xstate(:,2)==1,3));
                title(['Incremental --> Step: ',num2str(step),'.']);
                hold on;
                if ~isempty(Groundtruth)
%                     Groundtruth = Groundtruth(:,1:step);
                else
                    Groundtruth = [];
                end
                arrowLength = 0.3;
                
                pose_part = Xstate(Xstate(:,2)==1);
                reshape_pose = reshape(pose_part,3,step); % reshape pose_xstate
                fea_part = Xstate(Xstate(:,2)==2);
                num_fea = size(unique(Xstate(Xstate(:,2)==2,3)),1);
                Fea = reshape(fea_part,5,num_fea); % reshape Feature_xstate
                

                %
                if ~isempty(scan)
                    Thisscan = transScan(scan(2:3,:),pose_part(end-2:end));
                    obj(1).figLaser = plot(Thisscan(1,:),Thisscan(2,:),'g.');
                else
                    obj(1).figLaser = [];
                end
                obj(1).figFea = plotFeature(Fea,color,'notext'); % Draw estimated features
                hold on;
%                 for i = 1: length(obj(1).figFea)
%                     fill(obj(1).figFea(i).XData, obj(1).figFea(i).YData,color,'EdgeColor','none','FaceAlpha',0.7);
%                 end
%                 obj(1).figFea = [];
                % Instead of arrow-----draw groundtruth
                if ~isempty(Groundtruth)
                    theta = Groundtruth(3,:);
                    dx = Groundtruth(1,:) + arrowLength * cos(theta);
                    dy = Groundtruth(2,:) + arrowLength * sin(theta);
                    gt1 = plot([Groundtruth(1,:);dx],[Groundtruth(2,:);dy],'k','LineWidth',2);
                end
                % Instead of arrow-----draw estimated path
                theta1 = reshape_pose(3,:);
                dx1 = reshape_pose(1,:) + arrowLength*cos(theta1);
                dy1 = reshape_pose(2,:) + arrowLength*sin(theta1);
                ps1 = plot([reshape_pose(1,:);dx1],[reshape_pose(2,:);dy1],'Color',color,'LineWidth',1.2);
%                 ps1 = [];%for real exp
                if ~isempty(Groundtruth)
                    % connect
%                     obj(1).figLink = plot([reshape_pose(1,:);Groundtruth(1,:)],[reshape_pose(2,:);Groundtruth(2,:)],'b--');
                    % Draw path point for Groundtruth and estimation
                    gt2 = plot(Groundtruth(1,:),Groundtruth(2,:),'ko','MarkerFaceColor','k','MarkerSize',4); % Draw true positions
                    gt3 = plot(Groundtruth(1,:),Groundtruth(2,:),'k','LineWidth',2); % Draw true positions
                    obj(1).figGt=gt1;
                    obj(2).figGt=gt2;
                    obj(3).figGt=gt3;
                end
                ps2 = plot(reshape_pose(1,:),reshape_pose(2,:),'o','Color',color,'Marker','o','MarkerFaceColor',color,'MarkerSize',4);
                ps3 = plot(reshape_pose(1,:),reshape_pose(2,:),'Color',color,'LineStyle','-','LineWidth',linewidth);
                % for real exp
%                 ps2 = [];
%                 ps3 = plot((reshape_pose(1,:)-27.818277505778530),(reshape_pose(2,:)+13.562028973949493),'Color',color,'LineStyle','-','LineWidth',linewidth);%only for real data

                %--------------------------Draw COV----------------------%
                %--For the sake of speedup, we don't draw 3 sigma bound.---------%
                if nargin <7
                    colorCOV=[0.57,0.62,0.66];
                    sdwidth = 2;
                    filterSize = 1;
                end
                if ~isempty(cov)
                    cov_pose = cov(1:length(pose_part),1:length(pose_part));
                    cov_fea = cov(length(pose_part)+1:end,length(pose_part)+1:end);
                    for i = 1:filterSize:step
                        poxyCov = cov_pose(3*i-2:3*i-1,3*i-2:3*i-1);
                        obj(1).figCov{i}=drawCOV(pose_part(3*i-2:3*i-1,1),poxyCov,colorCOV,sdwidth);
                        % real exp only
%                         obj(1).figCov{i}=drawCOV(pose_part(3*i-2:3*i-1,1)-[27.818277505778530;-13.562028973949493],poxyCov,colorCOV,sdwidth);%for real experiment
%                        drawnow
                    end
                    for numfea = 1:num_fea
%                         xyCov = cov_fea(5*numfea-4:5*numfea-3,5*numfea-4:5*numfea-3);
%                         obj(2).figCov{numfea}=drawCOV(Fea(1:2,numfea),xyCov,colorCOV,sdwidth);
                        % real exp only
%                         obj(2).figCov{numfea}=drawCOV(Fea(1:2,numfea)- [27.818277505778530;-13.562028973949493],xyCov,colorCOV,sdwidth);
                    end
                end
                
                % warning on;
                
                obj(1).figPs=ps1;
                obj(2).figPs=ps2;
                obj(3).figPs=ps3;
            end
        end
        
    end
end
