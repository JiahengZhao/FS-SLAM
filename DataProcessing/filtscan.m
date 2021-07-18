function scan = filtscan(scan,numFilter)
% FILTSCAN filter some scan and odometry.

numLength = length(scan.scan);
if numFilter == 0 || numFilter == 1
    return
end
filterID = 1:numFilter:numLength;
if filterID(end) ~= numLength
    check = true;
    filterID = [filterID numLength];
else
    check = false;
end
scan.scan = scan.scan(1,filterID);

for i = filterID
    if i==1
        P_last = scan.odom(1,:);
        P_first = P_last;
        continue;
    end
    P_this = scan.odom(i,:);
    
    if i==filterID(end-1) && check
        tmp=zeros(1,3);
        for j = 1:(filterID(end)-filterID(end-1)+1)
            P_curr =  scan.odom(i+j-numFilter,:);
            R_last = theta2R(P_last(3));
            tmp(1,1:2) = tmp(1,1:2) + P_curr(1:2)*R_last';
            tmp(1,3) = tmp(1,3) + P_curr(3);
            P_last = P_curr;
        end
        R_first = theta2R(P_first(3));
        P_this(1:2) = tmp(1:2) * R_first';
        P_this(3) = tmp(3);
    else
        tmp=zeros(1,3);
        for j = 1:numFilter
            P_curr =  scan.odom(i+j-numFilter,:);
            R_last = theta2R(P_last(3));
            tmp(1,1:2) = tmp(1,1:2) + P_curr(1:2)*R_last';
            tmp(1,3) = tmp(1,3) + P_curr(3);
            P_last = P_curr;
        end
        R_first = theta2R(P_first(3));
        P_this(1:2) = tmp(1:2) * R_first';
        P_this(3) = tmp(3);
    end
    P_first = scan.odom(i,:);
    scan.odom(i,:) = P_this;
    
end
scan.odom = scan.odom(filterID,:);
scan.groundtruth = scan.groundtruth(:,filterID);