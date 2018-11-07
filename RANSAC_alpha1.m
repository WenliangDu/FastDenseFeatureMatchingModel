function [Best_F, inliers] = RANSAC_alpha1(Location1, Location2, Trials_N, Threshold)

if size(Location1,1) > size(Location1,2)
    Location1 = Location1';
    Location2 = Location2';
end

if size(Location1,1) < 3
    Location1(3,:) = 1;
    Location2(3,:) = 1;
end

Location_N = size(Location1,2);
inliers = false(1, Location_N);
Sample_N = 8;
DistType = 'Algebraic';%'sampson''Algebraic'

Best_F = [];
%integerClass = 'int32';
% Confidence % For update MaxTrials_N
if Location_N > Sample_N
    MaxTrials_N = Trials_N;
    CurTrials_N = int32(0);
    BestInliers_N = int32(0);
%     logOneMinusConf = log(1 - Confidence); % For update MaxTrials_N
%     oneOverNPts = double(1 / Location_N); % For update MaxTrials_N
    
    while CurTrials_N < MaxTrials_N
        Indices = randsample(Location_N, Sample_N);
        Current_F = RANSAC_Norm8Point(Location1(:,Indices), Location2(:,Indices));
        Current_Distance = Calculate_Dist_alpha1(DistType, Location1, Location2, Current_F);
        
        %[curInliers, curNInliers] = findInliers(Current_Distance, Sample_N, Threshold);
        %% Find inliers
        CurIinliers = Current_Distance <= Threshold;
        CurInliers_N = int32(sum(CurIinliers));
        
        if BestInliers_N < CurInliers_N
            BestInliers_N = CurInliers_N;
            inliers = CurIinliers;
            % Update the number of trials
%             maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, ...
%                 outputClass, integerClass, curNInliers, maxNTrials);
        end
        CurTrials_N = CurTrials_N + 1;
    end
    
    if BestInliers_N >= 8
        Best_F = RANSAC_Norm8Point(Location1(:, inliers), Location2(:, inliers));
    else
        Best_F = [];
    end

else
    % error
end