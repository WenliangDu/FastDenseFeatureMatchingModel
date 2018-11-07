function SearchingRange = ParFeatureMatching_CalcuSearchingRange_alpha2(Inliers_ProjectedLocation1,EpiLines1ab)
%{
2018/09/02
ParFeatureMatching_CalcuSearchingRange_alpha2
1. Calculate the searching range depending on the epipolar line which has
the maximum slope and on the projected features in I1.

2018/08/31
ParFeatureMatching_CalcuSearchingRange_alpha1
1. Calculate the searching range depending on the epipolar line which has
the maximum slope and on the projected size of I1.

%}
% numRows = SizeI1(1);
% numCols = SizeI1(2);
% inPts = [1, 1; 1, numRows; numCols, numRows; numCols, 1];
% outPts = transformPointsForward(tform1, inPts);

[MaxSlope,MaxSlopeN] = max(abs(EpiLines1ab(:,1)));
% xL = [Left_Most1 Right_Most1];
EL = polyval(EpiLines1ab(MaxSlopeN,:),[max(Inliers_ProjectedLocation1(:,1)) min(Inliers_ProjectedLocation1(:,1))]);

SearchingRange = abs(EL(1)-EL(2));