function [indexPairs,Toltalinliers,IfUniform] = ParFeatureMatching_FindGoodEpipolar_checkUniform_alpha1(features1,features2,validBlobs1,validBlobs2,I2,Indices)
%{
2018/04/14
ParFeatureMatching_FindGoodEpipolar_alpha2
1. For uniform distribution
2. Add block

%}
Samplefeatures1 = features1(:,Indices);
SamplevalidBlobs1 = validBlobs1(Indices,:);
%
[indexPairs, Scores] = vl_ubcmatch(Samplefeatures1,features2);
indexPairs = indexPairs';

%%

%%
SampleLocation1 = SamplevalidBlobs1.Location;
Location2 = validBlobs2.Location;
MatchedSampleLocation1 = SampleLocation1(indexPairs(:,1),:);
MatchedLocation2 = Location2(indexPairs(:,2),:);
[indexPairs,~] = SIFT_RemoveSameMatching_alpha1(MatchedSampleLocation1,MatchedLocation2,Scores,indexPairs);


MatchedSampleLocation1 = SampleLocation1(indexPairs(:,1),:);
MatchedLocation2 = Location2(indexPairs(:,2),:);
%%
IfUniform = false;
BlockNum = 5;
One2OneNum = size(indexPairs,1);
Toltalinliers = false(One2OneNum,1);

% [~, inliers] = RANSAC_alpha1(MatchedSampleLocation1, MatchedLocation2, One2OneNum, 0.00001);
% Toltalinliers(inliers) = true;

%
Mdl_MatchedSampleLocation1Row = KDTreeSearcher(MatchedSampleLocation1(:,2));
Min_MatchedSampleLocation1Row = min(Mdl_MatchedSampleLocation1Row.X);
Max_MatchedSampleLocation1Row = max(Mdl_MatchedSampleLocation1Row.X);

Interval = (Max_MatchedSampleLocation1Row - Min_MatchedSampleLocation1Row)/BlockNum;
% NumRef = 1:One2OneNum;
% Curr_NumRefDebug = [];
inliersNum = zeros(1,BlockNum);
for i = 1:BlockNum
    Curr_RowMin = Min_MatchedSampleLocation1Row + Interval*(i-1);
    if i ~= BlockNum
        Curr_RowMax = Min_MatchedSampleLocation1Row + Interval*i;
    else
        Curr_RowMax = Min_MatchedSampleLocation1Row + (Interval*i) + 1;
    end
    Curr_Num = (Mdl_MatchedSampleLocation1Row.X >= Curr_RowMin) & (Mdl_MatchedSampleLocation1Row.X < Curr_RowMax);
%     Curr_NumRef = NumRef(Curr_Num);
%     Curr_NumRefDebug = [Curr_NumRefDebug Curr_NumRef];
    [~, inliers] = RANSAC_alpha1(MatchedSampleLocation1(Curr_Num,:), MatchedLocation2(Curr_Num,:), sum(Curr_Num), 0.00001);
    inliersNum(i) = sum(inliers);
    Toltalinliers(Curr_Num) = inliers;
end

if sum(Toltalinliers) > 8
    IfUniform = CheckUniform_alpha2(double(MatchedSampleLocation1),Toltalinliers,0.75);
end
% 
% if IfUniform
%     RegistedMatchedSampleLocation1 = MatchedSampleLocation1(Toltalinliers,:);
%     RegistedMatchedLocation2 = MatchedLocation2(Toltalinliers,:);
%     Best_F = RANSAC_Norm8Point(RegistedMatchedSampleLocation1', RegistedMatchedLocation2');
%     [t1, t2] = estimateUncalibratedRectification(Best_F, RegistedMatchedSampleLocation1, RegistedMatchedLocation2, size(I2));
%     tform1 = projective2d(t1);
%     tform2 = projective2d(t2);
% else
%     tform1 = [];
%     tform2 = [];
% end

% Length1 = max(SampleLocation1(:,1)) - min(SampleLocation1(:,1));
% Width1 = max(SampleLocation1(:,2)) - min(SampleLocation1(:,2));
% 
% LengthMatched1 = max(MatchedSampleLocation1(Toltalinliers,1)) - min(MatchedSampleLocation1(Toltalinliers,1));
% WidthMatched1 = max(MatchedSampleLocation1(Toltalinliers,2)) - min(MatchedSampleLocation1(Toltalinliers,2));
% 
% LengthRatio1 = LengthMatched1/Length1;
% WidthRatio1 = WidthMatched1/Width1;
% if ~IfGood
%     InliersDifferSTD = 100;
%     tform1 = [];
%     tform2 = [];
%     indexPairs = [];
%     Toltalinliers = [];
%    return 
% end
% matchedPoints1 = SamplevalidBlobs1(indexPairs(Toltalinliers,1),:);
% matchedPoints2 = validBlobs2(indexPairs(Toltalinliers,2),:);
% figure,showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2);
%{
RegistedMatchedSampleLocation1 = MatchedSampleLocation1(Toltalinliers,:);
RegistedMatchedLocation2 = MatchedLocation2(Toltalinliers,:);
Best_F = RANSAC_Norm8Point(RegistedMatchedSampleLocation1', RegistedMatchedLocation2');
[t1, t2] = estimateUncalibratedRectification(Best_F, RegistedMatchedSampleLocation1, RegistedMatchedLocation2, size(I2));
tform1 = projective2d(t1);
tform2 = projective2d(t2);

RegistedMatchedSampleLocation1out = transformPointsForward(tform1, RegistedMatchedSampleLocation1);
RegistedMatchedLocation2out = transformPointsForward(tform2, RegistedMatchedLocation2);

InliersDiffer = RegistedMatchedSampleLocation1out(:,2) - RegistedMatchedLocation2out(:,2);
InliersDifferSTD = std(InliersDiffer);


% figure,
% plot(RegistedMatchedSampleLocation1out(:,1),RegistedMatchedSampleLocation1out(:,2),'yo');hold on
% plot(RegistedMatchedLocation2out(:,1),RegistedMatchedLocation2out(:,2),'g*');
% for i = 1:sum(inliers)
%     line([RegistedMatchedSampleLocation1out(i,1) RegistedMatchedLocation2out(i,1)],[RegistedMatchedSampleLocation1out(i,2) RegistedMatchedLocation2out(i,2)],'Color','k');
% end
%%
%{
[Best_F, inliers] = RANSAC_alpha1(MatchedSampleLocation1, MatchedLocation2, 2000, 0.00001);%0.00001
%[Accuracy,Recall,Precision,Specificity,ResultingInliersNum,ResultingOutliersNum,ResultingOutliersRatio] = Evaluate_alpha1(MatchedLocation1,inliers,TrueInliers,TrueOutiers);

% matchedPoints1 = SamplevalidBlobs1(indexPairs(inliers,1),:);
% matchedPoints2 = validBlobs2(indexPairs(inliers,2),:);
% figure,
% showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2);
%%
RegistedMatchedSampleLocation1 = MatchedSampleLocation1(inliers,:);
RegistedMatchedLocation2 = MatchedLocation2(inliers,:);


[t1, t2] = estimateUncalibratedRectification(Best_F, RegistedMatchedSampleLocation1, RegistedMatchedLocation2, size(I2));

tform1 = projective2d(t1);
tform2 = projective2d(t2);

RegistedMatchedSampleLocation1out = transformPointsForward(tform1, RegistedMatchedSampleLocation1);
RegistedMatchedLocation2out = transformPointsForward(tform2, RegistedMatchedLocation2);

InliersDiffer = RegistedMatchedSampleLocation1out(:,2) - RegistedMatchedLocation2out(:,2);
InliersDifferSTD = std(InliersDiffer);

%%
matchedPoints1 = SamplevalidBlobs1(indexPairs(inliers,1),:);
matchedPoints2 = validBlobs2(indexPairs(inliers,2),:);
figure,showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2);
% plot(RegistedMatchedSampleLocation1out(:,1),RegistedMatchedSampleLocation1out(:,2),'yo');hold on
% plot(RegistedMatchedLocation2out(:,1),RegistedMatchedLocation2out(:,2),'g*');
% for i = 1:sum(inliers)
%     line([RegistedMatchedSampleLocation1out(i,1) RegistedMatchedLocation2out(i,1)],[RegistedMatchedSampleLocation1out(i,2) RegistedMatchedLocation2out(i,2)],'Color','k');
% end
% 
% x =1;
%}

%}