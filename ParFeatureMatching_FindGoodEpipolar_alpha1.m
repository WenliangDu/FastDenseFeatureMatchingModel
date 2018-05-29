function [InliersDifferSTD,tform1,tform2,indexPairs,inliers] = ParFeatureMatching_FindGoodEpipolar_alpha1(features1,features2,validBlobs1,validBlobs2,I2,Indices)


Samplefeatures1 = features1(:,Indices);
SamplevalidBlobs1 = validBlobs1(Indices,:);
%
[indexPairs, Scores] = vl_ubcmatch(Samplefeatures1,features2);
indexPairs = indexPairs';


SampleLocation1 = SamplevalidBlobs1.Location;
Location2 = validBlobs2.Location;
MatchedSampleLocation1 = SampleLocation1(indexPairs(:,1),:);
MatchedLocation2 = Location2(indexPairs(:,2),:);
[indexPairs,~] = SIFT_RemoveSameMatching_alpha1(MatchedSampleLocation1,MatchedLocation2,Scores,indexPairs);


MatchedSampleLocation1 = SampleLocation1(indexPairs(:,1),:);
MatchedLocation2 = Location2(indexPairs(:,2),:);
%
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