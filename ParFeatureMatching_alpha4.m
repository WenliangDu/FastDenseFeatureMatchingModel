function [indexPairsPar,RecordIndiesNew,BestInliersDifferSTD,tform1,tform2,IndicesInitialRANSAC,indexPairsInitialRANSAC,inliersInitialRANSAC] = ParFeatureMatching_alpha4(features1,features2,validBlobs1,validBlobs2,I2,ThresholdD)
%{
2018/04/14
ParFeatureMatching_alpha4
1. ParFeatureMatching_FindGoodEpipolar_alpha1 = > ParFeatureMatching_FindGoodEpipolar_alpha2
2. For uniform distribution Add block


2018/04/14
ParFeatureMatching_alpha3_debug
1. For debuging thresholdD = 10


2018/04/08
ParFeatureMatching_alpha2
1. ParMatching for real stereo images
2. Add limit of iteration like ParFeatureMatching_FindPar_alpha1

2018/04/12
ParFeatureMatching_alpha3
1. Add ParFeatureMatching_LimitInliersDifferSTD_alpha1
2. ParFeatureMatching_ParMatching_alpha1 => ParFeatureMatching_ParMatching_alpha2
Add i=1 and i = Match_Num
3. Update RecordIndies

%}

SamplingNum = 1000;
features1Num = size(features1,2);
%%
IndicesInitialRANSAC = randsample(features1Num, SamplingNum);
[BestInliersDifferSTD,Besttform1,Besttform2,indexPairsInitialRANSAC,inliersInitialRANSAC] = ParFeatureMatching_FindGoodEpipolar_alpha2(features1,features2,validBlobs1,validBlobs2,I2,IndicesInitialRANSAC);
LimitInliersDifferSTD = ParFeatureMatching_LimitInliersDifferSTD_alpha1(length(indexPairsInitialRANSAC));

IterationNum = SamplingNum;
Curr_IterationNum = 1;
while (BestInliersDifferSTD >= LimitInliersDifferSTD) && (Curr_IterationNum <= IterationNum)
    
    IndicesInitialRANSAC = randsample(features1Num, SamplingNum);
    [InliersDifferSTD,tform1,tform2,indexPairsInitialRANSAC,inliersInitialRANSAC] = ParFeatureMatching_FindGoodEpipolar_alpha2(features1,features2,validBlobs1,validBlobs2,I2,IndicesInitialRANSAC);
    LimitInliersDifferSTD = ParFeatureMatching_LimitInliersDifferSTD_alpha1(length(indexPairsInitialRANSAC));
    if InliersDifferSTD < BestInliersDifferSTD
        BestInliersDifferSTD = InliersDifferSTD;
        Besttform1 = tform1;
        Besttform2 = tform2;
    end
    Curr_IterationNum = Curr_IterationNum + 1;
end
if isempty(Besttform1)
    indexPairsPar = [];
    RecordIndiesNew = [];
    return
else
    tform1 = Besttform1;
    tform2 = Besttform2;
end


Location1 = validBlobs1.Location;
Location2 = validBlobs2.Location;
% 
Location1Rect = transformPointsForward(tform1, Location1);
Location2Rect = transformPointsForward(tform2, Location2);
% 
% figure,
% plot(Location1Rect(:,1),Location1Rect(:,2),'ro');hold on
% plot(Location2Rect(:,1),Location2Rect(:,2),'b*');
% 
% Location1RectM = transformPointsForward(tform1, Location1(indexPairs(inliers,1),:));
% Location2RectM = transformPointsForward(tform2, Location2(indexPairs(inliers,2),:));
% 
% InliersDifferAll = Location1RectM(:,2) - Location2RectM(:,2);
% InliersDifferSTDAll = std(InliersDifferAll);
% InliersDifferAllabsSort = sort(abs(InliersDifferAll));
% ReNum = round(0.90*length(inliers));
% GoodThresholdD = InliersDifferAllabsSort(ReNum);
% % matchedPoints1 = validBlobs1(indexPairs(:,1),:);
% % matchedPoints2 = validBlobs2(indexPairs(:,2),:);
% 
% plot(Location1RectM(:,1),Location1RectM(:,2),'yo');hold on
% plot(Location2RectM(:,1),Location2RectM(:,2),'g*');
% for i = 1:length(inliers)
%     line([Location1RectM(i,1) Location2RectM(i,1)],[Location1RectM(i,2) Location2RectM(i,2)],'Color','k');
% end
% 
% Location1RectR = transformPointsForward(tform1, Location1(IndicesInitialRANSAC(indexPairsInitialRANSAC(inliersInitialRANSAC,1)),:));
% Location2RectR = transformPointsForward(tform2, Location2(indexPairsInitialRANSAC(inliersInitialRANSAC,2),:));
% figure,
% plot(Location1RectR(:,1),Location1RectR(:,2),'yo');hold on
% plot(Location2RectR(:,1),Location2RectR(:,2),'g*');
% for i = 1:sum(inliersInitialRANSAC)
%     line([Location1RectR(i,1) Location2RectR(i,1)],[Location1RectR(i,2) Location2RectR(i,2)],'Color','k');
% end
% figure,showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2);
%
[indexPairsPar,~,RecordIndies] = ParFeatureMatching_ParMatching_alpha2(features1,features2,Location1Rect,Location2Rect,ThresholdD,I2);
MatchedLocation1Temp = Location1(indexPairsPar(:,1),:);
MatchedLocation2Temp = Location2(indexPairsPar(:,2),:);
[indexPairsPar,RecordIndies] = SIFT_RemoveSameMatching_alpha1(MatchedLocation1Temp,MatchedLocation2Temp,RecordIndies,indexPairsPar);

%% RecordIndiesNew
% RecordIndiesNum = length(RecordIndies);
ValidRecordRef = unique(RecordIndies);
ValidRecordRefL = length(ValidRecordRef);
RecordIndiesNew = zeros(ValidRecordRefL,2);

%% j == 1
RecordIndiesNew(1,1) = 1;
RecordIndiesNew(1,2) = RecordIndiesNew(1,1) + find(RecordIndies(RecordIndiesNew(1,1):end) ~= ValidRecordRef(1),1) - 2;
for j = 2:ValidRecordRefL-1
    RecordIndiesNew(j,1) = RecordIndiesNew(j-1,2) + 1;
    RecordIndiesNew(j,2) = RecordIndiesNew(j,1) + find(RecordIndies(RecordIndiesNew(j,1):end) ~= ValidRecordRef(j),1) - 2;
end
RecordIndiesNew(ValidRecordRefL,1) = RecordIndiesNew(ValidRecordRefL-1,2) + 1;
RecordIndiesNew(ValidRecordRefL,2) = length(RecordIndies);
%}
