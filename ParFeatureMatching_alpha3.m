function [indexPairsPar,RecordIndiesNew,BestInliersDifferSTD,tform1,tform2,IndicesInitialRANSAC,indexPairsInitialRANSAC,inliersInitialRANSAC] = ParFeatureMatching_alpha3(features1,features2,validBlobs1,validBlobs2,I2,ThresholdD)
%{
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
IndicesInitialRANSAC = randsample(features1Num, SamplingNum);
[BestInliersDifferSTD,Besttform1,Besttform2,indexPairsInitialRANSAC,inliersInitialRANSAC] = ParFeatureMatching_FindGoodEpipolar_alpha1(features1,features2,validBlobs1,validBlobs2,I2,IndicesInitialRANSAC);
% xx = 1;
LimitInliersDifferSTD = ParFeatureMatching_LimitInliersDifferSTD_alpha1(length(indexPairsInitialRANSAC));
% ThresholdD = 3;

IterationNum = SamplingNum;
Curr_IterationNum = 1;
while (BestInliersDifferSTD >= LimitInliersDifferSTD) && (Curr_IterationNum <= IterationNum)
    
    IndicesInitialRANSAC = randsample(features1Num, SamplingNum);
    [InliersDifferSTD,tform1,tform2,indexPairsInitialRANSAC,inliersInitialRANSAC] = ParFeatureMatching_FindGoodEpipolar_alpha1(features1,features2,validBlobs1,validBlobs2,I2,IndicesInitialRANSAC);
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

Location1Rect = transformPointsForward(tform1, Location1);
Location2Rect = transformPointsForward(tform2, Location2);

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
