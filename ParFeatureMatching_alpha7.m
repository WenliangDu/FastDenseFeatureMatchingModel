function [indexPairsPar,RecordIndiesNew,tform1,tform2,IndicesInitialRANSAC,indexPairsInitialRANSAC,inliersInitialRANSAC,GoodClassidx,SearchingRange,RealSearchingRange] =ParFeatureMatching_alpha7(features1,features2,validBlobs1,validBlobs2,I1,I2)
%{
2018/09/11
ParFeatureMatching_alpha7
1. For getting numerous matching areas

2018/09/10
ParFeatureMatching_alpha7_checkUniformAndSR
1. Check if the searching range program could be work in the real images.
2. Find good sampling with accordant transmation angle by k-means

2018/09/02
ParFeatureMatching_alpha5_checkUniformAndSR
1. Check if the searching range program could be work in the real images.


2018/04/17
ParFeatureMatching_alpha5
1. ParFeatureMatching_LimitInliersDifferSTD_alpha1 =>
ParFeatureMatching_LimitInliersDifferSTD_alpha2
Dynamic P

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


Location1 = double(validBlobs1.Location);
Location2 = double(validBlobs2.Location);
RandomNum = 100;
IdealMeanRatios = zeros(1,RandomNum);
RealMeanRatios = zeros(1,RandomNum);
IfGoods = zeros(1,RandomNum);
% LengthRatios = zeros(1,RandomNum);
% WidthRatios = zeros(1,RandomNum);
IfGoodRatios = zeros(1,RandomNum);

% WrongSampledLocation1 = cell(1,1);
% WrongMatchedLocation1 = cell(1,1);
% WrongNum = 1;


IterationNum = 300;
Curr_IterationNum = 1;
QualifiedSampling = false;
while ~QualifiedSampling && (Curr_IterationNum <= IterationNum)
    IndicesInitialRANSAC = randsample(features1Num, SamplingNum);
%     SampledLocation1 = Location1(IndicesInitialRANSAC,:);
    [indexPairsInitialRANSAC,inliersInitialRANSAC,IfUniform] = ParFeatureMatching_FindGoodEpipolar_checkUniform_alpha1(features1,features2,validBlobs1,validBlobs2,I2,IndicesInitialRANSAC);
    if IfUniform
        %% Here!!!
        Matched_Location1 = Location1(IndicesInitialRANSAC(indexPairsInitialRANSAC(inliersInitialRANSAC,1)),:);
        Matched_Location2 = Location2(indexPairsInitialRANSAC(inliersInitialRANSAC,2),:);
        x = 1;
        [GoodClassidx,tform1,tform2,QualifiedSampling,Qualified_Matched_Location1,Qualified_Matched_Location2,SearchingRange,RealSearchingRange,k,EpiLines1ab,EpiLines2ab] = ParFeatureMatching_FindGoodEmission_alpha2(Matched_Location1,Matched_Location2,size(I2),Location1);
        x = 1;
    end
end

%%
Location1Rect = transformPointsForward(tform1, Location1);
Location2Rect = transformPointsForward(tform2, Location2);

Qualified_Matched_Location1Rect = transformPointsForward(tform1, Qualified_Matched_Location1);
Qualified_Matched_Location2Rect = transformPointsForward(tform2, Qualified_Matched_Location2);

% figure,hold on
% plot(Qualified_Matched_Location1Rect(:,1),Qualified_Matched_Location1Rect(:,2),'r*');hold on
% plot(Qualified_Matched_Location2Rect(:,1),Qualified_Matched_Location2Rect(:,2),'b*');
% line([Qualified_Matched_Location1Rect(:,1)';Qualified_Matched_Location2Rect(:,1)'],[Qualified_Matched_Location1Rect(:,2)';Qualified_Matched_Location2Rect(:,2)'],'Color','g');

% Diff = abs(Qualified_Matched_Location1Rect(:,2) - Qualified_Matched_Location2Rect(:,2));
% WrongIdx = Diff>=SearchingRange;
% Matched_Location1Rect(WrongIdx,:) = [];
% Matched_Location2Rect(WrongIdx,:) = [];
% 
% figure,hold on
% plot(Matched_Location1Rect(:,1),Matched_Location1Rect(:,2),'r*');hold on
% plot(Matched_Location2Rect(:,1),Matched_Location2Rect(:,2),'b*');
% line([Matched_Location1Rect(:,1)';Matched_Location2Rect(:,1)'],[Matched_Location1Rect(:,2)';Matched_Location2Rect(:,2)'],'Color','g');

[indexPairsPar,~,RecordIndies] = ParFeatureMatching_ParMatching_alpha3(features1,features2,Location1Rect,Location2Rect,RealSearchingRange,I2,Qualified_Matched_Location1Rect,Qualified_Matched_Location2Rect);


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
%%

