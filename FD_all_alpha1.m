
run('vlfeat-0.9.20/toolbox/vl_setup')
addpath(genpath('OthersAlgorithms/LLT_up/'));
addpath(genpath('OthersAlgorithms/WGTM/'));

ThresholdS = 10;
PeakThresh = 0;
EdgeThresh = 40;

method = 3;
% method 1: LLTA
% method 2: LLTR
% method 3: LLTV
% method 4: SOCBV
% method 5: WGTM

%% Feature extraction

t0 = clock;
% I1 = imread( 'Slope20180301_L.jpg' );
% I2 = imread( 'Slope20180301_R.jpg' );
% [validBlobs1,features1] = DetectSIFTFeaturesANDInfo_alpha3(I1,0,40);
% [validBlobs2,features2] = DetectSIFTFeaturesANDInfo_alpha3(I2,0,40);

load('Slope20180301_SIFTandOne2One.mat');
Location1 = validBlobs1.Location;
Location2 = validBlobs2.Location;

%% FD feature matching
t1 = clock;
[indexPairsPar,RecordIndies,BestInliersDifferSTD,tform1,tform2,IndicesInitialRANSAC,indexPairsInitialRANSAC,inliersInitialRANSAC] = ParFeatureMatching_alpha4(features1,features2,validBlobs1,validBlobs2,I2,ThresholdS);
inliersPar = ParFeatureMatching_RemoveOutliers_alpha2_github(Location1,Location2,indexPairsPar,RecordIndies,method);
t2 = clock;
FDMatchingTime = etime(t2,t1);
FeatureExtractionTime = etime(t1,t0);

matchedPoints1 = validBlobs1(indexPairsPar(inliersPar,1),:);
matchedPoints2 = validBlobs2(indexPairsPar(inliersPar,2),:);
figure,showMatchedFeatures(I1, I2, matchedPoints1, matchedPoints2);

%% DTM
ShowDisparity_alpha3(I1, I2, tform1, tform2,inliersPar,indexPairsPar,Location1,Location2);



