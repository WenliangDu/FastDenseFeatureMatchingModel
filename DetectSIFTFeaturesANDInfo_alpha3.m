function [validBlobs1,features1] = DetectSIFTFeaturesANDInfo_alpha3(I1,PeakThresh,edge_thresh)
%{
2016/11/20
DetectSIFTFeaturesANDInfo_alpha2
1. Only detect one image.

2018/04/02
DetectSIFTFeaturesANDInfo_alpha3
1. Add PeakThresh defalt=0; 0 10 20 30
%}
%% Read Image Pair
%peak_thresh = 0;%0 10 20 30


Color1 = size(I1,3);
if Color1 > 1,
	[validBlobs1V,features1] = vl_sift(im2single(rgb2gray(I1)),'edgethresh', edge_thresh);%,'peak_thresh',PeakThresh
else
    [validBlobs1V,features1] = vl_sift(im2single(I1),'edgethresh', edge_thresh) ;
end

%%
Location1(:,1) = validBlobs1V(1,:);
Location1(:,2) = validBlobs1V(2,:);

Scale1 = validBlobs1V(3,:);
Orientation1 = validBlobs1V(4,:);

validBlobs1 = SURFPoints;
Num1 = length(validBlobs1V);
for i = 1:Num1,
    SIFTTemp = SURFPoints(Location1(i,:));
    SIFTTemp.Orientation = Orientation1(i);
	SIFTTemp.Scale = Scale1(i);
    SIFTTemp.Metric = 1;

    validBlobs1 = [validBlobs1;SIFTTemp];
end