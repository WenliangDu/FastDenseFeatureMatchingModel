function IfUniform = CheckUniform_alpha2(Location1,Indices,Ratio)
%{
2018/09/06
CheckUniform_alpha2
1. add ratio to IdealMeanRatio


2018/08/29
CheckUniform_alpha1
1. sub-program for checking if uniform distribution
%}

Location1DT = delaunay(Location1);
TriDistance = [sqrt(sum((Location1(Location1DT(:,1),:)-Location1(Location1DT(:,2),:)).^2,2)) sqrt(sum((Location1(Location1DT(:,1),:)-Location1(Location1DT(:,3),:)).^2,2)) sqrt(sum((Location1(Location1DT(:,2),:)-Location1(Location1DT(:,3),:)).^2,2))];
MeanMatrixOriginal = mean(TriDistance(:));

Location1Sampled = Location1(Indices,:);
Location1SampledDT = delaunay(Location1Sampled);
TriDistanceSampled = [sqrt(sum((Location1Sampled(Location1SampledDT(:,1),:)-Location1Sampled(Location1SampledDT(:,2),:)).^2,2)) sqrt(sum((Location1Sampled(Location1SampledDT(:,1),:)-Location1Sampled(Location1SampledDT(:,3),:)).^2,2)) sqrt(sum((Location1Sampled(Location1SampledDT(:,2),:)-Location1Sampled(Location1SampledDT(:,3),:)).^2,2))];
MeanMatrixSampled = mean(TriDistanceSampled(:));

IdealMeanRatio = sqrt(size(Location1,1)./size(Location1Sampled,1));
RealMeanRatio = MeanMatrixSampled/MeanMatrixOriginal;

IfUniform = RealMeanRatio > (IdealMeanRatio*Ratio);