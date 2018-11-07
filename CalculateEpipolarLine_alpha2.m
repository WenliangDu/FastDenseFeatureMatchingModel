function [EpiLines1ab,EpiLines2ab] = CalculateEpipolarLine_alpha2(Location1,Location2,Location1All,Location2All)
%{
2018/08/29
CalculateEpipolarLine_alpha1
1. sub-program for calculating epipolar lines.

2018/08/31
CalculateEpipolarLine_alpha2
1. sub-program for calculating epipolar lines for all locations.

%}
if size(Location1,1) > size(Location1,2)
    Location1T = Location1';
    Location2T = Location2';
end
FundemantalMatrix = RANSAC_Norm8Point(Location1T,Location2T);

EpiLines1 = (FundemantalMatrix'*([Location2All ones(size(Location1All,1),1)])')';
EpiLines1ab = EpiLines1(:,[1 3]);
EpiLines1ab = (EpiLines1ab./repmat(EpiLines1(:,2),[1 2])).*(-1);

EpiLines2 = (FundemantalMatrix*([Location1All ones(size(Location2All,1),1)])')';
EpiLines2ab = EpiLines2(:,[1 3]);
EpiLines2ab = (EpiLines2ab./repmat(EpiLines2(:,2),[1 2])).*(-1);