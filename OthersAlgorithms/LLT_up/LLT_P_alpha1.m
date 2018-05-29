function [inliersLLTA,inliersLLTR,inliersLLTV,LLTAMatchingTime,LLTRMatchingTime,LLTVMatchingTime] = LLT_P_alpha1(Location1,Location2)
%{
2016/03/20
LLT_alpha1
1. For OtherAlgorightms_alpha5
2. For nonrigid

2016/11/26
LLT_P_alpha1
1. For processing in batch.

%}
inliersLLTA = [];
inliersLLTR = [];
inliersLLTV = [];
%% conf
conf.lambda = 9 * (10^9);
conf.Kn = 15;
if ~exist('conf'), conf = []; end
conf = LLT_init(conf);

%% Norm
% X = [Location1(:,1)';Location1(:,2)']';
% Y = [Location2(:,1)';Location2(:,2)']';
[nX, nY, ~]=norm2(Location1,Location2);

%% LLTA
% sprintf('LLTA:')
t1=clock;
VecFldA=LLTA(nX, nY, conf);
inliersLLTA = VecFldA.VFCIndex;
t2=clock;
LLTAMatchingTime = etime(t2,t1);
%% LLTR
% sprintf('LLTR:')
t1=clock;
VecFldR=LLTR(nX, nY, conf);
inliersLLTR = VecFldR.VFCIndex;
t2=clock;
LLTRMatchingTime = etime(t2,t1);
%% LLTV
% sprintf('LLTV:')
t1=clock;
VecFldV=LLTV(nX, nY, conf);
inliersLLTV = VecFldV.VFCIndex;
t2=clock;
LLTVMatchingTime = etime(t2,t1);