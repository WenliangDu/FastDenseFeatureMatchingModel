function [inliersLLTA,inliersLLTR,inliersLLTV] = LLT_alpha1(Location1,Location2)
%{
2016/03/20
LLT_alpha1
1. For OtherAlgorightms_alpha5
2. For nonrigid

%}
inliersLLTA = [];
inliersLLTR = [];
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
sprintf('LLTA:')
tic,
VecFldA=LLTA(nX, nY, conf);
inliersLLTA = VecFldA.VFCIndex;
toc,
%% LLTR
sprintf('LLTR:')
tic,
VecFldR=LLTR(nX, nY, conf);
inliersLLTR = VecFldR.VFCIndex;
toc,
%% LLTV
sprintf('LLTV:')
tic,
VecFldV=LLTV(nX, nY, conf);
inliersLLTV = VecFldV.VFCIndex;
toc,