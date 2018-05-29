function inliers = LLT_DMA_sub_alpha1(Location1,Location2,mode)
%{
2016/12/05
LLT_DMA_alpha1
1. LLT for DMA

mode 1: LLTA
mode 2: LLTR
mode 3: LLTV
%}

%% conf
conf.lambda = 9 * (10^9);
%% Dynamic conf.Kn and M 2016/12/6
if size(Location1,1) <= 15,
    conf.Kn = 5;
    conf.M = size(Location1,1); %For LLTV
else
    conf.Kn = 15;
end
%conf.Kn = 5;
%%
if ~exist('conf'), conf = []; end
conf = LLT_init(conf);

%% Norm
% X = [Location1(:,1)';Location1(:,2)']';
% Y = [Location2(:,1)';Location2(:,2)']';
[nX, nY, ~]=norm2(Location1,Location2);
if mode == 1,
    %% LLTA
    % sprintf('LLTA:')
    VecFldA=LLTA(nX, nY, conf);
    inliers = VecFldA.VFCIndex;
elseif mode == 2,
    %% LLTR
    % sprintf('LLTR:')
    VecFldR=LLTR(nX, nY, conf);
    inliers = VecFldR.VFCIndex;
else
    %% LLTV
    % sprintf('LLTV:')
    VecFldV=LLTV(nX, nY, conf);
    inliers = VecFldV.VFCIndex;
end