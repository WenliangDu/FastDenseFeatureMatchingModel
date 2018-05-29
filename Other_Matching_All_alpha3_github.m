function inliers = Other_Matching_All_alpha3_github(Locations,method)
%{
2018/04/13
Other_Matching_All_alpha3


2018/04/04
Other_Matching_All_alpha1
1. Based on Other_Modeling_All_alpha1
2. Only retain matching


2016/12/05
Other_Modeling_All_alpha1
1. Adapt DMA dynamic matching area to advanced algorithms.

method 1: LLTA
method 2: LLTR
method 3: LLTV
method 4: SOCBV
method 5: WGTM
%}
MatchedLocation1 = Locations(:,1:2);
MatchedLocation2 = Locations(:,3:4);
if method == 1 || method == 2 || method == 3,
    %% mode 1: LLTA; 2: LLTR; 3: LLTV.
    mode = method;
    inliers = LLT_DMA_sub_alpha1(MatchedLocation1,MatchedLocation2,mode);
elseif method == 4,
    %% SOCBV
    K = 20;
    [inliers,~] = SOCBV(MatchedLocation1,MatchedLocation2,K,0.2);
elseif method == 5,
    %% WGTM
    K = 20;
    [inliers,~] = WGTM_P_function_alpha1(MatchedLocation1,MatchedLocation2,K,0.2);
end


