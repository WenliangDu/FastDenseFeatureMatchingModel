function inliers = Other_Matching_All_alpha3(Locations,method)
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
method 6: RANSAC
method 7: GTM
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
elseif method == 6,
    %% RANSAC
    [~, inliersSub] = RANSAC_alpha1(MatchedLocation1, MatchedLocation2, 5, 0.00001);
    Num = 1:size(MatchedLocation1,1);
    inliers = Num(inliersSub);
elseif method == 7,
    %% GTM
    [inliers_Location1,~] = GTM_function_ForDMA_alpha1(MatchedLocation1,MatchedLocation2);
    [~,~,inliers] = GTM_LocationtoIndex_alpha1(inliers_Location1,MatchedLocation1,indexPairs);
end

% [Accuracy,Precision,Recall,Specificity,ResultingInliersNum,MissingInliersNum,ResultingOutliersNum,MissingOutliersNum,TrueInliers,TrueOutiers

%% ResultingNum and MissingNum
% ToltalNum = size(MatchedLocation1,1);
% RecordInliers = false(1,ToltalNum);
% RecordInliers(inliers) = true;
% 
% ResultingInliersNum = length(find(RecordInliers(TrueInliers)==1));
% MissingInliersNum = length(find(RecordInliers(TrueInliers)==0));
% 
% ResultingOutliersNum = length(find(RecordInliers(TrueOutiers)==1));
% MissingOutliersNum = length(find(RecordInliers(TrueOutiers)==0));
% 
% ResultingNumAll = ResultingInliersNum + ResultingOutliersNum;
% MissingNumAll = MissingInliersNum + MissingOutliersNum;
% %% Accuracy, Recall, Precision and Specificity
% Accuracy = (ResultingInliersNum + MissingOutliersNum) / (ResultingNumAll + MissingNumAll);
% Recall = ResultingInliersNum / (ResultingInliersNum + MissingInliersNum);
% if ResultingNumAll == 0,
%     Precision = 0;
% else
%     Precision = ResultingInliersNum / ResultingNumAll;
% end
% if MissingOutliersNum + ResultingOutliersNum == 0,
%     Specificity = 0;
% else
%     Specificity = MissingOutliersNum / (MissingOutliersNum + ResultingOutliersNum);
% end

%{
Center1 = [mean(MatchedLocation1(:,1)) mean(MatchedLocation1(:,2))];
Center2 = [mean(MatchedLocation2(:,1)) mean(MatchedLocation2(:,2))];
TransMatrix_ALL = [1 0 (Center2(:,1) - Center1(:,1));0 1 (Center2(:,2) - Center1(:,2));0 0 1];
%TransMatrix_ALL = Transformation_alpha1(MatchedLocation1(TrueInliers,:),MatchedLocation2(TrueInliers,:));
ProjectionSet_MatchedLocation1 = Tranversal_KNN_CalculateProjection_alpha1(MatchedLocation1,TransMatrix_ALL);

figure,
plot(ProjectionSet_MatchedLocation1(TrueInliers,1),ProjectionSet_MatchedLocation1(TrueInliers,2),'r*');hold on
plot(MatchedLocation2(TrueInliers,1),MatchedLocation2(TrueInliers,2),'b*');
plot(ProjectionSet_MatchedLocation1(TrueOutiers,1),ProjectionSet_MatchedLocation1(TrueOutiers,2),'mo');
plot(MatchedLocation2(TrueOutiers,1),MatchedLocation2(TrueOutiers,2),'co');
for i = 1:length(TrueInliers),
    line([ProjectionSet_MatchedLocation1(TrueInliers(i),1) MatchedLocation2(TrueInliers(i),1)],[ProjectionSet_MatchedLocation1(TrueInliers(i),2) MatchedLocation2(TrueInliers(i),2)],'color','g');
end
for i = 1:length(TrueOutiers),
    line([ProjectionSet_MatchedLocation1(TrueOutiers(i),1) MatchedLocation2(TrueOutiers(i),1)],[ProjectionSet_MatchedLocation1(TrueOutiers(i),2) MatchedLocation2(TrueOutiers(i),2)],'color','r');
end
title('Initial');

figure,
plot(ProjectionSet_MatchedLocation1(TrueInliers,1),ProjectionSet_MatchedLocation1(TrueInliers,2),'r*');hold on
plot(MatchedLocation2(TrueInliers,1),MatchedLocation2(TrueInliers,2),'b*');

for i = 1:length(ResultingInliers),
    line([ProjectionSet_MatchedLocation1(ResultingInliers(i),1) MatchedLocation2(ResultingInliers(i),1)],[ProjectionSet_MatchedLocation1(ResultingInliers(i),2) MatchedLocation2(ResultingInliers(i),2)],'color','g');
end
if ~isempty(MissingInliers),
    for i = 1:length(MissingInliers),
        line([ProjectionSet_MatchedLocation1(MissingInliers(i),1) MatchedLocation2(MissingInliers(i),1)],[ProjectionSet_MatchedLocation1(MissingInliers(i),2) MatchedLocation2(MissingInliers(i),2)],'color','r');
    end
end
title('Inliers');

figure,
plot(ProjectionSet_MatchedLocation1(TrueOutiers,1),ProjectionSet_MatchedLocation1(TrueOutiers,2),'mo');hold on
plot(MatchedLocation2(TrueOutiers,1),MatchedLocation2(TrueOutiers,2),'co');
for i = 1:length(MissingOutliers),
    line([ProjectionSet_MatchedLocation1(MissingOutliers(i),1) MatchedLocation2(MissingOutliers(i),1)],[ProjectionSet_MatchedLocation1(MissingOutliers(i),2) MatchedLocation2(MissingOutliers(i),2)],'color','g');
end
if ~isempty(ResultingOutliers),
    for i = 1:length(ResultingOutliers),
        line([ProjectionSet_MatchedLocation1(ResultingOutliers(i),1) MatchedLocation2(ResultingOutliers(i),1)],[ProjectionSet_MatchedLocation1(ResultingOutliers(i),2) MatchedLocation2(ResultingOutliers(i),2)],'color','r');
    end
end
title('Outliers');

x = 1;
%}
%CalculatingOutliers = 

% ResultingOutliers = find(RecordTrueOutliers(RecordProcessedOutliers));
% MissingOutliers = find(RecordProcessedOutliers(RecordTrueOutliers));


