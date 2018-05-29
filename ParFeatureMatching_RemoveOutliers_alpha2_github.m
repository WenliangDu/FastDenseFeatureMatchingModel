function inliersPar = ParFeatureMatching_RemoveOutliers_alpha2_github(Location1,Location2,indexPairsPar,RecordIndies,method)

%{
2018/04/07
ParFeatureMatching_Model_RemoveOutliers_alpha1
1. remove outliers by 4 methods

2018/04/12
ParFeatureMatching_RemoveOutliers_alpha1
1. Add parfor

method 1: LLTA
method 2: LLTR
method 3: LLTV
method 4: SOCBV
method 5: WGTM
method 6: RANSAC
method 7: GTM
%}
ParInNumber = (RecordIndies(:,2) - RecordIndies(:,1)) + 1;
RecordIndiesNum = size(RecordIndies,1);
ParInWrong = (ParInNumber) < 6;
IfallWrong = false;
ParInRight = ~ParInWrong;
RecordIndiesNew = RecordIndies;
if any(ParInWrong)
    ParInWrongNumber = find(ParInWrong);
    for j = 1:length(ParInWrongNumber)
        Curr_Num = ParInWrongNumber(j);
        NextRight = find(ParInRight(Curr_Num:RecordIndiesNum),1);
        if ~isempty(NextRight)
            NextRight = NextRight + Curr_Num - 1;
            if RecordIndies(Curr_Num,1) < RecordIndiesNew(NextRight,1)
                RecordIndiesNew(NextRight,1) = RecordIndies(Curr_Num,1);
            end
        else
            NextRight = find(ParInRight(1:Curr_Num),1,'last');
            if ~isempty(NextRight)
%                 NextRight = Curr_Num - NextRight + 1;
                RecordIndiesNew(NextRight,2) = RecordIndies(Curr_Num,2);
            else
                %% all wrong
                RecordIndiesNew = [RecordIndies(1,1) RecordIndies(RecordIndiesNum,2)];
                IfallWrong = true;
                break
            end
        end
        
    end
    if ~IfallWrong
        RecordIndiesNew(ParInWrong,:) = [];
    end
end


%%
TraverNum = size(RecordIndiesNew,1);
inliersPar = false(TraverNum,1);
ParNum = zeros(TraverNum,1);
inliersParCell = cell(TraverNum,1);
inliersParIndices = zeros(TraverNum,1);


MatchedLocation = cell(TraverNum,1);

for i = 1:TraverNum
    Curr_indexPairs = indexPairsPar(RecordIndiesNew(i,1):RecordIndiesNew(i,2),:);
    MatchedLocation{i,1} = [Location1(Curr_indexPairs(:,1),:) Location2(Curr_indexPairs(:,2),:)];
    ParNum(i) = RecordIndiesNew(i,2)-RecordIndiesNew(i,1)+1;
end
% for i = 1:TraverNum
parfor (i = 1:TraverNum,4) %% parallel

    
    inliers = Other_Matching_All_alpha3_github(MatchedLocation{i,1},method);
    
    Curr_inliersRecord = false(ParNum(i),1);
    Curr_inliersRecord(inliers) = true;
    inliersParCell{i} = Curr_inliersRecord;
    inliersParIndices(i) = i;
end

for j = 1:TraverNum
   Current_j = inliersParIndices(j);
   inliersPar(RecordIndiesNew(Current_j,1):RecordIndiesNew(Current_j,2)) = inliersParCell{Current_j};
end


