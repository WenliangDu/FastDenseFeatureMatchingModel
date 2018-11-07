function inliersPar = ParFeatureMatching_RemoveOutliers_alpha2_1(Location1,Location2,indexPairsPar,RecordIndies,method)

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
ThresholdNumber = 300;
% RecordIndiesNum = size(RecordIndies,1);
ParInWrong = (ParInNumber) < ThresholdNumber;
% IfallWrong = false;
% ParInRight = ~ParInWrong;
% RecordIndiesNew = RecordIndies;
%{
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
%}

%%



% RecordIndiesNew = RecordIndies;
% ParInNumberNew = ParInNumber;
if any(ParInWrong)
%     ParInWrongNumber = find(ParInWrong);
    while any(ParInWrong)
        
        Curr_Num = find(ParInWrong,1);
        Curr_Toltal_Num = ParInNumber(Curr_Num);
        CaseNum = false;
        if Curr_Num+1 < length(ParInNumber)
            for i = Curr_Num+1:length(ParInNumber)
                Curr_Toltal_Num = Curr_Toltal_Num + ParInNumber(i);
                if Curr_Toltal_Num >= ThresholdNumber
                    Curr_Next = i;
                    CaseNum = true;
                    break
                end
    %             if i == length(ParInNumber)
    %                 
    %             end
            end
        end
        
        if CaseNum
%             RecordIndiesNew(Curr_Num,2) = RecordIndies(Curr_Next,2);
%             RecordIndiesNew(Curr_Num+1:Curr_Next,:) = [];
%             ParInNumberNew(Curr_Num) = Curr_Toltal_Num;
%             ParInNumberNew(Curr_Num+1:Curr_Next) = [];
            ParInWrong(Curr_Num:Curr_Next) = false;
            ParInWrong(Curr_Num+1:Curr_Next) = [];
            
            ParInNumber(Curr_Num) = Curr_Toltal_Num;
            RecordIndies(Curr_Num,2) = RecordIndies(Curr_Next,2);
            RecordIndies(Curr_Num+1:Curr_Next,:) = [];
            ParInNumber(Curr_Num+1:Curr_Next) = [];
            
            x = 1;
        else
            x =1;
            RecordIndies(Curr_Num-1,2) = RecordIndies(end,2);
            ParInNumber(Curr_Num-1) = ParInNumber(Curr_Num-1) + Curr_Toltal_Num;
            RecordIndies(Curr_Num:end,:) = [];
            ParInNumber(Curr_Num:end) = [];
            
            ParInWrong(Curr_Num:end) = false;
            ParInWrong(Curr_Num:end) = [];
        end
        x = 1;
        
        
        
    end
    
    
    
end
RecordIndiesNew = RecordIndies;
ParInNumberNew = ParInNumber;



%%
%inliersPar = logical([]);
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
%     MatchedLocation{i,2} = Location2(Curr_indexPairs(:,2),:);
end
% for i = 1:TraverNum
parfor (i = 1:TraverNum,4)
%     Curr_indexPairs = indexPairsPar(RecordIndiesNew(i,1):RecordIndiesNew(i,2),:);
%     MatchedLocation1 = Location1(Curr_indexPairs(:,1),:);
%     MatchedLocation2 = Location2(Curr_indexPairs(:,2),:);
    
%     figure,plot(MatchedLocation1(:,1),MatchedLocation1(:,2),'r*');hold on
%     plot(MatchedLocation2(:,1),MatchedLocation2(:,2),'b+');
%     line([MatchedLocation1(:,1)';MatchedLocation2(:,1)'],[MatchedLocation1(:,2)';MatchedLocation2(:,2)'],'Color','y');
%     x = 1;
    
    inliers = Other_Matching_All_alpha3(MatchedLocation{i,1},method);
    
%     figure,plot(MatchedLocation1(inliers,1),MatchedLocation1(inliers,2),'r*');hold on
%     plot(MatchedLocation2(inliers,1),MatchedLocation2(inliers,2),'b+');
%     line([MatchedLocation1(inliers,1)';MatchedLocation2(inliers,1)'],[MatchedLocation1(inliers,2)';MatchedLocation2(inliers,2)'],'Color','b');
    
%     x = 1;
    Curr_inliersRecord = false(ParNum(i),1);
    Curr_inliersRecord(inliers) = true;
    inliersParCell{i} = Curr_inliersRecord;
    inliersParIndices(i) = i;
end

for j = 1:TraverNum
   Current_j = inliersParIndices(j);
   inliersPar(RecordIndiesNew(Current_j,1):RecordIndiesNew(Current_j,2)) = inliersParCell{Current_j};
end

%{
for i = 1:TraverNum
% parfor (i = 1:TraverNum,4)
    Curr_indexPairs = indexPairsPar(RecordIndiesNew(i,1):RecordIndiesNew(i,2),:);
    inliers = Other_Matching_All_alpha1(Location1,Location2,Curr_indexPairs,method);
    
    Curr_inliersRecord = false(RecordIndiesNew(i,2)-RecordIndiesNew(i,1)+1,1);
    Curr_inliersRecord(inliers) = true;
    
    inliersParCell{i} = Curr_inliersRecord;
end
for j = 1:TraverNum
   inliersPar(RecordIndiesNew(j,1):RecordIndiesNew(j,2)) = inliersParCell{j};
end
%}
%

% for i = 1:TraverNum
% % parfor (i = 1:TraverNum,4)
% 
%     Curr_indexPairs = indexPairsPar(RecordIndiesNew(i,1):RecordIndiesNew(i,2),:);
%     inliers = Other_Matching_All_alpha1(Location1,Location2,Curr_indexPairs,method);
%     
%     Curr_inliersRecord = false(RecordIndiesNew(i,2)-RecordIndiesNew(i,1)+1,1);
%     Curr_inliersRecord(inliers) = true;
%     inliersPar(RecordIndiesNew(i,1):RecordIndiesNew(i,2)) = Curr_inliersRecord;
% end
%}
% if ~isempty(inliersPar)
%     inliersAll = false(size(Location1,1),1);
%     inliersAll(indexPairsPar(inliersPar,1)) = true;
% end

