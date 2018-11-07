function [indexPairsParAll,ScoresParAll,RecordIndies] = ParFeatureMatching_ParMatching_alpha3(features1,features2,Location1Rect,Location2Rect,ThresholdD,I2,Qualified_Matched_Location1Rect,Qualified_Matched_Location2Rect)
%{
2018/09/11
ParFeatureMatching_ParMatching_alpha3
1. Get dynamic searching reange for real HiRISE images

2018/04/12
ParFeatureMatching_ParMatching_alpha2
1. Add i=1 and i = Match_Num
%}
Mdl_Location1Rect_Column = KDTreeSearcher(Location1Rect(:,2));
Mdl_Location2Rect_Column = KDTreeSearcher(Location2Rect(:,2));

Location1RectNum = size(Location1Rect,1);
Allindex_L = 1:Location1RectNum;
Allindex_L = Allindex_L';
ToltalNum_R = size(Location2Rect,1);
Allindex_R = 1:ToltalNum_R;
Allindex_R = Allindex_R';
% ThresholdD = 10;
ColumnNum = size(I2,1);
indexPairsParAll = uint32([]);
ScoresParAll = [];
RecordIndies = uint32([]);
% RecordRef = false(Match_Num,1);

%% ! here
Top = max(Mdl_Location1Rect_Column.X);%(:,2)
Bottom = min(Mdl_Location1Rect_Column.X);
TopR = max(Mdl_Location2Rect_Column.X);
BottomR = min(Mdl_Location2Rect_Column.X);

Interval = ThresholdD;
Interval_Half = Interval/2;


%%
SearchingLine_L = Top-Interval_Half:-Interval:Bottom+Interval_Half;
SearchingLine_L_Top = SearchingLine_L + Interval_Half;
SearchingLine_L_Bottom = SearchingLine_L - Interval_Half;
SearchingLine_L_Top = [SearchingLine_L_Top SearchingLine_L_Bottom(end)];
SearchingLine_L_Bottom = [SearchingLine_L_Bottom Bottom];
SearchingLine_L = [SearchingLine_L (SearchingLine_L_Top(end) - (SearchingLine_L_Top(end)-SearchingLine_L_Bottom(end))/2)];

%% Get the nearest features to searching line
SearchingLine_Num = length(SearchingLine_L);
Matched_Location_Num = size(Qualified_Matched_Location1Rect,1);
% SearchingLine_R = SearchingLine_L;

SearchingLine_Matrix = repmat(SearchingLine_L,[Matched_Location_Num 1]);
Matched_Location_L_Column_Matrix = repmat(Qualified_Matched_Location1Rect(:,2),[1 SearchingLine_Num]);
Matched_Location_R_Column_Matrix = repmat(Qualified_Matched_Location2Rect(:,2),[1 SearchingLine_Num]);
Distance_L_Matrix_ref = abs(Matched_Location_L_Column_Matrix - SearchingLine_Matrix);
Distance_L_Matrix = Matched_Location_L_Column_Matrix - SearchingLine_Matrix;
Distance_R_Matrix = Matched_Location_R_Column_Matrix - SearchingLine_Matrix;

[~,Min_Dis_Num] = min(Distance_L_Matrix_ref);
% [~,Sort_Dis_Num] = sort(Distance_L_Matrix_ref);
Distance_L_Matrix_Min = zeros(1,SearchingLine_Num);
Distance_R_Matrix_Min = zeros(1,SearchingLine_Num);
for i = 1:SearchingLine_Num
    Distance_L_Matrix_Min(i) = Distance_L_Matrix(Min_Dis_Num(i),i);
    Distance_R_Matrix_Min(i) = Distance_R_Matrix(Min_Dis_Num(i),i);
end
%%
% SearchingLine_R = SearchingLine_L;
SearchingLine_R_Top = SearchingLine_L_Top + Interval;
SearchingLine_R_Bottom = SearchingLine_L_Bottom - Interval;
Lift = (Distance_L_Matrix_Min - Distance_R_Matrix_Min);
LiftNegative = Lift;
LiftPositive = Lift;
LiftNegative(Lift > 0) = 0;
LiftPositive(Lift < 0) = 0;
SearchingLine_R_Bottom = SearchingLine_R_Bottom - LiftPositive;
SearchingLine_R_Top = SearchingLine_R_Top - LiftNegative;

Over_SearchingLine_R_Top = SearchingLine_R_Top > TopR;
SearchingLine_R_Top(Over_SearchingLine_R_Top) = TopR;
Over_SearchingLine_R_Bottom = SearchingLine_R_Bottom < BottomR;
SearchingLine_R_Bottom(Over_SearchingLine_R_Bottom) = BottomR;
%}

SearchingLine_R_Bottom(end) = SearchingLine_R_Bottom(end) - 1;
SearchingLine_L_Bottom(end) = SearchingLine_L_Bottom(end) - 1;
SearchingLine_L_Top(1) = Top;
SearchingLine_R_Top(1) = TopR;



% for i = 1:SearchingLine_Num
parfor (i = 1:SearchingLine_Num,4)
%     Curr_L_C = (1+2*ThresholdD)*i - ThresholdD;
%     Curr_L_T = Curr_L_C - ThresholdD -1;
%     Curr_L_B = Curr_L_C + ThresholdD;
%     
%     Curr_R_C = Curr_L_C;
%     Curr_R_T = Curr_R_C - 2*ThresholdD - 1;
%     Curr_R_B = Curr_R_C + 2*ThresholdD;
    
    Curr_L_idxB = (Mdl_Location1Rect_Column.X > SearchingLine_L_Bottom(i)) & (Mdl_Location1Rect_Column.X <= SearchingLine_L_Top(i));
    Curr_L_idx = Allindex_L(Curr_L_idxB);
    
    Curr_R_idxB = (Mdl_Location2Rect_Column.X > SearchingLine_R_Bottom(i)) & (Mdl_Location2Rect_Column.X <= SearchingLine_R_Top(i));
    Curr_R_idx = Allindex_R(Curr_R_idxB);
    
    if ~isempty(Curr_L_idx)
        Curr_features_L = features1(:,Curr_L_idx);
        Curr_features_R = features2(:,Curr_R_idx);
        [Curr_indexPairs_Temp, Curr_Scores] = vl_ubcmatch(Curr_features_L,Curr_features_R);
        if ~isempty(Curr_indexPairs_Temp)
%             if i==190
%                 x = 1; 
%             end
            Curr_indexPairs = [Curr_L_idx(Curr_indexPairs_Temp(1,:)) Curr_R_idx(Curr_indexPairs_Temp(2,:))];
            Curr_RecordIndies = ones(size(Curr_indexPairs,1),1);
            Curr_RecordIndies = Curr_RecordIndies*i;
            indexPairsParAll = [indexPairsParAll;Curr_indexPairs];
            ScoresParAll = [ScoresParAll Curr_Scores];
            RecordIndies = [RecordIndies;Curr_RecordIndies];
            
            %% For show
%             Left = min(Location1Rect(:,1));
%             Right = max(Location1Rect(:,1));
%             Curr_MatchedLocation1 = Location1Rect(Curr_indexPairs(:,1),:);
%             Curr_MatchedLocation2 = Location2Rect(Curr_indexPairs(:,2),:);
%             
%             figure,
%             line([Left;Right],[SearchingLine_L_Top(i);SearchingLine_L_Top(i)],'Color','r');hold on
%             line([Left;Right],[SearchingLine_L_Bottom(i);SearchingLine_L_Bottom(i)],'Color','k');
%             
%             line([Left;Right],[SearchingLine_R_Top(i);SearchingLine_R_Top(i)],'Color','m');
%             line([Left;Right],[SearchingLine_R_Bottom(i);SearchingLine_R_Bottom(i)],'Color','b');
%             
%             plot(Location1Rect(Curr_L_idx,1),Location1Rect(Curr_L_idx,2),'r*');
%             plot(Location2Rect(Curr_R_idx,1),Location2Rect(Curr_R_idx,2),'k*');
%             
%             plot(Curr_MatchedLocation1(:,1),Curr_MatchedLocation1(:,2),'mo','MarkerSize',12);
%             plot(Curr_MatchedLocation2(:,1),Curr_MatchedLocation2(:,2),'bo','MarkerSize',12);
%             line([Curr_MatchedLocation1(:,1)';Curr_MatchedLocation2(:,1)'],[Curr_MatchedLocation1(:,2)';Curr_MatchedLocation2(:,2)'],'Color','g');
%             x = 1;
%             RecordRef(i) = true;
        end
    end
end
%{
%% when i == 1
Curr_L_C = (1+2*ThresholdD) - ThresholdD;
Curr_L_T = Curr_L_C - ThresholdD - 1;
Curr_L_B = Curr_L_C + ThresholdD;

Curr_R_C = Curr_L_C;
Curr_R_T = Curr_L_T;
Curr_R_B = Curr_R_C + 2*ThresholdD;

Curr_L_idxB = (Mdl_Location1Rect.X(:,2) <= Curr_L_B) & (Mdl_Location1Rect.X(:,2) >= Curr_L_T);
Curr_L_idx = Allindex_L(Curr_L_idxB);

Curr_R_idxB = (Mdl_Location2Rect.X(:,2) <= Curr_R_B) & (Mdl_Location2Rect.X(:,2) > Curr_R_T);
Curr_R_idx = Allindex_R(Curr_R_idxB);

if ~isempty(Curr_L_idx)
    Curr_features_L = features1(:,Curr_L_idx);
    Curr_features_R = features2(:,Curr_R_idx);
    [Curr_indexPairs_Temp, Curr_Scores] = vl_ubcmatch(Curr_features_L,Curr_features_R);
    if ~isempty(Curr_indexPairs_Temp)
        Curr_indexPairs = [Curr_L_idx(Curr_indexPairs_Temp(1,:)) Curr_R_idx(Curr_indexPairs_Temp(2,:))];
        Curr_RecordIndies = ones(size(Curr_indexPairs,1),1);
        Curr_RecordIndies = Curr_RecordIndies*1;
        indexPairsParAll = [Curr_indexPairs;indexPairsParAll];
        ScoresParAll = [Curr_Scores ScoresParAll];
        RecordIndies = [Curr_RecordIndies;RecordIndies];
%         RecordRef(1) = true;
    end
end
%% when i == Match_Num
Curr_L_C = (1+2*ThresholdD)*Match_Num - ThresholdD;
Curr_L_T = Curr_L_C - ThresholdD - 1;
Curr_L_B = ColumnNum;

Curr_R_C = Curr_L_C;
Curr_R_T = Curr_R_C - 2*ThresholdD - 1;
Curr_R_B = ColumnNum;

Curr_L_idxB = (Mdl_Location1Rect.X(:,2) <= Curr_L_B) & (Mdl_Location1Rect.X(:,2) > Curr_L_T);
Curr_L_idx = Allindex_L(Curr_L_idxB);

Curr_R_idxB = (Mdl_Location2Rect.X(:,2) <= Curr_R_B) & (Mdl_Location2Rect.X(:,2) > Curr_R_T);
Curr_R_idx = Allindex_R(Curr_R_idxB);

if ~isempty(Curr_L_idx)
    Curr_features_L = features1(:,Curr_L_idx);
    Curr_features_R = features2(:,Curr_R_idx);
    [Curr_indexPairs_Temp, Curr_Scores] = vl_ubcmatch(Curr_features_L,Curr_features_R);
    if ~isempty(Curr_indexPairs_Temp)
        Curr_indexPairs = [Curr_L_idx(Curr_indexPairs_Temp(1,:)) Curr_R_idx(Curr_indexPairs_Temp(2,:))];
        Curr_RecordIndies = ones(size(Curr_indexPairs,1),1);
        Curr_RecordIndies = Curr_RecordIndies*Match_Num;
        indexPairsParAll = [indexPairsParAll;Curr_indexPairs];
        ScoresParAll = [ScoresParAll Curr_Scores];
        RecordIndies = [RecordIndies;Curr_RecordIndies];
%         RecordRef(Match_Num) = true;
    end
end
%}

%{
k = uint32(1);
RecordIndiesRef = RecordIndies;
RecordIndiesNew = [];
while ~isempty(RecordIndiesRef)
    RecordIndiesRef(1)
end
            %% ParFeatureMatching_ParMatching_alpha2
            Curr_indexPairsNum = size(Curr_indexPairs,1);
            Curr_RecordIndies(1) = k;
            Curr_RecordIndies(2) = k + Curr_indexPairsNum - 1;
            RecordIndies = [RecordIndies;Curr_RecordIndies];
            k = k + Curr_indexPairsNum;
%}