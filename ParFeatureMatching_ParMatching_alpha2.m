function [indexPairsParAll,ScoresParAll,RecordIndies] = ParFeatureMatching_ParMatching_alpha2(features1,features2,Location1Rect,Location2Rect,ThresholdD,I2)
%{
2018/04/12
ParFeatureMatching_ParMatching_alpha2
1. Add i=1 and i = Match_Num

%}
Mdl_Location1Rect = KDTreeSearcher(Location1Rect);
Mdl_Location2Rect = KDTreeSearcher(Location2Rect);

Location1RectNum = size(Location1Rect,1);
Allindex_L = 1:Location1RectNum;
Allindex_L = Allindex_L';
ToltalNum_R = size(Location2Rect,1);
Allindex_R = 1:ToltalNum_R;
Allindex_R = Allindex_R';
% ThresholdD = 10;
ColumnNum = size(I2,1);
Match_Num = floor(ColumnNum/(1+2*ThresholdD));
indexPairsParAll = [];
ScoresParAll = [];
RecordIndies = [];
% RecordRef = false(Match_Num,1);

% for i = 2:Match_Num-1
parfor (i = 2:Match_Num-1,4)
    Curr_L_C = (1+2*ThresholdD)*i - ThresholdD;
    Curr_L_T = Curr_L_C - ThresholdD -1;
    Curr_L_B = Curr_L_C + ThresholdD;
    
    Curr_R_C = Curr_L_C;
    Curr_R_T = Curr_R_C - 2*ThresholdD - 1;
    Curr_R_B = Curr_R_C + 2*ThresholdD;
    
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
            Curr_RecordIndies = Curr_RecordIndies*i;
            indexPairsParAll = [indexPairsParAll;Curr_indexPairs];
            ScoresParAll = [ScoresParAll Curr_Scores];
            RecordIndies = [RecordIndies;Curr_RecordIndies];
%             RecordRef(i) = true;
        end
    end
end
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