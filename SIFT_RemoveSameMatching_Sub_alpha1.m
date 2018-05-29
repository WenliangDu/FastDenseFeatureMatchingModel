function ToBeRemovedNodesSet = SIFT_RemoveSameMatching_Sub_alpha1(MatchedLocation,Scores)
%{
2016/03/25
SIFT_RemoveSameMatching_Sub_alpha1
1. For same matching caused by SIFT
2. The subprogram of SIFT_RemoveSameMatching_alpha1

%}




[MatchedLocationU,LocationU2O,LocationO2U] = unique(MatchedLocation,'rows','stable');

MatchedLocationL = size(MatchedLocation,1);
RecordNum = false(MatchedLocationL,1);
RecordNum(LocationU2O) = true;
EqualMatchedLocationNodes = find(RecordNum == false); % Original Order

EqualNumbers = LocationU2O( LocationO2U(EqualMatchedLocationNodes) ); % Original Order

ToBeRemovedNodesSet = [];
while ~isempty(EqualNumbers),
    SameEqualNumbers = EqualNumbers == EqualNumbers(1);
    ToBeRemovedNodes = [EqualNumbers(1);EqualMatchedLocationNodes( SameEqualNumbers )];
    [~,RightNodesOrder] = min( Scores(ToBeRemovedNodes) );
    ToBeRemovedNodes(RightNodesOrder) =[];
    
    ToBeRemovedNodesSet = [ToBeRemovedNodesSet;ToBeRemovedNodes];
    EqualNumbers( SameEqualNumbers ) = [];
    EqualMatchedLocationNodes( SameEqualNumbers ) = [];
end