function [indexPairs,Scores] = SIFT_RemoveSameMatching_alpha1(MatchedLocation1,MatchedLocation2,Scores,indexPairs)
%{
2016/03/25
SIFT_RemoveSameMatching_alpha1
1. For same matching caused by SIFT

%}

ToBeRemovedNodesSet1 = SIFT_RemoveSameMatching_Sub_alpha1(MatchedLocation1,Scores);


indexPairs(ToBeRemovedNodesSet1,:) = [];
MatchedLocation2(ToBeRemovedNodesSet1,:) = [];
Scores(ToBeRemovedNodesSet1) = [];

ToBeRemovedNodesSet2 = SIFT_RemoveSameMatching_Sub_alpha1(MatchedLocation2,Scores);
indexPairs(ToBeRemovedNodesSet2,:) = [];
Scores(ToBeRemovedNodesSet2) = [];