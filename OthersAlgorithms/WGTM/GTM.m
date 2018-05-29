function [inliers outliers]=GTM(P1,P2,K)
% Graph Transformation Matching Algorithm: 
%          The algorithm takes two sets of corresponding points (initial matches) and find matches
%          by finding two isomorphic graphs for them. 
% Input -  P1: N X 2 vactor containing the first set of points (it could be generally used for matching points of higher dimentional spaces).
%       -  P2: N X 2 vactor containing the corresponding points of the first set of points (the m-th point of P1 and P2 are initially matched).
%       -  K: the number of closet neighbores of a point that are connected with point in the graph.

% Output - inliers: the inlier matches
%        - outliers: the outlier matches


% Copyright (c) 2010-2011 by Mohammad Izadi (mia4@sfu.ca), 
% Laboratory of Robotic Vision, School of Engineering Science, Simon Fraser University

% Created: April 25 2010
% Last Modified July 04 2011  

%% Distance matrices
dist1=zeros(size(P1,1));
dist2=zeros(size(P2,1));
for i=1:size(P1,1)
    for j=1:size(P1,1)
        dist1(i,j)=sqrt(sum((P1(i,:)-P1(j,:)).^2));
        dist2(i,j)=sqrt(sum((P2(i,:)-P2(j,:)).^2));
    end
end
%% median distance
medianDist1=median(dist1(:));
medianDist2=median(dist2(:));

%% GTM
outliers=zeros(1,length(P1));
G1 = createKNNGraph(dist1,K,medianDist1);
G2 = createKNNGraph(dist2,K,medianDist2);
h = waitbar(0,'Please wait...');
while sum(sum(G1~=G2))>0
    waitbar(sum(outliers)/size(P1,1))
    R = abs(G1-G2);
    [m,ind]=max(sum(R'));
     outliers(ind(1))=1;
    G1(ind(1),:)=0;G1(:,ind(1))=0;
    G2(ind(1),:)=0;G2(:,ind(1))=0;
    dist1(ind(1),:)=medianDist1+1;dist1(:,ind(1))=medianDist1+1;
    dist2(ind(1),:)=medianDist2+1;dist2(:,ind(1))=medianDist2+1;
    G1 = createKNNGraph(dist1,K,medianDist1);
    G2 = createKNNGraph(dist2,K,medianDist2);
end
outliers=find(sum(G1,2)'<=1);
inliers=find(sum(G1,2)'>1);
close(h)
end
% function G=createKNNGraph(dist,K,medianDist)
% G = zeros(length(dist));
% for i=1:length(dist)
%     ind=find(dist(i,:)<=medianDist);
%     [d,ind1]=sort(dist(i,ind));
%     ind=ind(ind1(1:min(K,length(ind1))));
%     G(i,ind)=1;
% end
% end
% 
% 
%         
%     
function G=createKNNGraph(dist,K,medianDist)
G = zeros(length(dist));

% disconnecting all vertices that do not have K neighbors
flag=0;
while flag==0
    flag=1;
    for i=1:length(dist)
        ind=find(dist(i,:)<=medianDist);
        if length(ind)<K&& length(ind)>0
            dist(i,:)=Inf;dist(:,i)=Inf;
            flag=0;
        end
    end    
end

% creating K-NN graph
for i=1:length(dist)
    ind=find(dist(i,:)<=medianDist);
    [d,ind1]=sort(dist(i,ind));
    if length(ind)>=K
        ind=ind(ind1(1:K));
    else
        ind=[];
    end
    G(i,ind)=1;
end
end
