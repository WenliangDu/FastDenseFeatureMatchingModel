function [inliers outliers,x]=WGTM(P1,P2,K,emu)
% Weighted Graph Transformation Matching Algorithm: 
%          The algorithm takes two sets of corresponding points (initial matches) and find matches
%          using the geometrical relation between points.

% Input -  P1: an array with size N X 2 containing the first set of points.
%       -  P2: N X 2 vactor containing the corresponding points of the first set of points (the m-th point of P1 and P2 are initially matched).
%       -  K: The number of closet neighbores of a point that are connected with point in the graph.
%       -  emu: is the stopping threshold or the maximum change in the weighted graph after removing an
%               outlier.

% Output - inliers: the inlier matches
%        - outliers: the outlier matches


% Copyright (c) 2010-2011 by Mohammad Izadi (mia4@sfu.ca), 
% Laboratory of Robotic Vision, School of Engineering Science, Simon Fraser University

% Created: May 10 2010
% Last Modified October 18 2011  
%
%% Distance matrices
dist1=zeros(size(P1,1));
dist2=zeros(size(P2,1));
for i=1:size(P1,1)
    for j=1:size(P1,1)
        dist1(i,j)=norm(P1(i,:)-P1(j,:));
        dist2(i,j)=norm(P2(i,:)-P2(j,:));
    end
end
%% median distance
medianDist1=median(dist1(:));
medianDist2=median(dist2(:));

%% WGTM algorithm

outliers=[];

% creating KNN graphs
KNNG1 = createKNNGraph(dist1,K,medianDist1);
KNNG2 = createKNNGraph(dist2,K,medianDist2);

% Creating the weighted graph
WG = createWeightedGraph(KNNG1,KNNG2,P1,P2);

angDist=90;
flag=0;
mu=[];
Mu=[];
x=[];
y=[];
for i=1:length(WG)

    if isempty(WG(i,WG(i,:)>=0))
        x(i)=-1;
        y(i)=-1;
    else
        x(i)=mean(WG(i,WG(i,:)>=0));
        y(i)=std(WG(i,WG(i,:)>=0));
    end
end
inliers=1:length(P1);
inliers=inliers(~ismember(inliers,outliers));

Mu=[Mu;max(WG(find(WG>=0)))*180/3.14];
mu=[mu;mean(x(x>=0))];
while flag==0
    KNNG1_old=KNNG1;
    KNNG2_old=KNNG2;
    
    % finding nodes that have less than two edges and considering as outlier
    x=max(WG');
    ind=find(x<0);
    x1=sum(KNNG1');
    x2=sum(KNNG1);
    ind=[ind find(x1<=2&x2<=2)];
    ind=ind(~ismember(ind,outliers));
    flag2=0;
    if ~isempty(ind)
        flag2=1;
        for i=1:length(ind)
            outliers=[outliers ind(i)];
            dist1(ind(i),:)=medianDist1+1;dist1(:,ind(i))=medianDist1+1;
            dist2(ind(i),:)=medianDist2+1;dist2(:,ind(i))=medianDist2+1;
        end
        x=3.14;
        y=1;
    else
        
        %finding the worst node as outlier
        x=[];
        y=[];
        for i=1:length(WG)
            
            if isempty(WG(i,WG(i,:)>=0))
                x(i)=-1;
                y(i)=-1;
            else
                x(i)=mean(WG(i,WG(i,:)>=0));
                y(i)=std(WG(i,WG(i,:)>=0));
            end
        end
        [a,ind]=max(x);
        outliers=[outliers ind(1)];
       
        x(ind(1))=-1;
        dist1(ind(1),:)=medianDist1+1;dist1(:,ind(1))=medianDist1+1;
        dist2(ind(1),:)=medianDist2+1;dist2(:,ind(1))=medianDist2+1;
    end
    
    % updating the KNN graphs and the weeighted graph
    KNNG1 = createKNNGraph(dist1,K,medianDist1);
    KNNG2 = createKNNGraph(dist2,K,medianDist2);
    ind=find(sum(xor(KNNG1_old,KNNG1)')+sum(xor(KNNG2_old,KNNG2)'));
    [WG] = updateWeightedGraph(KNNG1,KNNG2,P1,P2,ind,WG);
    
    inliers=find(sum(KNNG1')>1);
    inliers=1:length(P1);
    inliers=inliers(~ismember(inliers,outliers));
    if isempty(mu)
        Mu=[360 ;max(WG(find(WG>=0)))*180/3.14 ];
        mu=[2*pi ;mean(x(x>=0))];
    else
        Mu=[Mu;max(WG(find(WG>=0)))*180/3.14];
        mu=[mu;mean(x(x>=0))];
    end
    
     [length(inliers) abs(mu(end)-mu(end-1)) Mu(end) angDist]
    
    % checking the stop criteria
    if  length(inliers)<=2||(flag2==0 && Mu(end)<angDist && abs(mu(end)-mu(end-1))<emu)
        flag=1;
    end    
end

% finding inliers
inliers=1:length(P1);
inliers=inliers(~ismember(inliers,outliers));
if length(inliers)<=2
    inliers=[];
    outliers=1:length(P1);    
end
end
%%
% the function creates KNN graph that is not symetric.
function G=createKNNGraph(dist,K,medianDist)
G = zeros(length(dist));
for i=1:length(dist)
    ind=find(dist(i,:)<=medianDist);
    [d,ind1]=sort(dist(i,ind));
    ind=ind(ind1(1:min(K,length(ind1))));
    G(i,ind)=1;
end
%  G=(G+G')>0;
end
%%
% the function creates weighted graph.
function [Gs]=createWeightedGraph(G1,G2,P1,P2)
    Gs = -ones(size(P1,1));
    for i=1:size(P1,1)    
        ind=find(G1(i,:)>0);
        ind=ind(ind~=i);
        ind1=find(G2(i,ind)==0);

        if ~isempty(ind)

            t=getAngularDistance(P1(ind,:)-ones(length(ind),1)*P1(i,:),P2(ind,:)-ones(length(ind),1)*P2(i,:));

            if length(ind1)/length(ind)>0.5
                t(ind1)=3.14; 
            end
            if any(isnan(t))
                t(isnan(t))=-1;
            end
             Gs(i,ind)=t;
        end
    end
end
%%
% the function updated weighted graph fastly.
function [Gs]=updateWeightedGraph(G1,G2,P1,P2,outliers,Gs)
    for j=1:length(outliers)
        i=outliers(j);
        ind=find(G1(i,:)>0);
        ind=ind(ind~=i);
        ind1=find(G2(i,ind)==0);

        if ~isempty(ind)

            t=getAngularDistance(P1(ind,:)-ones(length(ind),1)*P1(i,:),P2(ind,:)-ones(length(ind),1)*P2(i,:));

            if length(ind1)/length(ind)>0.5
                t(ind1)=3.14; 
            end
            if any(isnan(t))
                t(1:end)=-1;
            end
            Gs(i,:)=-1;
             Gs(i,ind)=t;
        else
            Gs(i,:)=-1;
        end
    end
end
%%
% the function computes the angular distance between two nodes.
function t=getAngularDistance(r1,r2)
    t=zeros(size(r1,1));
    for j=1:size(r1,1)
        t1=acos(dot(r1(j,:),r2(j,:))/(norm(r1(j,:))*norm(r2(j,:))));
        a=r2(j,:)*[cos(t1) -sin(t1);sin(t1) cos(t1)];
        if abs(a(2)/a(1)-r1(j,2)/r1(j,1))>abs(a(1)/a(2)-r1(j,2)/r1(j,1))
            t1=-t1;
        end

        r3=r2*[cos(t1) -sin(t1);sin(t1) cos(t1)];
        for i=1:size(r1,1)
            t(j,i)=abs(acos(dot(r1(i,:),r3(i,:))/(norm(r1(i,:))*norm(r3(i,:)))));
        end
    end
    [a,ind]=min(mean(t,2));
    t=t(ind,:);
end 
