function [inliers,Inter]= SOCBV(P1,P2,K,ea)
% Weighted Graph Transformation Matching Algorithm: 
%          The algorithm takes two sets of corresponding points (initial matches) and find matches
%          using the geometrical relation between points.

% Input -  P1: an array with size N X 2 containing the first set of points.
%       -  P2: N X 2 vactor containing the corresponding points of the first set of points (the m-th point of P1 and P2 are initially matched).
%       -  K: The number of closet neighbores of a point that are connected with point in the graph.
%       -  ea: is the stopping threshold or the min value of remaining matches in the weighted graph after removing an
%               outlier.

% Output - inliers: the inlier matches
%        - outliers: the outlier matches


% Copyright (c) 2014-2015 by Fanyang Meng (mengfanyang@email.szu.edu.cn)
% the college of Information Engineering, Shenzhen University, Shenzhen, China.
% we  would like to thank M.Izadi for his help with codes of GTM and WGTM
% algorithm.

%
%% Distance matrices
dist1=zeros(size(P1,1));angle1 = zeros(size(P1,1));
dist2=zeros(size(P2,1));angle2 = zeros(size(P2,1));
for i=1:size(P1,1)
    for j=1:size(P1,1)
        dist1(i,j)=norm(P1(i,:)-P1(j,:));
        dist2(i,j)=norm(P2(i,:)-P2(j,:));
        angle1(i,j) = atan2(P1(i,1)-P1(j,1),P1(i,2)-P1(j,2));  
        angle2(i,j) = atan2(P2(i,1)-P2(j,1),P2(i,2)-P2(j,2));  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
medianDist1 = median(dist1(:));medianDist2 = median(dist2(:));
% Creating the weighted
W = exp(-abs(dist1/medianDist1 - dist2/medianDist2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOCBV algorithm
% creating KNN graphs
    KNNG1 = createKNNGraph(dist1,K,medianDist1,W);
    KNNG2 = createKNNGraph(dist2,K,medianDist2,W);
% Creating the vote graph
    [VGX,VGY] = createVoteGraph(KNNG1,KNNG2,angle1,angle2,W);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Inneighbor and Outneighbor with weighted
     NeighborVoteX = sum((VGX'+VGX))./(sum((KNNG1+KNNG1'))+eps);
     NeighborVoteY = sum((VGY'+VGY))./(sum((KNNG2+KNNG2'))+eps);
     NeighborVote  = (NeighborVoteX + NeighborVoteY)/2;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flag=0;Inter = 1;
    inliers =1:length(P1); outliers = [];
while flag==0
     KNNG1_old=KNNG1;    KNNG2_old=KNNG2;    %finding the worst node as outlier
     [~,ind] = min(NeighborVote); 
     outliers = [outliers ind(1)];  
     dist1(ind(1),:)= medianDist1 + 1;dist1(:,ind(1))= medianDist1 + 1;
     dist2(ind(1),:)= medianDist2 + 1;dist2(:,ind(1))= medianDist2 + 1;
    % updating the KNN graphs and the InVote graph    
     KNNG1 = createKNNGraph(dist1,K,medianDist1,W);
     KNNG2 = createKNNGraph(dist2,K,medianDist2,W); 
     ind = find(sum(xor(KNNG1_old,KNNG1),2)+sum(xor(KNNG2_old,KNNG2),2))';      
     [VGX,VGY] = updateVoteGraph(KNNG1,KNNG2,angle1,angle2,ind,VGX,VGY,W);
     
     VGX(outliers(end),:) = 0;        VGX(:,outliers(end)) = 0;
     VGY(outliers(end),:) = 0;        VGY(:,outliers(end)) = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % updating Inneighbor and Outneighbor with weighted
     NeighborVoteX = sum((VGX'+VGX))./(sum((KNNG1+KNNG1'))+eps);
     NeighborVoteY = sum((VGY'+VGY))./(sum((KNNG2+KNNG2'))+eps);
     NeighborVote  = (NeighborVoteX + NeighborVoteY)/2;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inliers  =1:length(P1);
    inliers  = inliers(~ismember(inliers,outliers));  
    NeighborVote(outliers) = Inf;
    EA       =  1-  min(NeighborVote(inliers));  
     % checking the stop criteria
       if  length(inliers)<= 3||EA< ea 
           flag = 1;
       end 
       Inter = Inter +1;  
end

end
%%
% the function creates KNN graph that is not symetric.
function G=createKNNGraph(dist,K,medianDist,W)
    G = zeros(length(dist));
    for i = 1:length(dist)
    ind = find(dist(i,:)<medianDist);  
       if length(ind)>3
            [~,ind1] = sort(dist(i,ind)); 
            ind = ind(ind1(1:min(K,length(ind1))));  
            G(i,ind) = W(i,ind); 
       end
            G(i,i) = 1; W(i,i) =1; 
    end
end
%%
% the function creates weighted graph.
function [VGX,VGY]= createVoteGraph(G1,G2,A1,A2,W)
    VGX = zeros(size(A1));    VGY = zeros(size(A2));
    for i=1:size(A1,1)  
        % VGX
        ind = find(G1(i,:)>0);        ind = ind(ind~=i);
        if ~isempty(ind)
            knnA1X = A1(i,ind);knnA2X = A2(i,ind);
            [~,IndxAX] = sort(knnA1X);     [~,IndxAY] = sort(knnA2X);
            CirLcs = CirclLCS(IndxAX,IndxAY); 
            VGX(i,[ind(CirLcs) i])= length(CirLcs)/length(ind)*W(i,[ind(CirLcs) i]);    
        else
             VGX(i,[i ind])= 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % VGY
            ind = find(G2(i,:)>0);        ind = ind(ind~=i);           
            if ~isempty(ind)
                knnA1Y = A1(i,ind);knnA2Y = A2(i,ind);
                [~,IndxAX] = sort(knnA1Y);     [~,IndxAY] = sort(knnA2Y);
                CirLcs = CirclLCS(IndxAY,IndxAX);  
                VGY(i,[ind(CirLcs) i])= length(CirLcs)/length(ind)*W(i,[ind(CirLcs) i]);     
            else
               VGY(i,[i ind])= 0; 
            end
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the function updated weighted graph fastly.
function [VGX,VGY] = updateVoteGraph(G1,G2,A1,A2,outliers,VGX,VGY,W)
    for j=1:length(outliers)
        i = outliers(j);   
        ind = find(G1(i,:)>0);ind=ind(ind~=i); 
        if ~isempty(ind)  
            VGX(i,[i ind])= 0; 
            knnA1X = A1(i,ind);knnA2X = A2(i,ind);
            [~,IndxAX] = sort(knnA1X);     [~,IndxAY] = sort(knnA2X);
            CirLcs = CirclLCS(IndxAX,IndxAY);  
            VGX(i,[ind(CirLcs) i])= length(CirLcs)/length(ind)*W(i,[ind(CirLcs) i]); 
        else
             VGX(i,[i ind])= 0; 
       end   
        %%%%%%%%%%%%%%%%%%%%%%%%
         ind = find(G2(i,:)>0);     ind = ind(ind~=i);        
         if ~isempty(ind)
             VGY(i,[i ind])= 0; 
             knnA1Y = A1(i,ind);knnA2Y = A2(i,ind);
             [~,IndxAX] = sort(knnA1Y);     [~,IndxAY] = sort(knnA2Y);       
             CirLcs = CirclLCS(IndxAY,IndxAX);   
             VGY(i,[i ind(CirLcs)])= length(CirLcs)/length(ind)*W(i,[ind(CirLcs) i]);
          else
             VGY(i,[i ind])= 0; 
         end             
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lcstr,B] = CirclLCS(a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Db = [b b]; 
    M = length(a);   N = length(Db);
    C = zeros(M+1,N+1);  B = C;
    for i = 2:M+1  
        for j = 2:N+1 
              if a(i-1) == Db(j-1)             
                         C(i,j) = C(i-1,j-1) + 1;                 
                         B(i,j) = 0;  
              else
                      if C(i-1,j) >= C(i,j-1)
                             C(i,j) = C(i-1,j);
                             B(i,j)  = 1;
                    else         
                         C(i,j) = C(i,j-1);
                         B(i,j) = -1;                   
                      end          
              end   
        end    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    LCSL = 0;  LIndx= [];
    for i = 1: (ceil(N/2)-1)
          TM = M+1;   TN = ceil(N/2)+i;
          Indx = zeros(1,M);         
           while TM>1 && TN>i   
                  if B(TM,TN) == 0  
                      Indx(TM-1) = 1;  
                      TN  = TN-1; TM = TM-1; 
                  else
                       if B(TM,TN)==-1  
                                TN = TN-1;        
                       else
                               TM = TM-1;
                       end    
                  end          
           end
            if LCSL<sum(Indx)
                LCSL = sum(Indx);
                LIndx = Indx;
            end
    end  
    lcstr = a(LIndx>0);
end
