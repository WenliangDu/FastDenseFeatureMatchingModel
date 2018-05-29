function [Inliers,Inter] = RSOC_OP(X,Y,K)
% Restricted spatial order constraints Algorithm: 
%          The algorithm takes two sets of corresponding points (initial matches) and find matches
%          using the geometrical relation between points.

% Input -  P1: an array with size N X 2 containing the first set of points.
%       -  P2: N X 2 vactor containing the corresponding points of the first set of points (the m-th point of P1 and P2 are initially matched).
%       -  K: The number of closet neighbores of a point that are connected with point in the graph.
%       -  ea: is the stopping threshold or the mean value of remaining matches in the weighted graph after removing an
%               outlier.

% Output - inliers:  the indx of inlier matches
%        - Inter  :  the number of iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off all;
%
 Dthres = K/2;   
 Ala = 1;     Beta = 5;    
 E_maxExchange =Inf;    E_cost = Inf;
Inliers   = 1:length(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算相关的距离序关系和角度序关系矩阵
len           =    size(X,1);
 DisMatX   = zeros(len,len);        DisMatY   = zeros(len,len); 
 AngleMatX   = zeros(len,len);  AngleMatY   = zeros(len,len); 
%% Distance and angle matrices
    for i=1:size(X,1)
        for j=1:size(X,1)
            DisMatX(i,j)   = norm(X(i,:)-X(j,:));
            DisMatY(i,j)   = norm(Y(i,:)-Y(j,:));
            AngleMatX(i,j) = atan2(X(i,1)-X(j,1),X(i,2)-X(j,2));  
            AngleMatY(i,j) = atan2(Y(i,1)-Y(j,1),Y(i,2)-Y(j,2));  
        end
    end
   KNN1 = createKNNGraph(DisMatX,K); 
   KNN2 = createKNNGraph(DisMatX,K);
 %%%%%%%%%%%%%%%%%%%
    NewX = X;   NewY = Y;    
    Inter = 1;E_Exchange = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while ((E_maxExchange(1)>Ala) || (E_cost>Beta))
         KNN1_old = KNN1;    KNN2_old=KNN2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [Dif1,Dif2,C_inlier,C_outlier]= GetCandidate(KNN1,KNN2,AngleMatX,AngleMatY,K);
     %先判断是否找到外点，如果没有则跳出循环 
            if isempty(C_outlier) 
                break;
            end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if sum(Dif1(C_outlier)<Dthres) || sum(Dif2(C_outlier)<Dthres)
                 %
               [EpostOutlier,E_Exchange]  = Estimate_eChange(NewX,NewY,C_inlier,C_outlier);                 
               [~,Num] = min(E_Exchange); 
               outlier = C_outlier(Num(1));                 
            else    
                Temp = [Dif1(C_outlier);Dif2(C_outlier)]; 
               [~,MaxIndx] = find(Temp==max(Temp(:)));
                outlier = C_outlier(MaxIndx);
            end 
            
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          Indx = logical(ones(length(NewX),1)>0);  
          Indx(outlier) = 0;           TTT = logical(Indx>0);
          Inliers = Inliers(TTT);  
        if length(Inliers)<K || isempty(E_Exchange)
             break;
        end
    % updating the KNN graphs and the InVote graph  
     DisMatX(outlier,:)= -1;DisMatX(:,outlier)= -1;
     DisMatY(outlier,:)= -1;DisMatY(:,outlier)= -1;
     KNN1 = createKNNGraph(DisMatX,K);
     KNN2 = createKNNGraph(DisMatY,K); 
     ind = find(sum(xor(KNN1_old,KNN1),2)+sum(xor(KNN2_old,KNN2),2))';
     KNN1 = updateGraph(KNN1,DisMatX,ind,K);
     KNN2 = updateGraph(KNN2,DisMatY,ind,K);
     
    DisMatX  = DisMatX(Indx,Indx);     
    DisMatY  = DisMatY(Indx,Indx);
    AngleMatX = AngleMatX(Indx,Indx);     
    AngleMatY = AngleMatY(Indx,Indx);
    KNN1 = KNN1(Indx,Indx);
    KNN2 = KNN2(Indx,Indx); 
    NewX = NewX(Indx,:);NewY =  NewY(Indx,:);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            E_maxExchange = max(abs(E_Exchange));  
            E_cost = max(EpostOutlier); 
            Inter = Inter+1;
     end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    [EpostOutlier,E_Exchange] = Estimate_eChange(NewX,NewY,C_inlier,C_outlier)
      Nr = length(NewX);
      Tform    = cp2tform(NewX,NewY+eps,'affine');
      TformX   =   tformfwd(Tform,NewX) ;
%       Epre      =   norm(TformX-NewY);
      Epre      =   sqrt(sum(sum((TformX-NewY).^2,2))/Nr);
% Epre = sum(sum(abs(TformX-NewY)))/Nr;
  %计算误差EpostOutlieri  
      EpostOutlier = zeros(1,length(C_outlier));      
        for i = 1:length(C_outlier)  
            Indx = [1:i-1 i+1:length(C_outlier)]';
            LIndx  = [C_inlier ;C_outlier(Indx)];
            TempX_Outi = NewX(LIndx,:);  TempY_Outi = NewY(LIndx,:); 
%          TempX_Outi = NewX;   TempX_Outi(C_outlier(i),:) = [0 0];   
%          TempY_Outi = NewY;    TempY_Outi(C_outlier(i),:) = [0 0];   
%        TempX_Outi(C_outlier(i),:) = [0 0];   TempY_Outi(C_outlier(i),:) = [0 0];   
           Tform     = cp2tform(TempX_Outi,TempY_Outi+eps,'affine');
           TformX  =   tformfwd(Tform,TempX_Outi) ;    
%            EpostOutlier(i) = norm(TformX-TempY_Outi);
        EpostOutlier(i) = sqrt(sum(sum((TformX-TempY_Outi).^2,2))/(Nr-1));
%        EpostOutlier(i) = sum(sum(abs(TformX-TempY_Outi)))/(Nr-1);
        end
         %计算误差E_Exchange
       E_Exchange = EpostOutlier - Epre;
       %计算对应的阈值
%        E_maxExchange = max(abs(E_Exchange));
%        E_cost = max(EpostOutlier);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dif1,Dif2,C_inlier,C_outlier] = GetCandidate(KNN1,KNN2,AngleMatX,AngleMatY,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ddif = K/2;
    len           =    size(KNN1,1); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dif1 = zeros(1,len);  Dif2 = zeros(1,len);
    for i = 1:len    
        %对于第i个特征点    
        %先计算X的K紧邻约束下的空间序关系描述子对
         XKNN   = find(KNN1(i,:)>0); 
          %对应的角度
         AX_KNN = AngleMatX(i,XKNN);       AY_KNN = AngleMatY(i,XKNN);     
         [~,IndxAX] = sort(AX_KNN);     [~,IndxAY] = sort(AY_KNN);
         %对角度进行排序，并计算出对应的最短编辑距离
         Dif1(i) = EditDistance(IndxAX,[IndxAY IndxAY]) - K;
          %先计算X的K紧邻约束下的空间序关系描述子对
         YKNN   = find(KNN2(i,:)>0); 
          %对应的角度
         AX_KNN = AngleMatX(i,YKNN);       AY_KNN = AngleMatY(i,YKNN);
         %对角度进行排序，并计算出对应的最长公共子序列
         [~,IndxAX] = sort(AX_KNN);     [~,IndxAY] = sort(AY_KNN);
         Dif2(i)  = EditDistance(IndxAY,[IndxAX IndxAX]) - K;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    CandidateIndx  = ones(len,1);
%      MaxDif = zeros(len,1);
      for i = 1:len    
              if Dif1(i)==0 ||  Dif2(i)==0
                     CandidateIndx(i) = 1;
              else
                     if abs(Dif1(i)-Dif2(i))< Ddif
                         CandidateIndx(i) = 0;
                     end         
              end
%                MaxDif(i) = max(Dif1(i),Dif2(i));
      end
%        MaxIndx = find(MaxDif==max(MaxDif));
    %%%%%%%%%%%%%%%%%%%
    if prod(CandidateIndx)
           MaxDif1 = max(Dif1);        MaxDif2 = max(Dif2);   
        if MaxDif1==0 && MaxDif2==0
            CandidateIndx =  zeros(len,1);
        else
%             MaxDif = zeros(len,1);
            for i = 1:len
                  if Dif1(i)==MaxDif1 && Dif2(i)==MaxDif2
                         CandidateIndx(i) = 0;
                  end
%                    MaxDif(i) = max(Dif1(i),Dif2(i));
            end  
%             MaxIndx = find(MaxDif==max(MaxDif));
        end
    end
    C_inlier     =  find(CandidateIndx>0);
    C_outlier    =  find(CandidateIndx<=0);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function [V,v] = EditDistance(string1,string2)
% Edit Distance is a standard Dynamic Programming problem. Given two strings s1 and s2, the edit distance between s1 and s2 is the minimum number of operations required to convert string s1 to s2. The following operations are typically used:
% Replacing one character of string by another character.
% Deleting a character from string
% Adding a character to string
% Example:
% s1='article'
% s2='ardipo'
% EditDistance(s1,s2)
% > 4
% you need to do 4 actions to convert s1 to s2
% replace(t,d) , replace(c,p) , replace(l,o) , delete(e)
% using the other output, you can see the matrix solution to this problem
%
%
% by : Reza Ahmadzadeh (seyedreza_ahmadzadeh@yahoo.com - reza.ahmadzadeh@iit.it)
% 14-11-2012

m=length(string1);
n=length(string2);
v=zeros(m+1,n+1);
for i=1:1:m
    v(i+1,1)=i;
end
for j=1:1:n
    v(1,j+1)=j;
end
for i=1:m
    for j=1:n
        if (string1(i) == string2(j))
            v(i+1,j+1)=v(i,j);
        else
            v(i+1,j+1)=1+min(min(v(i+1,j),v(i,j+1)),v(i,j));
        end
    end
end
V=v(m+1,n+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function G = createKNNGraph(dist,K)
    G = zeros(length(dist));
    for i = 1:length(dist)
    ind = find(dist(i,:)>0); 
       if length(ind)>3
        [~,ind1] = sort(dist(i,ind)); 
        ind = ind(ind1(1:min(K,length(ind1))));  
        G(i,ind) = 1; 
       end
       G(i,i) = 1; 
    end
 end

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function G1 = updateGraph(G1,dist1,outliers,K)
    for j=1:length(outliers)
        i = outliers(j);  
        ind = find(dist1(i,:)>0);
        [~,ind1] = sort(dist(i,ind)); 
        ind = ind(ind1(1:min(K,length(ind1))));   
        G1(i,ind) = 1;
    end
end
 