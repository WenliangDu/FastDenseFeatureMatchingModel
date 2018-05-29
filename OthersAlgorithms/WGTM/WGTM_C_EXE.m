function [inliers,outliers]=WGTM_C_EXE(P1,P2,K,mu)
% Weighted Graph Transformation Matching Algorithm implemented in C++: 
%          The algorithm takes two sets of corresponding points (initial matches) and find matches
%          using the geometrical relation between points.
% The WGTM algorithm is implemented in C++ and its execution file is used in 
% the following code for speeding up and utilizing in Matlab functions. 

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
fid=fopen('input.txt','w');
fprintf(fid,'%f %f %f %f %f ',double(length(P1)),double(K),double(mu),double(P1'),double(P2'));
fclose(fid);
dos('WGTM_C_Function.exe');    
fid=fopen('output.txt','r');
y=fscanf(fid,'%f');
outliers=find(y');
inliers=find(y'==0);
fclose(fid);
