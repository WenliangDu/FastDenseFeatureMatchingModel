
function [disparityMapNew,matchedLocation2New] = ShowDisparity_sub_alpha2(matchedLocation2,disparityMap)
% [hm01a,hm01b] = size(hm01);
% 
% Newhm01 = hm01;
% WindowSize = 5;
% for i = WindowSize:hm01a-WindowSize
%     for j = WindowSize:hm01b-WindowSize
%         
%         CurrentWindow = hm01(i-WindowSize+1:i+WindowSize-1,j-WindowSize+1:j+WindowSize-1);
%         
%         
%     end
% end

matchedLocation2L = size(matchedLocation2,1);
Wrongrecord = false(1,matchedLocation2L);

WindowSize = 5;
matchedLocation2KDT = KDTreeSearcher(matchedLocation2);
IdxNN = knnsearch(matchedLocation2KDT,matchedLocation2,'K',WindowSize+1);

for i = 1:matchedLocation2L
    CurrentdisparityWindow = disparityMap(IdxNN(i,2:WindowSize+1));
    if disparityMap(i) > 5+median(CurrentdisparityWindow)
        Wrongrecord(i) = true;
    elseif disparityMap(i)+5 < median(CurrentdisparityWindow)
        Wrongrecord(i) = true;
    end
end
% sum(Wrongrecord)
disparityMapNew = disparityMap;
disparityMapNew(Wrongrecord) = [];

matchedLocation2New = matchedLocation2;
matchedLocation2New(Wrongrecord,:) = [];

% [xm,ym] = meshgrid(1:size(I2,2),1:size(I2,1));
% hm = interpolate_alpha1(double(matchedLocation2New(:,1)), double(matchedLocation2New(:,2)), double(disparityMapNew), xm, ym);
% 
% figure,mesh(xm, ym, hm);
% colormap jet