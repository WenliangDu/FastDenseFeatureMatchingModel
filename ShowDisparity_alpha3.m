function ShowDisparity_alpha3(I1, I2, tform1, tform2,inliersPar,indexPairsPar,Location1,Location2)

[I1Rect, I2Rect] = rectifyStereoImages(I1, I2, tform1, tform2);
%% Get disparityMap
% Location1 = validBlobs1.Location;
% Location2 = validBlobs2.Location;
Location1Rect = transformPointsForward(tform1, Location1);
Location2Rect = transformPointsForward(tform2, Location2);
matchedLocation1 = Location1Rect(indexPairsPar(inliersPar,1),:);
matchedLocation2 = Location2Rect(indexPairsPar(inliersPar,2),:);

disparityMap = matchedLocation1(:,1) - matchedLocation2(:,1); 

[disparityMapNew,matchedLocation2New] = ShowDisparity_sub_alpha2(matchedLocation2,disparityMap);

%% Interpolate
[xm,ym] = meshgrid(1:size(I2Rect,2),1:size(I2Rect,1));
hm = interpolate_alpha1(double(matchedLocation2New(:,1)), double(matchedLocation2New(:,2)), double(disparityMapNew), xm, ym);
%%
MaxDisp = max(hm(:));
MinDisp = min(hm(:));
ToltalInterval = MaxDisp - MinDisp;

hm01 = ((hm - MinDisp)./ToltalInterval);
% hm01 = hm;
% hm01IMJ = imadjust(hm01);
%%
% hm01FL = hm01;
hm01FL = fliplr(hm01);
figure,mesh(xm, ym, hm01FL);
colormap jet
title('DTM-3D');
% axis equal
% axis(xy,'equal')

%% histogram equalization hm01IMJ
% MaxDisp = max(hm(:));
% MinDisp = min(hm(:));
% ToltalInterval = MaxDisp - MinDisp;
% 
% hm01 = ((hm - MinDisp)./ToltalInterval);
% hm01IMJ = imadjust(hm01);
% 
% %% Mean filter
% hm01IMJM = medfilt1(hm01IMJ,140);
% %% flip to the real rotation
% hm01IMJMF = fliplr(hm01IMJM);
% 
% %% Uint16
% hmUint16 = uint16(hm01IMJMF.*65536);
% 
% figure,mesh(xm, ym, hm);
% colormap jet
% axis equal
%% Fusion
%
% tform1 = projective2d(t1);
% tform2 = projective2d(t2);
% [I1Rect, I2Rect] = rectifyStereoImages(I1, I2, tform1, tform2);

hmUint16 = uint16(hm01FL.*65536);
I2d = double(I2Rect);
% I2d = double(I2); % Fresh
I2Uint16 = uint16((I2d./255).*65536);
hmUint16Fu = fliplr(hmUint16)+I2Uint16.*0.1;
hmUint16FuR = rot90(hmUint16Fu,2); %% rotate 180 degree
figure,imshow(hmUint16FuR,'Colormap',jet(65536));
title('DTM-2D');