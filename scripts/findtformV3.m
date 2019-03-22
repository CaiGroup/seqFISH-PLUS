function [tform, controlimage] = findtformV3(filename1,filename2)

%first use imageJ to do a 10 pixel background subtraction
%save individual channels
%Filename2 is reference image

img = loadtiff(filename1);
img2 = loadtiff(filename2);

tform = imregcorr(img,img2,'similarity');
controlimage = imwarp(img,tform,'OutputView',imref2d(size(img2)));

falsecolorOverlay1 = imfuse(img2,controlimage);
falsecolorOverlay2 = imfuse(img2,img);
figure;
imshow(falsecolorOverlay1,'InitialMagnification','fit');
figure;
imshow(falsecolorOverlay2,'InitialMagnification','fit');

LinkFigures(1:2)
