% load image file
I=imread('shaved_ice.jpg');
RGB=imresize(I, 1/5);
[R,G,B] = imsplit(RGB);
% imshow(RGB)
% whos R, whos G, whos B

% resolution of R
[RC,RS]=wavedec2(R,5,'haar');
MRL1=appcoef2(RC,RS,'haar',1);
MRL2=appcoef2(RC,RS,'haar',2);
MRL3=appcoef2(RC,RS,'haar',3);
MRL4=appcoef2(RC,RS,'haar',4);
MRL5=appcoef2(RC,RS,'haar',5);
RL1=mat2gray(MRL1);
RL2=mat2gray(MRL2);
RL3=mat2gray(MRL3);
RL4=mat2gray(MRL4);
RL5=mat2gray(MRL5);
% resolution of G
[GC,GS]=wavedec2(G,5,'haar');
MGL1=appcoef2(GC,GS,'haar',1);
MGL2=appcoef2(GC,GS,'haar',2);
MGL3=appcoef2(GC,GS,'haar',3);
MGL4=appcoef2(GC,GS,'haar',4);
MGL5=appcoef2(GC,GS,'haar',5);
GL1=mat2gray(MGL1);
GL2=mat2gray(MGL2);
GL3=mat2gray(MGL3);
GL4=mat2gray(MGL4);
GL5=mat2gray(MGL5);
% resolution of B
[BC,BS]=wavedec2(B,5,'haar');
MBL1=appcoef2(BC,BS,'haar',1);
MBL2=appcoef2(BC,BS,'haar',2);
MBL3=appcoef2(BC,BS,'haar',3);
MBL4=appcoef2(BC,BS,'haar',4);
MBL5=appcoef2(BC,BS,'haar',5);
BL1=mat2gray(MBL1);
BL2=mat2gray(MBL2);
BL3=mat2gray(MBL3);
BL4=mat2gray(MBL4);
BL5=mat2gray(MBL5);

% R+G+B
% RGBL1
RGBL1(:,:,1)=(RL1);
RGBL1(:,:,2)=(GL1);
RGBL1(:,:,3)=(BL1);
% RGBL2
RGBL2(:,:,1)=(RL2);
RGBL2(:,:,2)=(GL2);
RGBL2(:,:,3)=(BL2);
% RGBL3
RGBL3(:,:,1)=(RL3);
RGBL3(:,:,2)=(GL3);
RGBL3(:,:,3)=(BL3);
% RGBL4
RGBL4(:,:,1)=(RL4);
RGBL4(:,:,2)=(GL4);
RGBL4(:,:,3)=(BL4);
% RGBL5
RGBL5(:,:,1)=(RL5);
RGBL5(:,:,2)=(GL5);
RGBL5(:,:,3)=(BL5);

figure;
subplot(2,3,1);
imshow(RGB);
title('Original RGB Image');
subplot(2,3,2);
imshow(RGBL1);
title('Level 1 Approximation');
subplot(2,3,3);
imshow(RGBL2);
title('Level 2 Approximation');
subplot(2,3,4);
imshow(RGBL3);
title('Level 3 Approximation');
subplot(2,3,5);
imshow(RGBL4);
title('Level 4 Approximation');
subplot(2,3,6);
imshow(RGBL5);
title('Level 5 Approximation');
