% wavelet_gbb2.m
%%%%%%%%%%%%%%%%%%%%
% load image file
RGB=imread('eckhart.jpg');
A=im2gray(RGB);
% multiresolution decomposition
[C,S]=wavedec2(A,3,'coif1');
A1=appcoef2(C,S,'coif1',1);
A2=appcoef2(C,S,'coif1',2);
A3=appcoef2(C,S,'coif1',3);
% plot A(n)
subplot(2,2,1);
imagesc(A);
colormap('gray');
title('Level 0, i.e., Original Omage');
subplot(2,2,2);
imagesc(A1);
colormap('gray');
title('Approximation of Level 1');
subplot(2,2,3);
imagesc(A2);
colormap('gray');
title('Approximation of Level 2');
subplot(2,2,4);
imagesc(A3);
colormap('gray');
title('Approximation of Level 3');



