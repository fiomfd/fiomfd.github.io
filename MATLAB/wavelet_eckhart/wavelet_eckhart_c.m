% wavelet_gbb2c.m
%%%%%%%%%%%%%%%%%%%%
% load image file
RGB=imread('eckhart.jpg');
A=im2gray(RGB);
% multiresolution decomposition
[C,S]=wavedec2(A,2,'coif1');
A2=appcoef2(C,S,'coif1',2);
[H2,V2,D2] = detcoef2('all',C,S,2);
% plot A(n)
subplot(2,2,1);
imagesc(A2);
colormap('gray');
title('Approximation of Level 2');
subplot(2,2,2);
imagesc(H2);
title('Horizontal Detail of Level 2');
subplot(2,2,3);
imagesc(V2);
title('Vertical Detail of Level 2');
subplot(2,2,4);
imagesc(D2);
title('Diagonal Detail of Level 2');



