% wavelet_gbb2b.m
%%%%%%%%%%%%%%%%%%%%
% load image file
RGB=imread('eckhart.jpg');
A=im2gray(RGB);
% multiresolution decomposition
[C,S]=wavedec2(A,1,'coif1');
A1=appcoef2(C,S,'coif1',1);
[H1,V1,D1] = detcoef2('all',C,S,1);
% plot A(n)
subplot(2,2,1);
imagesc(A1);
colormap('gray');
title('Approximation of Level 1');
subplot(2,2,2);
imagesc(H1);
title('Horizontal Detail of Level 1');
subplot(2,2,3);
imagesc(V1);
title('Vertical Detail of Level 1');
subplot(2,2,4);
imagesc(D1);
title('Diagonal Detail of Level 1');



