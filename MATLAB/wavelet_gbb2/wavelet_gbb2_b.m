% wavelet_gbb2b.m
%%%%%%%%%%%%%%%%%%%%
% load image file
RGB=imread('fl2.jpg');
G=im2gray(RGB);
A=imresize(G, 1/5);
% multiresolution decomposition
[C,S]=wavedec2(A,1,'haar');
A1=appcoef2(C,S,'haar',1);
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



