% wavelet_gbb2b.m
%%%%%%%%%%%%%%%%%%%%
% load image file
RGB=imread('marina_bay.jpg');
G=im2gray(RGB);
A=imresize(G, 1/5);
% multiresolution decomposition
[C,S]=wavedec2(A,3,'db2');
A3=appcoef2(C,S,'db2',3);
[H3,V3,D3] = detcoef2('all',C,S,3);
% plot A(n)
subplot(2,2,1);
imagesc(A3);
colormap('gray');
title('Approximation of Level 3');
subplot(2,2,2);
imagesc(H3);
title('Horizontal Detail of Level 3');
subplot(2,2,3);
imagesc(V3);
title('Vertical Detail of Level 3');
subplot(2,2,4);
imagesc(D3);
title('Diagonal Detail of Level 3');



