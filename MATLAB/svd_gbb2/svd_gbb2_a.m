% svd_gbb2.m
%%%%%%%%%%%%%%%%%%%%
% load image file
RGB=imread('fl2.jpg');
G=im2gray(RGB);
A=imresize(G, 1/5);
% uint8 -> float64
B=double(A);
% A20
[U20,S20,V20]=svds(B,20);
B20=U20*S20*V20';
A20=uint8(B20);
% A50
[U50,S50,V50]=svds(B,50);
B50=U50*S50*V50';
A50=uint8(B50);
% A80
[U80,S80,V80]=svds(B,80);
B80=U80*S80*V80';
A80=uint8(B80);
% plot A(n)
subplot(2,2,1);
imshow(A);
title('A, 365 by 548, rank(A)=365');
subplot(2,2,2);
imshow(A20);
title('A(20)');
subplot(2,2,3);
imshow(A50);
title('A(50)');
subplot(2,2,4);
imshow(A80);
title('A(80)');


