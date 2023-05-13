% svd_eckhart.m
%%%%%%%%%%%%%%%%%%%%
% load image file
RGB=imread('eckhart.jpg');
A=im2gray(RGB);
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
% plot A
imshow(A);
title('A, 480 by 640, rank(A)=480');
% plot A(20)
imshow(A20);
title('A(20)');
% plot A(50)
imshow(A50);
title('A(50)');
% plot A(80)
imshow(A80);
title('A(80)');


