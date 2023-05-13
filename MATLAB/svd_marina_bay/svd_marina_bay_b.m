% svd_marina_bay.m
%%%%%%%%%%%%%%%%%%%%
% load image file
RGB=imread('marina_bay.jpg');
G=im2gray(RGB);
A=imresize(G, 1/5);
% uint8 -> float64
B=double(A);
% A10
[U10,S10,V10]=svds(B,10);
B10=U10*S10*V10';
A10=uint8(B10);
% A50
[U50,S50,V50]=svds(B,50);
B50=U50*S50*V50';
A50=uint8(B50);
% A90 
[U90,S90,V90]=svds(B,90);
B90=U90*S90*V90';
A90=uint8(B90);
% plot A
imshow(A);
title('A, 448 by 794, rank(A)=448');
% plot A(10)
imshow(A10);
title('A(10)');
% plot A(50)
imshow(A50);
title('A(50)');
% plot A(90)
imshow(A90);
title('A(90)');


