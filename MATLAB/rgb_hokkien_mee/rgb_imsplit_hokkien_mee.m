% load image file
I=imread('hokkien_mee.jpg');
RGB=imresize(I, 1/5);
[R,G,B] = imsplit(RGB);
figure;
subplot(2,2,1)
imshow(RGB);
title('Original RGB Image');
subplot(2,2,2)
imshow(R);
title('R');
subplot(2,2,3)
imshow(G);
title('G');
subplot(2,2,4)
imshow(B);
title('B');