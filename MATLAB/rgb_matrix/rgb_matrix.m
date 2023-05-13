A=mat2gray(rand(10,15));
B=mat2gray(rand(10,15));
C=mat2gray(rand(10,15));

X(:,:,1)=A;
X(:,:,2)=B;
X(:,:,3)=C;

subplot(2,2,1)
imshow(A);
title('Grayscale Image (RED)');
subplot(2,2,2)
imshow(B);
title('Grayscale Image (GREEN)');
subplot(2,2,3)
imshow(C);
title('Grayscale Image (BLUE)');
subplot(2,2,4)
imshow(X);
title('RGB Image');
