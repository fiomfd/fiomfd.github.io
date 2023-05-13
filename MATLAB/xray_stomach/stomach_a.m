% define f(x,y)
f=imread('stomach.png');
% xray transform of f
theta=0:180;
[Xf,t]=radon(f,theta);
% backprojection of Xf
F=iradon(Xf,theta,'linear','Ram-Lak');
U=iradon(Xf,theta,'linear','none');
% plotting
subplot(2,2,1);
imshow(f);
title('Original Grayscale Image');
p2=subplot(2,2,2);
imagesc(theta,t,Xf);
title('R_1f(\theta,t) the X-ray transform of f(x,y)');
xlabel('\theta [degrees]');
ylabel('t');
set(gca,'XTick',0:20:180);
colormap(p2,hot);
colorbar;
subplot(2,2,3);
imshow(U,[]);
title('Unfiltered Backprojection');
subplot(2,2,4);
imshow(F,[]);
title('Filtered Backprojection');
