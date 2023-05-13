% define f(x,y)
f=zeros(400);
f(1:400, 201:400)=1;
% xray transform of f
theta=0:180;
[Xf,t]=radon(f,theta);
% backprojection of Xf
F=iradon(Xf,theta);
U=iradon(Xf,theta,'linear','none');
% plotting
subplot(2,2,1);
imshow(f);
title('Original Grayscale Image');
subplot(2,2,2);
imagesc(theta,t,Xf);
title('Xf(\theta,t) the X-ray transform of f(x,y)');
xlabel('\theta [degrees]');
ylabel('t');
set(gca,'XTick',0:20:180);
colormap(hot);
colorbar;
subplot(2,2,3);
imshow(U,[]);
title('Unfiltered Backprojection');
subplot(2,2,4);
imshow(F,[]);
title('Filtered Backprojection');
