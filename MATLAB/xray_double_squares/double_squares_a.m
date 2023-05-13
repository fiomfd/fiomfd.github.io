% define f(x,y)
f=zeros(500);
f(101:400, 101:400)=1;
f(201:300, 201:300)=0;
% xray transform of f
theta=0:180;
[Xf,t]=radon(f,theta);
% backprojection of Xf
F=iradon(Xf,theta,'Linear','Cosine');
U=iradon(Xf,theta,'linear','none');
% plotting
p2=subplot(2,2,1);
imshow(f);
title('Original Grayscale Image');
subplot(2,2,2);
imagesc(theta,t,Xf);
title('Xf(\theta,t) the X-ray transform of f(x,y)');
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
title('Filtered Backprojection (Cosine)');
