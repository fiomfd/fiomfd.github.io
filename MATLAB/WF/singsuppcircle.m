% semiclassical parameter
h=1/1000;

% f(x,y)

load circleGS;
A=circleGS/255;

% FBI transform of f
F0=zeros(120,120);
F1=zeros(120,120);
F2=zeros(120,120);
F3=zeros(120,120);
F4=zeros(120,120);
F5=zeros(120,120);
F6=zeros(120,120);
F7=zeros(120,120);
F8=zeros(120,120);

for j=1:120
    for k=1:120
        for q=-8:8
            for p=-8:8
                E3=exp(-(p/50)^2/(2*h));
                E4=exp(-(q/50)^2/(2*h));
                F0(j,k)=F0(j,k)+p*q*A(2*j+7+p,2*k+7+q)*E3*E4;
                F1(j,k)=F1(j,k)+p*A(2*j+7+p,2*k+7+q)*E3*E4;
                F2(j,k)=F2(j,k)+q*A(2*j+7+p,2*k+7+q)*E3*E4;
                F3(j,k)=F3(j,k)+(p+q)*A(2*j+7+p,2*k+7+q)*E3*E4;
                F4(j,k)=F4(j,k)+(2*p+q)*A(2*j+7+p,2*k+7+q)*E3*E4;
                F5(j,k)=F5(j,k)+(p+2*q)*A(2*j+7+p,2*k+7+q)*E3*E4;
                F6(j,k)=F6(j,k)+(p-q)*A(2*j+7+p,2*k+7+q)*E3*E4;
                F7(j,k)=F7(j,k)+(2*p-q)*A(2*j+7+p,2*k+7+q)*E3*E4;
                F8(j,k)=F8(j,k)+(p-2*q)*A(2*j+7+p,2*k+7+q)*E3*E4;
          end
        end
    end
end

F0=abs(F0);
M0=max(max(F0));
F1=abs(F1);
M1=max(max(F1));
F2=abs(F2);
M2=max(max(F2));
F3=abs(F3);
M3=max(max(F3));
F4=abs(F4);
M4=max(max(F4));
F5=abs(F5);
M5=max(max(F5));
F6=abs(F6);
M6=max(max(F6));
F7=abs(F7);
M7=max(max(F7));
F8=abs(F8);
M8=max(max(F8));

X1=zeros(1,1);
Y1=zeros(1,1);
Z1=zeros(2,1);

for j=1:120
    for k=1:120
        if (F0(j,k)>80*M0/100)
           X1=[X1 (k-60)/50];
           Y1=[Y1 (60-j)/50];
           Z1=[Z1 [j;k]];
        end
    end
end

for j=1:120
    for k=1:120
        if (F1(j,k)>70*M1/100)
           X1=[X1 (k-60)/50];
           Y1=[Y1 (60-j)/50];
           Z1=[Z1 [j;k]];
        end
    end
end

for j=1:120
    for k=1:120
        if (F2(j,k)>70*M2/100)
           X1=[X1 (k-60)/50];
           Y1=[Y1 (60-j)/50];
           Z1=[Z1 [j;k]];
        end
    end
end

for j=1:120
    for k=1:120
        if (F3(j,k)>80*M3/100)
           X1=[X1 (k-60)/50];
           Y1=[Y1 (60-j)/50];
           Z1=[Z1 [j;k]];
        end
    end
end

for j=1:120
    for k=1:120
        if (F4(j,k)>80*M4/100)
           X1=[X1 (k-60)/50];
           Y1=[Y1 (60-j)/50];
           Z1=[Z1 [j;k]];
        end
    end
end

for j=1:120
    for k=1:120
        if (F5(j,k)>80*M5/100)
           X1=[X1 (k-60)/50];
           Y1=[Y1 (60-j)/50];
           Z1=[Z1 [j;k]];
        end
    end
end

for j=1:120
    for k=1:120
        if (F6(j,k)>80*M6/100)
           X1=[X1 (k-60)/50];
           Y1=[Y1 (60-j)/50];
           Z1=[Z1 [j;k]];
        end
    end
end

for j=1:120
    for k=1:120
        if (F7(j,k)>80*M7/100)
           X1=[X1 (k-60)/50];
           Y1=[Y1 (60-j)/50];
           Z1=[Z1 [j;k]];
        end
    end
end

for j=1:120
    for k=1:120
        if (F8(j,k)>80*M8/100)
           X1=[X1 (k-60)/50];
           Y1=[Y1 (60-j)/50];
           Z1=[Z1 [j;k]];
        end
    end
end

L=length(X1);
X=X1(1,2:L);
Y=Y1(1,2:L);
Z=int8([Z1(1,2:L);Z1(2,2:L)]);

dx1=subplot(1,2,1);
imagesc(A);
title('Original Grayscale Image');
xlabel('x');
ylabel('y');
xticks([28 128 228]);
xticklabels({'-1','0','1'});
yticks([28 128 228]);
yticklabels({'1','0','-1'});
pbaspect([1 1 1]);
colormap(dx1,bone);
caxis([0 1]);
%
subplot(1,2,2)
scatter(X(1,:),Y(1,:),'MarkerFaceColor',[1 0 0])
title('Singular support of a circle');
xlabel('x');
xticks([-1 0 1]);
xticklabels({'-1','0','1'});
ylabel('y');
yticks([-1 0 1]);
yticklabels({'-1','0','1'});
pbaspect([1 1 1]);
set(gcf,'Position',[300,200,900,400]);
saveas(gcf,'scatter-circle.png');
