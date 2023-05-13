% FBI transform of f(x,y)

% semiclassical parameter
h=1/1000;

% f(x,y)
A=zeros(501);
for j=51:451
    A(j,51:j)=1;
    A(j,j)=1/2;
end

% FBI transform of f
F0=zeros(47,47);
F1=zeros(47,47);
F2=zeros(47,47);

for j=1:47
    for k=1:47
        for q=-10:10
            for p=-10:10
                E3=exp(-(-p/20)^2/(2*h));
                E4=exp(-(q/20)^2/(2*h));
                F0(j,k)=F0(j,k)+p*q*A(10*j+1+p,10*k+1+q)*E3*E4;
                F1(j,k)=F1(j,k)+p*A(10*j+1+p,10*k+1+q)*E3*E4;
                F2(j,k)=F2(j,k)+q*A(10*j+1+p,10*k+1+q)*E3*E4;
            end
        end
    end
end

% The absolute value of FBI transform of f

F0=abs(F0);
M0=max(max(F0));
F1=abs(F1);
M1=max(max(F1));
F2=abs(F2);
M2=max(max(F2));

Z10=zeros(2,1);
Z11=zeros(2,1);
Z12=zeros(2,1);

for j=1:47
    for k=1:47
        if (F0(j,k)>50*M0/100)
            Z10=[Z10 [j;k]];
        end
    end
end

for j=1:47
    for k=1:47
        if (F1(j,k)>60*M1/100)
            Z11=[Z11 [j;k]];
        end
    end
end

for j=1:47
    for k=1:47
        if (F2(j,k)>60*M2/100)
            Z12=[Z12 [j;k]];
        end
    end
end

L0=length(Z10);
Z0=Z10(:,2:L0);
L1=length(Z11);
Z1=Z11(:,2:L1);
L2=length(Z12);
Z2=Z12(:,2:L2);

% FBI transform of f
G0=zeros(L0-1,32);
G1=zeros(L1-1,32);
G2=zeros(L2-1,32);

% semiclassical parameter
h0=1/100;
lambda0=10*sqrt(2)*h0;
for l=1:32
    for n=1:L0-1
        for q=-10:10
            for p=-10:10
                E1=exp(1i*lambda0*(p/20)*sin(pi*(l-1)/16)/h0);
                E2=exp(-1i*lambda0*(q/20)*cos(pi*(l-1)/16)/h0);
                E3=exp(-(-p/20)^2/(2*h0));
                E4=exp(-(q/20)^2/(2*h0));
                E5=p/20-1i*lambda0*sin(pi*(l-1)/16);
                E6=q/20+1i*lambda0*cos(pi*(l-1)/16);
                G0(n,l)=G0(n,l)+A(10*Z0(1,n)+1+p,10*Z0(2,n)+1+q)*E1*E2*E3*E4*E5*E6;
            end
        end
    end
end

% semiclassical parameter
h1=1/10;
lambda1=16*sqrt(2)*h1;
for l=1:32
    for n=1:L1-1
        for q=-10:10
            for p=-10:10
                E1=exp(1i*lambda1*(p/20)*sin(pi*(l-1)/16)/h1);
                E2=exp(-1i*lambda1*(q/20)*cos(pi*(l-1)/16)/h1);
                E3=exp(-(-p/20)^2/(2*h1));
                E4=exp(-(q/20)^2/(2*h1));
                E5=p/20-1i*lambda1*sin(pi*(l-1)/16);
                G1(n,l)=G1(n,l)+A(10*Z1(1,n)+1+p,10*Z1(2,n)+1+q)*E1*E2*E3*E4*E5;
            end
        end
    end
end

% semiclassical parameter
h2=1/10;
lambda2=16*sqrt(2)*h2;
for l=1:32
    for n=1:L2-1
        for q=-10:10
            for p=-10:10
                E1=exp(1i*lambda2*(p/20)*sin(pi*(l-1)/16)/h2);
                E2=exp(-1i*lambda2*(q/20)*cos(pi*(l-1)/16)/h2);
                E3=exp(-(p/20)^2/(2*h2));
                E4=exp(-(q/20)^2/(2*h2));
                E6=q/20+1i*lambda2*cos(pi*(l-1)/16);
                G2(n,l)=G2(n,l)+A(10*Z2(1,n)+1+p,10*Z2(2,n)+1+q)*E1*E2*E3*E4*E6;
            end
        end
    end
end


% The absolute value of FBI transform of f

G0=abs(G0);
M0=max(G0);
m0=max(max(M0));
G1=abs(G1);
M1=max(G1);
m1=max(max(M1));
G2=abs(G2);
M2=max(G2);
m2=max(max(M2));

U1=zeros(1,1);
V1=zeros(1,1);
W1=zeros(1,1);

for l=1:32
    for n=1:L0-1
        if (G0(n,l)>40*m0/1000)
                U1=[U1 (Z0(2,n)-24)/20];
                V1=[V1 (24-Z0(1,n))/20];
                W1=[W1 pi*(l-1)/16];
        end
    end
end

for l=1:32
    for n=1:L1-1
        if (G1(n,l)>60*m1/100)
                U1=[U1 (Z1(2,n)-24)/20];
                V1=[V1 (24-Z1(1,n))/20];
                W1=[W1 pi*(l-1)/16];
        end
    end
end

for l=1:32
    for n=1:L2-1
        if (G2(n,l)>60*m2/100)
                U1=[U1 (Z2(2,n)-24)/20];
                V1=[V1 (24-Z2(1,n))/20];
                W1=[W1 pi*(l-1)/16];
        end
    end
end

L3=length(U1);
U=U1(1,2:L3);
V=V1(1,2:L3);
W=W1(1,2:L3);

dx1=subplot(1,2,1);
imagesc(A);
title('Original Grayscale Image');
xlabel('x');
ylabel('y');
xticks([51 251 451]);
xticklabels({'-1','0','1'});
yticks([51 251 451]);
yticklabels({'1','0','-1'});
pbaspect([1 1 1]);
colormap(dx1,bone);
caxis([0 1]);
%
subplot(1,2,2)
scatter3(U(1,:),V(1,:),W(1,:),'MarkerFaceColor',[1 0 0])
title('The wave front set of the characteristic function of a triangle','detected by the FBI transform');
xlabel('x');
xticks([-1 0 1]);
xticklabels({'-1','0','1'});
ylabel('y');
yticks([-1 0 1]);
yticklabels({'-1','0','1'});
zlabel('\theta');
zticks([0 pi/2 pi 3*pi/2]);
zticklabels({'0','\pi/2','\pi','3\pi/2'});
pbaspect([1 1 1]);
view(-30,10)
set(gcf,'Position',[1200,100,900,400]);
saveas(gcf,'wf_triangle.png');
