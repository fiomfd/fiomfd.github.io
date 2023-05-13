% computer_sound_d.m
%%%%%%%%%%%%%%%%%%%%%%
% signal 
fs=1000;
Ts=1/fs;
N=12000;
t=(0:N-1)*Ts;
x=transpose(sin(2*t.*(t-3).*(t-6).*(t-9).*(t-12)));

% notation
L=256;
M=256;
N0=250;
n0=1+floor((N-L)/(L-N0));
s=L-N0;
l=transpose(1:L);
k=transpose(1:M/2+1);
m=1:M;
j=1:n0;

% window
w=hann(L);

% sampling and locarization
Y=zeros(L,n0);
Y(l,j)=w(l,1).*x(s.*(j-1)+l);

% sqrt(M) times Fourier matrix
omega=exp(2*pi*i/M);
F=zeros(M);
F(l,m)=power(omega,-(l-1).*(m-1));

% discrete short-time Fourier transform
Z=F*Y;

% spectrogram
X=zeros(M/2+1,n0);
X(k,j)=abs(Y(k,j));

% plotting
imagesc(X);
colorbar;
colormap(cool(50))
title('Spectrogram of x[n]=x(n/fs)','x(t) = sin(2t(t-3)(t-6)(t-9)(t-12))')
xlabel('n = 0,1,2..., N-1')
xticks([0 1958/2 1958])
xticklabels({'0','6000','12000'})
ylabel('k = 0, 1, ..., M/2+1')
yticks([1 129/2 129])
yticklabels({'1000','500','0'})
