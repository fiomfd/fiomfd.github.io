% size N
N=16;

% the numbers of rows and columns 
m=transpose(1:N);
n=1:N;

% the N-th root of unity
omega=exp(2*pi*i/N);

% sqrt(N) times the Fourier matrix
F=zeros(N);
F(m,n)=power(omega,(m-1).*(n-1));

% the haar wavelet 
u0=zeros([N 1]);
u1=zeros([N 1]);
u0(1)=1/sqrt(2);
u0(2)=1/sqrt(2);
u1(1)=1/sqrt(2);
u1(2)=-1/sqrt(2);

% the Daubechies 2wavelet 
v0=zeros([N 1]);
v1=zeros([N 1]);
v0(1)=(3+sqrt(3))/4/sqrt(2);
v0(2)=(3-sqrt(3))/4/sqrt(2);
v0(3)=(1-sqrt(3))/4/sqrt(2);
v0(4)=(1+sqrt(3))/4/sqrt(2);
v1(1)=(-3+sqrt(3))/4/sqrt(2);
v1(2)=(3+sqrt(3))/4/sqrt(2);
v1(3)=(-1-sqrt(3))/4/sqrt(2);
v1(4)=(1-sqrt(3))/4/sqrt(2);

% FFT of haar
Fu0=F*u0;
Fu1=F*u1;
AbsFu0(m)=abs(Fu0(m));
AbsFu1(m)=abs(Fu1(m));

% FFT of db2
Fv0=F*v0;
Fv1=F*v1;
AbsFv0(m)=abs(Fv0(m));
AbsFv1(m)=abs(Fv1(m));

% Plotting
subplot(2,4,1)
stem(u0)
title('Scaling filter of Haar')
xlabel('n=0,1,2,...,N-1')
xlim auto
subplot(2,4,2)
stem(u1)
title('Wavelet filter of Haar')
xlabel('n=0,1,2,...,N-1')
subplot(2,4,[3,4])
plot(m-1,AbsFu0,'-o',m-1,AbsFu1,'-o')
legend('|Fu_0|','|Fu_1|','Location','southeast')
title('Discrete Fourier transform of Haar')
xlabel('k=0,1,2,...,N-1')
subplot(2,4,5)
stem(v0)
title('Scaling filter of Daubechies 2')
xlabel('n=0,1,2,...,N-1')
subplot(2,4,6)
stem(v1)
title('Wavelet filter of Daubechies 2')
xlabel('n=0,1,2,...,N-1')
subplot(2,4,[7,8])
plot(m-1,AbsFv0,'-o',m-1,AbsFv1,'-o')
legend('|Fu_0|','|Fu_1|','Location','southeast')
title('Discrete Fourier transform of Daubechies 2')
xlabel('k=0,1,2,...,N-1')
