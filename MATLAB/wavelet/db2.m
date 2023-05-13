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

% the wavelet 
u0=zeros([N 1]);
u1=zeros([N 1]);
u0(1)=(3+sqrt(3))/4/sqrt(2);
u0(2)=(3-sqrt(3))/4/sqrt(2);
u0(3)=(1-sqrt(3))/4/sqrt(2);
u0(4)=(1+sqrt(3))/4/sqrt(2);
u1(1)=(-3+sqrt(3))/4/sqrt(2);
u1(2)=(3+sqrt(3))/4/sqrt(2);
u1(3)=(-1-sqrt(3))/4/sqrt(2);
u1(4)=(1-sqrt(3))/4/sqrt(2);

% Fourier transform of the wavelet and absolute values
Fu0=F*u0;
Fu1=F*u1;
AbsFu0(m)=abs(Fu0(m));
AbsFu1(m)=abs(Fu1(m));

% Plotting
subplot(1,4,1)
stem(u0)
title('Scaling filter of Daubechies 2')
xlabel('n=0,1,2,...,N-1')
subplot(1,4,2)
stem(u1)
title('Wavelet filter of Daubechies 2')
xlabel('n=0,1,2,...,N-1')
subplot(1,4,[3,4])
plot(m-1,AbsFu0,'-o',m-1,AbsFu1,'-o')
legend('|Fu_0|','|Fu_1|','Location','southeast')
title('Discrete Fourier transform of Daubechies 2')
xlabel('k=0,1,2,...,N-1')