using FileIO: load, save
import LibSndFile

# audio signal
w=load("cheese.wav");
fs=Int(w.samplerate);
N=length(w);
t = (0 :N-1)/fs;
x=w.data[:, 1];
T=N/fs;

# setting integers
L=256;
M1=128;
M=2*M1;
N0=250;
n0=1+floor(Int, (N-L)/(L-N0));
s=L-N0;

using DSP, FFTW

# window
h=hanning(L);

using LinearAlgebra

# the N-th root of unity
omega=exp(2*pi*im/M);

# sqrt(N) times the Fourier matrix
F=zeros(Complex, M, M);
for m=1:M
    for n=1:M
        F[m,n]=omega^(-(m-1)*(n-1))
    end
end

# sampling and locarization
Y=zeros(Complex, L,n0);
for l=1:L
    for j=1:n0
        Y[l,j]=h[l]*x[s*(j-1)+l]
    end
end

# discrete short-time Fourier transform
Z=F*Y;

# spectrogram
X=zeros(Float64, M1+1,n0);
for k=1:M1+1
    for j=1:n0
        X[k,j]=abs.(Y[M1+2-k,j])
    end
end

using Plots

# plotting
plot(X, 
    title="Spectrogram of human voice", 
    st=:heatmap, 
    color=:coolwarm, 
    xaxis="n = 1, 1, 2, ..., n0-1", 
    yaxis="k = 0, 1, 2, ..., M/2")
plot!(xticks = ([0 5000 10000;], [0 5000 10000]))
plot!(yticks = ([0 50 100;], [0 50 100]))

plot(X, 
    title="Spectrogram of human voice", 
    st=:heatmap, 
    color=:coolwarm, 
    xaxis="time t [sec]", 
    yaxis="frequency ξ [Hz]")
plot!(xticks = ([0 n0/T 2*n0/T 3*n0/T;], [0, 1, 2, 3]))
plot!(yticks = ([0 129*5000/11025 129*10000/11025;], [0, 5000, 10000]))

