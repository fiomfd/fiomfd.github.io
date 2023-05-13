using LinearAlgebra

# size N
N=16;

# the N-th root of unity
omega=exp(2*pi*im/N);

# sqrt(N) times the Fourier matrix
F=zeros(Complex, N, N);
for m=1:N
    for n=1:N
        F[m,n]=omega^((m-1)*(n-1))
    end
end

# wavelet filter
u0=zeros(Complex, N, 1);
u0[1,1]=(3+sqrt(3))/4/sqrt(2);
u0[2,1]=(3-sqrt(3))/4/sqrt(2);
u0[3,1]=(1-sqrt(3))/4/sqrt(2);
u0[4,1]=(1+sqrt(3))/4/sqrt(2);
u1=zeros(Complex, N, 1);
u1[1,1]=(-3+sqrt(3))/4/sqrt(2);
u1[2,1]=(3+sqrt(3))/4/sqrt(2);
u1[3,1]=(-1-sqrt(3))/4/sqrt(2);
u1[4,1]=(1-sqrt(3))/4/sqrt(2);

# Fourier transform of the wavelet and absolute values
Fu0=F*u0;
Fu1=F*u1;
AbsFu0=abs.(Fu0);
AbsFu1=abs.(Fu1);

using Plots

# plotting the scaling filter
plot(u0, 
    grid=false, 
    title=("Scaling filter of the Daubechies 2 wavelet"), 
    xaxis="n = 0, 1, 2, ..., N-1",
    st=:steppost,
    linewidth=3, 
    legend=false, 
    c=[6], 
    maker="o")
plot!(xticks = ([1:1:16;], [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]))   
plot!(yticks = ([-0.5 0 0.5], [-0.5 0 0.5]))

# plotting the wavelet filter
plot(u1, 
    grid=false, 
    title=("Wavelet filter of the Daubechies 2 wavelet"), 
    xaxis="n = 0, 1, 2, ..., N-1",    
    st=:steppost,
    linewidth=3, 
    legend=false, 
    c=[6], 
    maker="o")
plot!(xticks = ([1:1:16;], [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]))   
plot!(yticks = ([-0.5 0 0.5], [-0.5 0 0.5]))

# plot the Fourier transform
plot([AbsFu0,AbsFu1],
    grid=false,
    linewidth=3, 
    marker="o",
    title="Discrete Fourier transform of the Daubechies 2 filters",
    xaxis="k= 0, 1, 2, ..., N-1",
    label=["|Fu0|" "|Fu1|"],
    legend=:bottomright,
    maker="o")
plot!(xticks = ([1:1:16;], [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]))
plot!(yticks = ([0 1 sqrt(2)], [0 1 1.414]))
