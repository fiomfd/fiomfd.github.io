# computer_sound_f.jl
##########################
fs=Int(1000);
Ts=1/fs;
N=Int(12000);
T=N/fs;
#t=(0:Int(N-1))*Ts;

x=zeros(Float64, N)
for m=1:N
    x[m]=sin(2*Ts*(m-1)*(Ts*(m-1)-3)*(Ts*(m-1)-6)*(Ts*(m-1)-9)*(Ts*(m-1)-12))
end

using DSP, Wavelets

x1=dwt(x, wavelet(WT.haar), 1);
x2=dwt(x, wavelet(WT.haar), 2);
x3=dwt(x, wavelet(WT.haar), 3);
x4=dwt(x, wavelet(WT.haar), 4);
x5=dwt(x, wavelet(WT.haar), 5);

x1a=zeros(6000);
x1d=zeros(6000);
x2a=zeros(3000);
x2d=zeros(3000);
x3a=zeros(1500);
x3d=zeros(1500);
x4a=zeros(750);
x4d=zeros(750);
x5a=zeros(375);
x5d=zeros(375);
x1a=x1[1:6000];
x1d=x1[6001:12000];
x2a=x2[1:3000];
x2d=x2[3001:6000];
x3a=x3[1:1500];
x3d=x3[1501:3000];
x4a=x4[1:750];
x4d=x4[751:1500];
x5a=x5[1:375];
x5d=x5[375:750];


using Plots

# x1
plot(x1, 
    xaxis="approximation and detail", 
    yaxis="amplitude", 
    title="Level 1 Decomposition", 
    grid=false, 
    legend=false)
plot!(xticks =([6000;],[]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x1a
plot(x1a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 1 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000;], [0 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x1d
plot(x1d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 1 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 6000;], [0 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))

# x2
plot(x2, 
    xaxis="approximation and detail", 
    yaxis="amplitude", 
    title="Level 2 Decomposition", 
    grid=false, 
    legend=false)
plot!(xticks =([6000;],[]))
plot!(yticks = ([-2 0 2;], [-2 0 2]))
# x2a
plot(x2a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 2 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 3000;], [0 12000]))
plot!(yticks = ([-2 0 2;], [-2 0 2]))
# x2d
plot(x2d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 2 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 3000;], [0 12000]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))

# x3
plot(x3, 
    xaxis="approximation and detail", 
    yaxis="amplitude", 
    title="Level 3 Decomposition", 
    grid=false, 
    legend=false)
plot!(xticks =([6000;],[]))
plot!(yticks = ([-3 0 3;], [-3 0 3]))
# x3a
plot(x3a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 3 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 1500;], [0 12000]))
plot!(yticks = ([-3 0 3;], [-3 0 3]))
# x3d
plot(x3d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 3 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 1500;], [0 12000]))
plot!(yticks = ([-2 0 2;], [-2 0 2]))

# x4
plot(x4, 
    xaxis="approximation and detail", 
    yaxis="amplitude", 
    title="Level 4 Decomposition", 
    grid=false, 
    legend=false)
plot!(xticks =([6000;],[]))
plot!(yticks = ([-3 0 3;], [-3 0 3]))
# x4a
plot(x4a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 4 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 750;], [0 12000]))
plot!(yticks = ([-3 0 3;], [-3 0 3]))
# x4d
plot(x4d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 4 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 750;], [0 12000]))
plot!(yticks = ([-3 0 3;], [-3 0 3]))

# x5
plot(x5, 
    xaxis="approximation and detail", 
    yaxis="amplitude", 
    title="Level 5 Decomposition", 
    grid=false, 
    legend=false)
plot!(xticks =([6000;],[]))
plot!(yticks = ([-5 0 5;], [-5 0 5]))
# x5a
plot(x5a, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 375;], [0 12000]))
plot!(yticks = ([-5 0 5;], [-5 0 5]))
# x5d
plot(x5d, 
    xaxis="n = 0, 1, 2, ..., N-1", 
    yaxis="amplitude", 
    title="Level 5 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 375;], [0 12000]))
plot!(yticks = ([-4 0 4;], [-4 0 4]))

