# cheese_multiresolution.jl
##########################
using FileIO: load, save
import LibSndFile

# audio signal
w=load("cheese.wav");
fs=Int(w.samplerate);
N=length(w);
t = (0 :N-1)/fs;
x=w.data[:, 1];
T=N/fs;
X=[t x];

using DSP, Wavelets

y1=dwt(X[:,2], wavelet(WT.db2), 1);
y2=dwt(X[:,2], wavelet(WT.db2), 2);
y3=dwt(X[:,2], wavelet(WT.db2), 3);
y4=dwt(X[:,2], wavelet(WT.db2), 4);

y1a=zeros(N);
y2a=zeros(N);
y3a=zeros(N);
y4a=zeros(N);
y1a[1:Int(N/2)]=y1[1:Int(N/2)];
y2a[1:Int(N/4)]=y2[1:Int(N/4)];
y3a[1:Int(N/8)]=y3[1:Int(N/8)];
y4a[1:Int(N/16)]=y4[1:Int(N/16)];
y1d=y1-y1a;
y2d=y2-y2a;
y3d=y3-y3a;
y4d=y4-y4a;

x1a=idwt(y1a, wavelet(WT.haar), 1);
x1d=idwt(y1d, wavelet(WT.haar), 1);
x2a=idwt(y2a, wavelet(WT.haar), 2);
x2d=idwt(y2d, wavelet(WT.haar), 2);
x3a=idwt(y3a, wavelet(WT.haar), 3);
x3d=idwt(y3d, wavelet(WT.haar), 3);
x4a=idwt(y4a, wavelet(WT.haar), 4);
x4d=idwt(y4d, wavelet(WT.haar), 4);

using Plots

# x1a
p1a=plot(x1a, 
    xaxis="time t [sec]", 
    yaxis="amplitude", 
    title="Level 1 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 fs 2*fs 3*fs;], [0 1 2 3]))
plot!(yticks = ([-0.8 0 0.8;], [-0.8 "0" 0.8]))
# x1d
p1d=plot(x1d, 
    xaxis="time t [sec]", 
    yaxis="amplitude", 
    title="Level 1 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 fs 2*fs 3*fs;], [0 1 2 3]))
plot!(yticks = ([-0.1 0 0.1;], [-0.1 "0" 0.1]))
# x1
plot(p1a, p1d, layout = (1,2))

# x2a
p# x2a
p2a=plot(x2a, 
    xaxis="time t [sec]", 
    yaxis="amplitude", 
    title="Level 2 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 fs 2*fs 3*fs;], [0 1 2 3]))
plot!(yticks = ([-0.8 0 0.8;], [-0.8 "0" 0.8]))
# x2d
p2d=plot(x2d, 
    xaxis="time t [sec]", 
    yaxis="amplitude", 
    title="Level 2 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 fs 2*fs 3*fs;], [0 1 2 3]))
plot!(yticks = ([-0.4 0 0.5;], [-0.4 "0" 0.5]))
# x2
plot(p2a, p2d, layout = (1,2))

# x3a
p3a=plot(x3a, 
    xaxis="time t [sec]", 
    yaxis="amplitude", 
    title="Level 3 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 fs 2*fs 3*fs;], [0 1 2 3]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x3d
p3d=plot(x3d, 
    xaxis="time t [sec]", 
    yaxis="amplitude", 
    title="Level 3 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 fs 2*fs 3*fs;], [0 1 2 3]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# plot!(yticks = ([-0.8 0 1;], [-0.8 "0" "1"]))
# x3
plot(p3a, p3d, layout = (1,2))

# x4a
p4a=plot(x4a, 
    xaxis="time t [sec]", 
    yaxis="amplitude", 
    title="Level 4 Approximation", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 fs 2*fs 3*fs;], [0 1 2 3]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x4d
p4d=plot(x4d, 
    xaxis="time t [sec]", 
    yaxis="amplitude", 
    title="Level 4 Detail", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 fs 2*fs 3*fs;], [0 1 2 3]))
plot!(yticks = ([-1 0 1;], [-1 0 1]))
# x4
plot(p4a, p4d, layout = (1,2))
