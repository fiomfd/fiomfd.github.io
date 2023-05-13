using FileIO: load, save
import LibSndFile

x=load("cheese.wav");
fs=x.samplerate
s=x.data[:, 1];
t = (0 : length(s)-1) / x.samplerate;

using Plots

# analogue signal 
plot(t,s, xaxis="time t [sec]", yaxis="x(t)", title="Analogue Signal", grid=false, legend=false)
# discrete time signal
plot(s, xaxis="n=0,1,2,3....", yaxis="x[n]", title="Discrete Time Signal", grid=false, legend=false)

using DSP, FFTW

hats=fft(s);
k = (0:length(hats)-1)*50/length(hats);
G=broadcast(abs, hats);

# DFT
plot(G,
    xaxis="k = 0, 1, 2, ..., N-1", 
    yaxis="|F(x)[k]|", 
    title="Discrete Fourier Transform", 
    grid=false, 
    legend=false)
plot!(xticks = ([0 30000 60000;], [0 30000 60000]))
plot!(yticks = ([0 200 400;], [0 200 400]))

# FT
plot(k,G,xaxis="frequency k", yaxis="|F(x)[k]|", title="Sampled Fourier Transform", grid=false, legend=false)


kshift = (-length(s)/2:length(s)/2-1)*(50/length(s));
l = (-length(hats)/2:length(hats)/2-1);
Gshift = fftshift(G);
# sampled FFT shifted
plot(kshift,Gshift,xaxis="frequency k", yaxis="|F(x)[k]|", title="Sampled Fourier Transform Shifted", grid=false, legend=false)
# FFT shifted
plot(l,Gshift,xaxis="discrete frequency k", yaxis="|F(x)[k]|", title="Discrete Fourier Transform Shifted", grid=false, legend=false)