using Images, TestImages, ImageView, LinearAlgebra, Plots, FFTW, DSP
include("xray.jl")

z=zeros(401,802);

for p=101:301
    for q=201-Int64(floor(sqrt(100^2-(p-201)^2))):201+Int64(floor(sqrt(100^2-(p-201)^2)))
        z[p,q]=1;
    z[p,q+401]=1;
    end
end

w,h = size(z);
Nr_div2 = floor(sqrt(w*w+h*h)/2);
Nr = Int64(Nr_div2*2+1);

angles = (0:1:179)/180*pi; # [rad]
sinogram = radon(z,angles);

theta = (0:1:359)/360*2*pi; # [rad]
xray = radon(z,theta);

p1=plot(Gray.(z),
    title="Original Grayscale Image", 
#    right_margin=Plots.Measures.Length(:mm, 10.0),
    size=(600,400), 
    xticks = ([0 h/2  h;], [-h/2 0 h/2]),
    xlabel="x",
    yticks = ([0 w/2  w;], [w/2 0 -w/2]),
    ylabel="y")
savefig("./xray_implants_1.png") 

p2=plot(sinogram, 
    size=(600,400), 
    title="Sinogram: X-ray transform Xf(θ,t)", 
    st=:heatmap, 
    color=:jet1, 
    xaxis="θ [rad]", 
    yaxis="t")
plot!(xticks = ([1 90 180;], [0 "π/2" "π"]))
plot!(yticks = ([0 Nr/2 Nr;], [-Nr/2 0 Nr/2]))
savefig("xray_implants_2.png")

#
(a,b)=size(sinogram);
p=Int64(floor(w-2));
q=Int64(floor(h-2));

U=zeros(p,q,b);
V=zeros(p,q);
t=zeros(p,q,b);
s=zeros(p,q,b);
    
for i=1:p
    for j=1:q
        for k=1:b
            t[i,j,k]=floor((j-q/2-1/2)*cos(pi*(k-1)/b)+(p/2+1/2-i)*sin(pi*(k-1)/b)+a/2+1/2);
            s[i,j,k]=ceil((j-q/2-1/2)*cos(pi*(k-1)/b)+(p/2+1/2-i)*sin(pi*(k-1)/b)+a/2+1/2);
            U[i,j,k]=sinogram[Int64(t[i,j,k]),k]/2+sinogram[Int64(s[i,j,k]),k]/2;
       end
        V[i,j]=sum(U[i,j,:])/(2*b);
    end
end

(M,)=findmax(V);
(N,)=findmin(V);
W=V-N*ones(p,q);

p3=plot(Gray.(W/(M-N)),
    title="Unfiltered Back-Projection", 
#    right_margin=Plots.Measures.Length(:mm, 10.0),
    size=(600,400), 
    xticks = ([0 q/2  q;], [-q/2 0 q/2]),
    xlabel="x",
    yticks = ([0 p/2  p;], [p/2 0 -p/2]),
    ylabel="y")
savefig("./xray_implants_3.png") 

# 
filter=zeros(Complex, a);
hf=ones(a)-hanning(a);
for k=1:a
    filter[k]=hf[k]*min(k-1,a+1-k);
end

omega=exp(2*pi*im/a)
K=zeros(Complex, a, a);
for m=1:a
    for n=1:a
        for k=1:a
        K[m,n]=K[m,n]+2*pi*filter[k]*omega^((m-n)*(k-1))/(a^2);
        end
    end
end

SINOGRAM=zeros(Complex, a, b);
for k=1:b
    SINOGRAM[:,k]=K*sinogram[:,k];
end

#

U=zeros(Complex, p,q,b);
V=zeros(Complex, p,q);
t=zeros(p,q,b);
s=zeros(p,q,b);
    
for i=1:p
    for j=1:q
        for k=1:b
            t[i,j,k]=floor((j-q/2-1/2)*cos(pi*(k-1)/b)+(p/2+1/2-i)*sin(pi*(k-1)/b)+a/2+1/2);
            s[i,j,k]=ceil((j-q/2-1/2)*cos(pi*(k-1)/b)+(p/2+1/2-i)*sin(pi*(k-1)/b)+a/2+1/2);
            U[i,j,k]=SINOGRAM[Int64(t[i,j,k]),k]/2+SINOGRAM[Int64(s[i,j,k]),k]/2;
       end
        V[i,j]=sum(U[i,j,:])/(2*b);
    end
end

W=real.(V);

p4=plot(Gray.(W),    
    title="Filtered Back-Projection", 
#    right_margin=Plots.Measures.Length(:mm, 10.0),
    size=(600,400), 
    xticks = ([0 q/2  q;], [-q/2 0 q/2]),
    xlabel="x",
    yticks = ([0 p/2  p;], [p/2 0 -p/2]),
    ylabel="y")
savefig("./xray_implants_4.png") 
#

plot(p1, p2, p3, p4,
     layout=4, 
     size=(1200,800), 
     left_margin=Plots.Measures.Length(:mm, 5.0),
     right_margin=Plots.Measures.Length(:mm, 15.0),
     top_margin=Plots.Measures.Length(:mm, 5.0),
     bottom_margin=Plots.Measures.Length(:mm, 5.0))
savefig("./xray_implants.png") 

###
art=zeros(Complex, a, b);
for k=1:b
    art[:,k]=K*(sinogram[:,k].^2);
end

U=zeros(Complex, p,q,b);
V=zeros(Complex, p,q);
t=zeros(p,q,b);
s=zeros(p,q,b);
    
for i=1:p
    for j=1:q
        for k=1:b
            t[i,j,k]=floor((j-q/2-1/2)*cos(pi*(k-1)/b)+(p/2+1/2-i)*sin(pi*(k-1)/b)+a/2+1/2);
            s[i,j,k]=ceil((j-q/2-1/2)*cos(pi*(k-1)/b)+(p/2+1/2-i)*sin(pi*(k-1)/b)+a/2+1/2);
            U[i,j,k]=art[Int64(t[i,j,k]),k]/2+art[Int64(s[i,j,k]),k]/2;
       end
        V[i,j]=sum(U[i,j,:])/(2*b);
    end
end

W=real.(V);
(M,)=findmax(W);
(N,)=findmin(W);
X=W-N*ones(p,q);

p5=plot(Gray.(X/(M-N)),    
    title="Streaking artifacts", 
#    right_margin=Plots.Measures.Length(:mm, 10.0),
    size=(600,400), 
    xticks = ([0 q/2  q;], [-q/2 0 q/2]),
    xlabel="x",
    yticks = ([0 p/2  p;], [p/2 0 -p/2]),
    ylabel="y")
savefig("./xray_implants_5.png") 

####################################
plot(p1, p2, p4, p5,
     layout=4, 
     size=(1200,800), 
     left_margin=Plots.Measures.Length(:mm, 5.0),
     right_margin=Plots.Measures.Length(:mm, 15.0),
     top_margin=Plots.Measures.Length(:mm, 5.0),
     bottom_margin=Plots.Measures.Length(:mm, 5.0))
savefig("./analysis_square.png") 