using Plots

Y1=zeros(3,32);
Y2=zeros(3,32);
Y3=zeros(3,32);
Y4=zeros(3,32);

for j=1:32
    Y1[:,j]=[1;1;2*pi*(j-1)/32];
    Y2[:,j]=[-1;1;2*pi*(j-1)/32];   
    Y3[:,j]=[1;-1;2*pi*(j-1)/32];
    Y4[:,j]=[-1;-1;2*pi*(j-1)/32];
end

Y51=zeros(3,41);
Y52=zeros(3,41);
Y61=zeros(3,41);
Y62=zeros(3,41);
Y71=zeros(3,41);
Y72=zeros(3,41);
Y81=zeros(3,41);
Y82=zeros(3,41);

for j=1:41
    Y51[:,j]=[(j-21)/20;1;pi/2];
    Y52[:,j]=[(j-21)/20;-1;3*pi/2];
    Y61[:,j]=[(j-21)/20;1;3*pi/2];
    Y62[:,j]=[(j-21)/20;-1;pi/2];    
    Y71[:,j]=[1;(j-21)/20;0];
    Y72[:,j]=[-1;(j-21)/20;pi];
    Y81[:,j]=[1;(j-21)/20;pi];
    Y82[:,j]=[-1;(j-21)/20;0];
end

l01=scatter(Y51[1,:], Y51[2,:], Y51[3,:], color="purple")
l02=scatter(l01, Y61[1,:], Y61[2,:], Y61[3,:], color="magenta")
l03=scatter(l02, Y72[1,:], Y72[2,:], Y72[3,:], color="cyan")
l04=scatter(l03, Y82[1,:], Y82[2,:], Y82[3,:], color="gold")
l05=scatter(l04, Y2[1,:], Y2[2,:], Y2[3,:], color="lime")
l06=scatter(l05, Y52[1,:], Y52[2,:], Y52[3,:], color="purple")
l07=scatter(l06, Y62[1,:], Y62[2,:], Y62[3,:], color="magenta")
l08=scatter(l07, Y71[1,:], Y71[2,:], Y71[3,:], color="cyan")
l09=scatter(l08, Y81[1,:], Y81[2,:], Y81[3,:], color="gold")
l10=scatter(l09, Y1[1,:], Y1[2,:], Y1[3,:], color="red")
l11=scatter(l10, Y4[1,:], Y4[2,:], Y4[3,:], color="orange")
p5=scatter(l11, Y3[1,:], Y3[2,:], Y3[3,:], color="blue",
    title="The wave front set of a square", 
    size=(600,400),
    xaxis="x",
    xticks=([-1 0 1;],[-1 0 1]),
    yaxis="y", 
    yticks=([-1 0 1;],[-1 0 1]),
    zaxis="θ",
    zticks=([0 pi/2 pi 3*pi/2;],["0" "π/2" "π" "3π/2"]),
    legend = :false,
    camera=(40,35))
savefig("./xray_square_5.png") 

X1=zeros(3,64);
X2=zeros(3,64);
X3=zeros(3,64);
X4=zeros(3,64);
for j=1:64
    X1[:,j]=[pi*(j-1)/64;cos(pi*(j-1)/64)+sin(pi*(j-1)/64);sin(pi*(j-1)/64)-cos(pi*(j-1)/64)];
    X2[:,j]=[pi*(j-1)/64;-cos(pi*(j-1)/64)+sin(pi*(j-1)/64);-sin(pi*(j-1)/64)-cos(pi*(j-1)/64)];
    X3[:,j]=[pi*(j-1)/64;cos(pi*(j-1)/64)-sin(pi*(j-1)/64);sin(pi*(j-1)/64)+cos(pi*(j-1)/64)];
    X4[:,j]=[pi*(j-1)/64;-cos(pi*(j-1)/64)-sin(pi*(j-1)/64);-sin(pi*(j-1)/64)+cos(pi*(j-1)/64)];
end

X5=zeros(3,21);
X6=zeros(3,21);
X7=zeros(3,21);
X8=zeros(3,21);
for j=1:21
    X5[:,j]=[pi/2;1;(j-11)/10];
    X6[:,j]=[pi/2;-1;(j-11)/10];
    X7[:,j]=[0;1;(j-11)/10];
    X8[:,j]=[0;-1;(j-11)/10];
end

L1=scatter(X7[1,:], X7[2,:], X7[3,:], color="cyan")
L2=scatter(L1, X5[1,:], X5[2,:], X5[3,:], color="purple")
L3=scatter(L2, X8[1,:], X8[2,:], X8[3,:], color="gold") 
L4=scatter(L3, X6[1,:], X6[2,:], X6[3,:], color="magenta")
L5=scatter(L4, X1[1,:], X1[2,:], X1[3,:], color="red")
L6=scatter(L5, X2[1,:], X2[2,:], X2[3,:], color="lime")
L7=scatter(L6, X3[1,:], X3[2,:], X3[3,:], color="blue")
p6=scatter(L7, X4[1,:], X4[2,:], X4[3,:], color="orange", 
    title="The wave front set of the sinogram of a square", 
    size=(600,400),
    xaxis="θ",
    xticks=([0 pi/2 pi;],["0" "π/2" "π"]),
    yaxis="t", 
    yticks=([-sqrt(2) -1 0 1 sqrt(2);],[-1.41 -1 0 1 1.41]),
    zaxis="μ",
    zticks=([-sqrt(2) -1 0 1 sqrt(2);],[-1.41 -1 0 1 1.41]),
    right_margin=Plots.Measures.Length(:mm, 10.0),
    legend = :false,
    camera=(10,40))
savefig("./xray_square_6.png")