using LinearAlgebra, LowRankApprox, Formatting, ImageFiltering, Images, Plots, ImageMagick, Colors, TestImages, ImageView, ImageTransformations

A=rand(80:255, 20, 20);
B=rand(80:255, 20, 20);
C=rand(80:255, 20, 20);
X=zeros(3,20,20);
X[1,:,:]=A/255;
X[2,:,:]=B/255;
X[3,:,:]=C/255;
P=colorview(RGB,X);

p1=plot(Gray.(A/255), 
    xaxis=false, 
    xticks=false, 
    yaxis=false, 
    yticks=false, 
    grid=false, 
    title="Grayscale Image (RED)")

p2=plot(Gray.(B/255), 
    xaxis=false, 
    xticks=false, 
    yaxis=false, 
    yticks=false, 
    grid=false, 
    title="Grayscale Image (GREEN)")

p3=plot(Gray.(C/255), 
    xaxis=false, 
    xticks=false, 
    yaxis=false, 
    yticks=false, 
    grid=false, 
    title="Grayscale Image (BLUE)")

p4=plot(P, 
    xaxis=false, 
    xticks=false, 
    yaxis=false, 
    yticks=false, 
    grid=false, 
    title="RGB Image")

plot(p1, p2, p3, p4, 
     layout=(2,2), 
     size=(1000,1000), 
     margin=Plots.Measures.Length(:mm, 2.0))
savefig("rgb_matrix.png") 