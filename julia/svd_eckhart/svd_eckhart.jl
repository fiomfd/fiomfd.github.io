using GR, LinearAlgebra, LowRankApprox, Formatting, ImageFiltering, Images, Plots, ImageMagick, Colors, TestImages, ImageView, ImageTransformations

I=load("eckhart.jpg")

G=Gray.(I)
# size(J)

A=imresize(G, ratio=1/5)
# rank(A), size(A), eltype(A)

(p,q)=size(A)
B=Array{Float64}(A)
U, S, V=psvd(B)
B010=sum(S[n]*U[1:p,n]*(V[1:q,n])' for n=1:10);
B040=sum(S[n]*U[1:p,n]*(V[1:q,n])' for n=1:40);
B070=sum(S[n]*U[1:p,n]*(V[1:q,n])' for n=1:70);
B100=sum(S[n]*U[1:p,n]*(V[1:q,n])' for n=1:100);
A010=Gray.(B010);
A040=Gray.(B040);
A070=Gray.(B070);
A100=Gray.(B100);