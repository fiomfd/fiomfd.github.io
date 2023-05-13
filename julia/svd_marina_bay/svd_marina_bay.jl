# svd_marina_bay.jp 
using GR, LinearAlgebra, LowRankApprox, Formatting, ImageFiltering, Images, Plots, ImageMagick, Colors, TestImages, ImageView, ImageTransformations

I=load("marina_bay.jpg")

G=Gray.(I)
# size(J)

A=imresize(G, ratio=1/5)
# rank(A), size(A), eltype(A)
(p,q)=size(A)
B=Array{Float64}(A)
U, S, V=psvd(B)
B50=sum(S[n]*U[1:p,n]*(V[1:q,n])' for n=1:50);
A50=Gray.(B50)