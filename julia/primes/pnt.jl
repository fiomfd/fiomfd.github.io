using Primes, Plots, LinearAlgebra

P=zeros(100000);
P[1]=1;
P[2]=1;
L=zeros(100000);
L[1]=3/log(3);
L[2]=3/log(3);
N=zeros(100000);
N[1]=1;
N[2]=2;
Q=ones(100000);

for i=3:100000
    A=zeros(i);
    for j=1:i
        A[j]=Float64(isprime(j));
    end
    P[i]=sum(A);
    L[i]=i/log(i);
    N[i]=i;
end

p1=plot([P L],  
    grid=false,
    linewidth=2, 
    title="Prime Number Theorem: Pi(x) and x/log(x)", 
    right_margin=Plots.Measures.Length(:mm, 10.0),
    left_margin=Plots.Measures.Length(:mm, 5.0),
    xlabel="x",
    xticks = ([1 50000 100000;], ["e" "50000" "100000"]),
    yaxis="Number of Primes",
    legendfont=font(10), 
    label=["Pi(x)" "x/log(x)"],
    palette = :seaborn_bright, 
    legend = :topleft)
savefig("pnt1.png") 

p2=plot([P./L Q],  
    grid=false,
    linewidth=2, 
    title="Prime Number Theorem: Pi(x)/(x/log(x))", 
    right_margin=Plots.Measures.Length(:mm, 10.0),
    left_margin=Plots.Measures.Length(:mm, 5.0),
    xlabel="x",
    xticks = ([1 50000 100000;], ["e" "50000" "100000"]),
    yaxis="",
    yticks = ([0.4 0.8 1.0 1.2;], ["0.4" "0.8" "1.0" "1.2"]),
    label = :false,
    palette = :seaborn_bright, 
    legend = :topright)
savefig("pnt2.png") 

plot(p1, p2,
     layout=(1,2), 
     size=(1200,400), 
     left_margin=Plots.Measures.Length(:mm, 5.0),
     right_margin=Plots.Measures.Length(:mm, 15.0),
     top_margin=Plots.Measures.Length(:mm, 5.0),
     bottom_margin=Plots.Measures.Length(:mm, 5.0))
savefig("pnt.png") 
