using CSV, Primes, DataFrames

kappa=1.4;
Q=[1 2 3];

for a=1:2000
    for b=a+1:100000
        c=a+b;
        X=factor(Vector,a*b*c);
        L=length(X);
        rad=X[1];
        if [gcd(a,b) gcd(a,c) gcd(b,c)]==[1 1 1]
            for i=2:L
                if X[i]==X[i-1]
                    rad=rad*1;
                else
                    rad=rad*X[i];
                end
            end
            if (c >= rad^kappa)
                q=[a b c];
                Q=[Q;q];
            end
        end
    end
end

R=DataFrame(Q,:auto);
CSV.write("abc14.csv", R);