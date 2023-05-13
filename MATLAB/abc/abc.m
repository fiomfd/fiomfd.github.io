% Examples of abc theorem

kappa=1.4;
Q=[1 2 3];

for a=1:5
    for b=a+1:180000
        c=a+b;
        X=factor(a*b*c);
        L=length(X);
        rad=X(1,1);
        if [gcd(a,b) gcd(a,c) gcd(b,c)]==[1 1 1]
            for i=2:L
                if X(1,i)==X(1,i-1)
                    rad=rad*1;
                else
                    rad=rad*X(1,i);
                end
            end
            if (c >= rad^kappa)
                q=[a b c]
                Q=[Q;q];
            end
        end
    end
end

writematrix(Q,'abc14.csv')