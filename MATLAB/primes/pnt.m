P=zeros(1,10000);
P(1,1)=1;
P(1,2)=1;
L=zeros(1,10000);
L(1,1)=3/log(3);
L(1,2)=3/log(3);
N=zeros(1,10000);
N(1,1)=1;
N(1,2)=2;
Q=zeros(1,10000);
Q=Q+1;

for i=3:10000
    A=zeros(1,i);
    for j=1:i
        A(1,j)=double(isprime(j));
    end
    P(1,i)=sum(A(1,:));
    L(1,i)=i/log(i);
    N(1,i)=i;
end

newcolors = [0 0 1; 
             0 1 0; 
             1 0 0];
colororder(newcolors)

subplot(1,2,1)
plot(N,P,N,L,'LineWidth',2)
title('The Prime Number Theorem: Pi(x) and x/log(x)')
xlabel('x');
xticks([1 5000 10000]);
xticklabels({'e','5000','10000'});
legend('Pi(x)','x/log(x)','Location','northwest');

subplot(1,2,2)
plot(N,P./L,N,Q,'LineWidth',2)
title('The Prime Number Theorem: Pi(x)/(x/log(x))')
xlabel('x');
xticks([1 5000 10000]);
xticklabels({'e','5000','10000'});

set(gcf,'Position',[600,200,1200,400]);
saveas(gcf,'pnt.png');