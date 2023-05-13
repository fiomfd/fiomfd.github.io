x = linspace(-1,1.5,1000);

% (1+x)^{-1}
y=zeros(1,length(x));
for n=0:20
    y=y+power(-1,n)*power(x,n);
    plot(x,y,x,power(1+x,-1),'LineWidth',2)
    ylim([0 10])
    title('Convergence and divergence of the Taylor expansion of (1+x)^{-1}')
    xticks([-1 0 1])
    xticklabels({'-1','0','-1'})
    yticks([0 1 2 3 4 5 6 7 8 9 10])
    yticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    drawnow nocallbacks
    pause(1)
end

% (1+x)^{-1/2}
y=zeros(1,length(x));
for n=0:30
    y=y+power(-1,n)*factorial(2*n)*power(x,n)/power(2,2*n)/factorial(n)/factorial(n);
    plot(x,y,x,power(1+x,-1/2),'LineWidth',2)
    ylim([0 5])
    title('Convergence and divergence of the Taylor expansion of (1+x)^{-1/2}')
    xticks([-1 0 1])
    xticklabels({'-1','0','-1'})
    yticks([0 1 2 3 4 5])
    yticklabels({'0','1','2','3','4','5'})
    drawnow nocallbacks
    pause(1)
end

% log(1+x)
y=zeros(1,length(x));
for n=1:20
    y=y+power(-1,n-1)*power(x,n)/n;
    plot(x,y,x,log(1+x),'LineWidth',2)
    ylim([-3 1.5])
    title('Convergence and divergence of the Taylor expansion of log(1+x)')
    xticks([-1 0 1])
    xticklabels({'-1','0','-1'})
    yticks([-3 -2 -1 0 1])
    yticklabels({'-3','-2','-1','0','1'})
    drawnow nocallbacks
    pause(1)
end

% (1+x)^{1/2}
y=ones(1,length(x));
for n=1:20
    y=y+power(-1,n-1)*factorial(2*n-2)*power(x,n)/power(2,2*n-1)/factorial(n-1)/factorial(n);
    plot(x,y,x,power(1+x,1/2),'LineWidth',2)
    ylim([-0.2 2])
    title('Convergence and divergence of the Taylor expansion of (1+x)^{1/2}')
    xticks([-1 0 1])
    xticklabels({'-1','0','-1'})
    yticks([0 1 sqrt(2)])
    yticklabels({'0','1','√2'})
    drawnow nocallbacks
    pause(1)
end
