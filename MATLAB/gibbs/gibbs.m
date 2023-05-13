x = linspace(-1/8,5/8,500);
y=(1+sign(x)).*(1+sign(1/2-x))/4;
z=ones(1,length(x))/2;

for k=0:50
    z=z+2*sin(2*pi*(2*k+1)*x)/(pi*(2*k+1));
    plot(x,y,x,z,'LineWidth',2)
    ylim([-0.3 1.3])
    title('Gibbs Phenomenon for a step function')
    xticks([0 1/2])
    xticklabels({'0','1/2'})
    yticks([0 1])
    yticklabels({'0','1'})
    drawnow nocallbacks
    pause(0.5)
end

x = linspace(-1/4,5/4,500);
y=x.*(1+sign(x)).*(1+sign(1-x))/4+(x-1).*(1+sign(x-1)).*(1+sign(2-x))/4+(x+1).*(1+sign(x+1)).*(1+sign(-x))/4;
z=ones(1,length(x))/2;

for k=1:50
    z=z-sin(2*pi*k*x)/(pi*k);
    plot(x,y,x,z,'LineWidth',2)
    ylim([-0.3 1.3])
    xlim([-1/4 5/4])
    title('Gibbs Phenomenon for a periodic ramp function')
    xticks([0 1])
    xticklabels({'0','1'})
    yticks([0 1])
    yticklabels({'0','1'})
    drawnow nocallbacks
    pause(0.5)
end