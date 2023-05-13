% fluid_flow.m
%%%%%%%%%%%%%%%%%%%
%
% solving 4 IVPs
%
a=[1 2;2 -1];
[t1,x1]=ode45(@(t,x) a*x,[0,10], [-7;10]);
[t2,x2]=ode45(@(t,x) a*x,[0,10], [-4;10]);
[t3,x3]=ode45(@(t,x) a*x,[0,10], [4;-10]);
[t4,x4]=ode45(@(t,x) a*x,[0,10], [7;-10]);
%
% vector fields
%
[x,y]=meshgrid(-10:1:10, -10:1:10);
X=x+2*y; Y=2*x-y;
quiver(x,y,X,Y);
hold on 
plot(x1(:,1),x1(:,2),x2(:,1),x2(:,2),x3(:,1),x3(:,2),x4(:,1),x4(:,2),'LineWidth',2)
axis equal, axis([-10 10 -10 10])
xlabel x, ylabel y
title('a planar vector field (x+2y, 2x-y) and integral curves')
hold off
%
% 
%
