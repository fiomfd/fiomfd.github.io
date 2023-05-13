% sir_model.m
%%%%%%%%%%%%%%%%%%%%
% Solving IVP with
% total population S+I+R=1000000
% S(0)=999900
% I(0)=100
% R(0)=0
%%%%%%%%%%%%%%%%%%%%
% The case of  
% infection rate beta=1/30000000
% recovery rate gamma=1/100
%%%%%%%%%%%%%%%%%%%%
[t1,u1]=ode45(@(t,u) [-u(1)*u(2)/180000000; u(1)*u(2)/200000000-u(2)/100; u(2)/100], [0,1000], [900000;100000;0]);
%%%%%%%%%%%%%%%%%%%%
% The case of  
% infection rate beta=1/20000000
% recovery rate gamma=1/100
%%%%%%%%%%%%%%%%%%%%
[t2,u2]=ode45(@(t,u) [-u(1)*u(2)/45000000; 2*u(1)*u(2)/100000000-u(2)/100; u(2)/100], [0,1000], [900000;100000;0]);
%%%%%%%%%%%%%%%%%%%%
% ploting
%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1);
plot(t1,u1(:,1),t1,u1(:,2),t1,u1(:,3),'LineWidth',2); 
title('SIR model with S(0)=9*10^5, I(0)=10^5, R(0)=0, \beta=18^{-1}*10^{-7}, \gamma=1/100, i.e., R_0=1/2'); 
xlabel('Time t'); 
ylabel('population');
legend('S','I','R','Location','east');

subplot(2,1,2);
plot(t2,u2(:,1),t2,u2(:,2),t2,u2(:,3),'LineWidth',2);
title('SIR model with S(0)=9*10^5, I(0)=10^5, R(0)=0, \beta=45^{-1}*10^{-6}, \gamma=1/100, i.e., R_0=2');
xlabel('Time t');
ylabel('population');
legend('S','I','R','Location','east');
%%%%%%%%%%%%%%%%%%%%%%
