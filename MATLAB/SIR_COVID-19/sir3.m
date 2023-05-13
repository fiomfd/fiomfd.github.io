% sir3.m
%%%%%%%%%%%%%%%%%%%%%%
% Simulation of SIR model for COVID-19 in the next D days 
% for A1, A2, A3, ....
%
% An: area, n=1,2,3,... 
% t: time [days]
% d0: staring date
% d1: ending date
% d2: starting date of taking the average t=-d
% d2 < d0 < d1
% D=d1-d0 [days]
% N: total population
%
% S=u(1): the number of susceptible individuals
% I=u(2): active cases
% R=u(3): the number of dischrged individuals
% S+I+R=N always holds.
% I+R: total cases
% 
% T0: total cases at t=0
% T2: daily cases at t=-d
% R0: total dischrged at t=0
% R2: total discharged at t=-d
% beta: infection rate
% gamma: recovery rate
%
% Input d0, D. 
% For n=1,2,3,..., input And2, 
% And2, An, AnT0, AnT1, AnR0, AnR1 and AnN.
% AnS(0), AnI(0), AnR(0), Anbeta, Angamma,
% are determined by AnT0, AnT1, AnR0, AnR1, AnN d0, d1 and And2. 
% Run the program in MATLAB. 
% MATLAB will show you the result.
%%%%%%%%%%%%%%%%%%%%%%%
% d0 starting date
d0=datetime(2021,1,2);
D=240
% d1 ending date
d1=d0+days(D)

% initial data
%%%%%%%%%%%%%%%%%%%%%%%
% A1
A1="Tokyo";
A1T0=60960;
A1R0=50861;
A1d2=datetime(2020,12,27);
A1T2=56559;
A1R2=48442;
A1N=13999568;
% A2
A2="South Korea";
A2T0=63244;
A2R0=44507;
A2d2=datetime(2020,12,20);
A2T2=45475;
A2R2=35155;
A2N=51779000;
% A3
A3="Okinawa";
A3T0=5424;
A3R0=4967;
A3d2=datetime(2020,12,29);
A3T2=5260;
A3R2=4847;
A3N=8815082;

% Computing
%%%%%%%%%%%%%%%%%%%%%
% period D [days]
D=days(d1-d0);
l0=datestr(d0,'dd mmm yyyy');
l1=datestr(d0+days(100),'dd mmm yyyy');
l2=datestr(d0+days(200),'dd mmm yyyy');
% A1
A1T1=(A1T0-A1T2)/days(d0-A1d2);
A1R1=(A1R0-A1R2)/days(d0-A1d2);
A1beta=(A1T1)/(A1T0-A1R0)/(A1N-A1T0);
A1gamma=A1R1/(A1T0-A1R0);
[A1t,A1u]=ode45(@(t,u) [-A1beta*u(1)*u(2); A1beta*u(1)*u(2)-A1gamma*u(2); A1gamma*u(2)], [0,D], [A1N-A1T0;A1T0-A1R0;A1R0]);
% A2
A2T1=(A2T0-A2T2)/days(d0-A2d2);
A2R1=(A2R0-A2R2)/days(d0-A2d2);
A2beta=(A2T1)/(A2T0-A2R0)/(A2N-A2T0);
A2gamma=A2R1/(A2T0-A2R0);
[A2t,A2u]=ode45(@(t,u) [-A2beta*u(1)*u(2); A2beta*u(1)*u(2)-A2gamma*u(2); A2gamma*u(2)], [0,D], [A2N-A2T0;A2T0-A2R0;A2R0]);
% A3
A3T1=(A3T0-A3T2)/days(d0-A3d2);
A3R1=(A3R0-A3R2)/days(d0-A3d2);
A3beta=(A3T1)/(A3T0-A3R0)/(A3N-A3T0);
A3gamma=A3R1/(A3T0-A3R0);
[A3t,A3u]=ode45(@(t,u) [-A3beta*u(1)*u(2); A3beta*u(1)*u(2)-A3gamma*u(2); A3gamma*u(2)], [0,D], [A3N-A3T0;A3T0-A3R0;A3R0]);

% Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%
% A1
subplot(1,3,1);
plot(A1t,A1u(:,2)+A1u(:,3),A1t,A1u(:,2),A1t,A1u(:,3),'LineWidth',2);
title({'SIR model for COVID-19',[A1]});
xticks([0 100 200]);
xticklabels({[l0],[l1],[l2]});
ylabel('population');
legend('Total Cases','Active Cases','Discharged','Location','northwest');
% A2
subplot(1,3,2);
plot(A2t,A2u(:,2)+A2u(:,3),A2t,A2u(:,2),A2t,A2u(:,3),'LineWidth',2);
title({'SIR model for COVID-19',[A2]});
xticks([0 100 200]);
xticklabels({[l0],[l1],[l2]});
ylabel('population');
legend('Total Cases','Active Cases','Discharged','Location','northwest');
% A3
subplot(1,3,3);
plot(A3t,A3u(:,2)+A3u(:,3),A3t,A3u(:,2),A3t,A3u(:,3),'LineWidth',2);
title({'SIR model for COVID-19',[A3]});
xticks([0 100 200]);
xticklabels({[l0],[l1],[l2]});
ylabel('population');
legend('Total Cases','Active Cases','Discharged','Location','northwest');
