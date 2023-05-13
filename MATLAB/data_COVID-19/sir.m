websave("./csv/mhlw_new_cases.csv","https://covid19.mhlw.go.jp/public/opendata/newly_confirmed_cases_daily.csv");
websave("./csv/mhlw_cases.csv","https://covid19.mhlw.go.jp/public/opendata/confirmed_cases_cumulative_daily.csv");
websave("./csv/mhlw_deaths.csv","https://covid19.mhlw.go.jp/public/opendata/deaths_cumulative_daily.csv");
websave("./csv/mhlw_recovered.csv","https://covid19.mhlw.go.jp/public/opendata/requiring_inpatient_care_etc_daily.csv");
XN=importdata("./csv/mhlw_new_cases.csv");
XC=importdata("./csv/mhlw_cases.csv");
XD=importdata("./csv/mhlw_deaths.csv");
XR=importdata("./csv/mhlw_recovered.csv");
l=length(XC.data(:,1));
VN=XN.data(:,:);
VC=XC.data(:,:);
VD=XD.data(:,:);
VR=XR.data(:,:);

d0=datetime('2020-05-09');
xd1=strrep(XC.textdata(l+1,1),'/','-');
d1=datetime(xd1);
D=1+days(d1-d0);
l0=datestr(d0,'yyyy-mm-dd');
l1=datestr(d0+days(floor((D-1)/3)),'yyyy-mm-dd');
l2=datestr(d0+days(floor(2*(D-1)/3)),'yyyy-mm-dd');
lf=datestr(d1,'yyyy-mm-dd');

d2=d1+days(28);
T=28;
ll0=datestr(d1,'yyyy-mm-dd');
ll1=datestr(d1+days(14),'yyyy-mm-dd');
ll2=datestr(d2,'yyyy-mm-dd');

% Japan
PJPN=125.845010;
XJPN=VN(:,1)/PJPN;
CJPN=VC(:,1)/PJPN;
DJPN=VD(:,1)/PJPN;
RJPN=VR(:,2)/PJPN;
AJPN=CJPN-RJPN-DJPN;

% Tokyo
PTKY=13.988129;
XTKY=VN(:,14)/PTKY;
CTKY=VC(:,14)/PTKY;
DTKY=VD(:,14)/PTKY;
RTKY=VR(:,41)/PTKY;
ATKY=CTKY-RTKY-DTKY;

% Okinawa
POKNW=1.469335;    
XOKNW=VN(:,48)/POKNW;
COKNW=VC(:,48)/POKNW;
DOKNW=VD(:,48)/POKNW;
ROKNW=VR(:,143)/POKNW;
AOKNW=max(0,COKNW-ROKNW-DOKNW);

newcolors = [0 0 1; 
             1 0 0; 
             0 1 0;];
colororder(newcolors)         
%c=hsv(6);
%colororder(c);

subplot(2,3,1)
plot([CJPN,AJPN,RJPN],'LineWidth',1)
title('COVID-19: Japan (125.36M)','data sourced by MOH of Japan')
xlabel('date')
ylabel('cases/1M')
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[lf]})
legend('confirmed','active','recovered','Location','northwest');
% plot
subplot(2,3,2)
plot([CTKY,ATKY,RTKY],'LineWidth',1)
title('COVID-19: Tokyo (14.049146M)','data sourced by MOH of Japan')
xlabel('date')
ylabel('cases/1M')
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[lf]})
legend('confirmed','active','recovered','Location','northwest');
% plot
subplot(2,3,3)
plot([COKNW,AOKNW,ROKNW],'LineWidth',1)
title('COVID-19: Okinawa (1.458870M)','data sourced by MOH of Japan')
xlabel('date')
ylabel('cases/1M')
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[lf]})
legend('confirmed','active','recovered','Location','northwest');
%
d2=d1+days(28);
W=28;
w0=datestr(d1,'yyyy-mm-dd');
w1=datestr(d1+days(14),'yyyy-mm-dd');
w2=datestr(d2,'yyyy-mm-dd');
POP=1000000;

% Japan
T0=CJPN(D);
R0=RJPN(D);
T2=CJPN(D-6);
R2=RJPN(D-6);
T1=(T0-T2)/7;
R1=(R0-R2)/7;
beta=T1/(T0-R0)/(POP-T0);
gamma=R1/(T0-R0);
[t,u]=ode45(@(t,u) [beta*(POP-u(1))*u(2); beta*(POP-u(1))*u(2)-gamma*u(2)], [0,W], [T0;T0-R0]);
%
subplot(2,3,4)
plot(t,u(:,1),t,u(:,2),t,u(:,1)-u(:,2),'LineWidth',1); 
title('SIR model for COVID-19: Japan (125.36M)','data sourced by MOH of Japan')
xticks([1 15 29]);
xticklabels({[w0],[w1],[w2]});
ylabel('cases/1M');
legend('confirmed','active','recovered','Location','east');

% Tokyo
T0=CTKY(D);
R0=RTKY(D);
T2=CTKY(D-6);
R2=RTKY(D-6);
T1=(T0-T2)/7;
R1=(R0-R2)/7;
beta=T1/(T0-R0)/(POP-T0);
gamma=R1/(T0-R0);
[t,u]=ode45(@(t,u) [beta*(POP-u(1))*u(2); beta*(POP-u(1))*u(2)-gamma*u(2)], [0,W], [T0;T0-R0]);
%
subplot(2,3,5)
plot(t,u(:,1),t,u(:,2),t,u(:,1)-u(:,2),'LineWidth',1); 
title('SIR model for COVID-19: Tokyo (14.049146M)','data sourced by MOH of Japan')
xticks([1 15 29]);
xticklabels({[w0],[w1],[w2]});
ylabel('cases/1M');
legend('confirmed','active','recovered','Location','east');

% Okinawa
T0=COKNW(D);
R0=ROKNW(D);
T2=COKNW(D-6);
R2=ROKNW(D-6);
T1=(T0-T2)/7;
R1=(R0-R2)/7;
beta=T1/(T0-R0)/(POP-T0);
gamma=R1/(T0-R0);
[t,u]=ode45(@(t,u) [beta*(POP-u(1))*u(2); beta*(POP-u(1))*u(2)-gamma*u(2)], [0,W], [T0;T0-R0]);
%
subplot(2,3,6)
plot(t,u(:,1),t,u(:,2),t,u(:,1)-u(:,2),'LineWidth',1); 
title('SIR model for COVID-19: Okinawa (1.458870M)','data sourced by MOH of Japan')
xticks([1 15 29]);
xticklabels({[w0],[w1],[w2]});
ylabel('cases/1M');
legend('confirmed','active','recovered','Location','east');
set(gcf,'Position',[600,200,1500,800]);
saveas(gcf,'sir.png');
