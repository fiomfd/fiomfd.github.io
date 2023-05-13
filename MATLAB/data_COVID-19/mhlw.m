websave("./csv/mhlw_new_cases.csv","https://covid19.mhlw.go.jp/public/opendata/newly_confirmed_cases_daily.csv");
websave("./csv/mhlw_cases.csv","https://covid19.mhlw.go.jp/public/opendata/confirmed_cases_cumulative_daily.csv");
websave("./csv/mhlw_deaths.csv","https://covid19.mhlw.go.jp/public/opendata/deaths_cumulative_daily.csv");
websave("./csv/mhlw_recovered.csv","https://covid19.mhlw.go.jp/public/opendata/requiring_inpatient_care_etc_daily.csv");
XN=importdata("./csv/mhlw_new_cases.csv");
XC=importdata("./csv/mhlw_cases.csv");
XD=importdata("./csv/mhlw_deaths.csv");
XR=importdata("./csv/mhlw_recovered.csv");
l=length(XC.data(:,1));
m=length(XN.data(:,1));
n=floor(l);
k=floor(m);
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

dd0=datetime('2021-07-01');
DD=1+days(d1-dd0);
ll0=datestr(dd0,'yyyy-mm-dd');
ll1=datestr(dd0+days(floor((DD-1)/2)),'yyyy-mm-dd');
ll2=datestr(d1,'yyyy-mm-dd');

% Japan
PJPN=125.845010;
XJPN=VN(:,1)/PJPN;
CJPN=VC(:,1)/PJPN;
DJPN=VD(:,1)/PJPN;
RJPN=VR(:,2)/PJPN;
AJPN=CJPN-RJPN-DJPN;
NJPN=zeros(D,1);
for j=1:7
    NJPN(j,1)=(XJPN(j+104,1)+XJPN(j+103,1)+XJPN(j+102,1)+XJPN(j+101,1)+XJPN(j+100,1)+XJPN(j+99,1)+XJPN(j+98,1))/7;
end
for j=8:D
    NJPN(j,1)=(CJPN(j,1)-CJPN(j-7,1))/7;
end
NDJPN=zeros(D,1);
for j=1:7
    NDJPN(j,1)=(DJPN(j,1)-DJPN(1,1))/j;
end
for j=8:D
    NDJPN(j,1)=(DJPN(j,1)-DJPN(j-7,1))/7;
end

% Tokyo
PTKY=13.988129;
XTKY=VN(:,14)/PTKY;
CTKY=VC(:,14)/PTKY;
DTKY=VD(:,14)/PTKY;
RTKY=VR(:,41)/PTKY;
ATKY=CTKY-RTKY-DTKY;
NTKY=zeros(D,1);
for j=1:7
    NTKY(j,1)=(XTKY(j+104,1)+XTKY(j+103,1)+XTKY(j+102,1)+XTKY(j+101,1)+XTKY(j+100,1)+XTKY(j+99,1)+XTKY(j+98,1))/7;
end
for j=8:D
    NTKY(j,1)=(CTKY(j,1)-CTKY(j-7,1))/7;
end
NDTKY=zeros(D,1);
for j=1:7
    NDTKY(j,1)=(DTKY(j,1)-DTKY(1,1))/j;
end
for j=8:D
    NDTKY(j,1)=(DTKY(j,1)-DTKY(j-7,1))/7;
end

% Okinawa
POKNW=1.469335;    
XOKNW=VN(:,48)/POKNW;
COKNW=VC(:,48)/POKNW;
DOKNW=VD(:,48)/POKNW;
ROKNW=VR(:,143)/POKNW;
AOKNW=COKNW-ROKNW-DOKNW;
NOKNW=zeros(D,1);
for j=1:7
    NOKNW(j,1)=(XOKNW(j+104,1)+XOKNW(j+103,1)+XOKNW(j+102,1)+XOKNW(j+101,1)+XOKNW(j+100,1)+XOKNW(j+99,1)+XOKNW(j+98,1))/7;
end
for j=8:D
    NOKNW(j,1)=(COKNW(j,1)-COKNW(j-7,1))/7;
end
NDOKNW=zeros(D,1);
for j=1:7
    NDOKNW(j,1)=(DOKNW(j,1)-DOKNW(1,1))/j;
end
for j=8:D
    NDOKNW(j,1)=(DOKNW(j,1)-DOKNW(j-7,1))/7;
end

% Osaka
POSK=8.797153;
XOSK=VN(:,28)/POSK;
COSK=VC(:,28)/POSK;
DOSK=VD(:,28)/POSK;
ROSK=VR(:,83)/POSK;
AOSK=COSK-ROSK-DOSK;
NOSK=zeros(D,1);
for j=1:7
    NOSK(j,1)=(XOSK(j+104,1)+XOSK(j+103,1)+XOSK(j+102,1)+XOSK(j+101,1)+XOSK(j+100,1)+XOSK(j+99,1)+XOSK(j+98,1))/7;
end
for j=8:D
    NOSK(j,1)=(COSK(j,1)-COSK(j-7,1))/7;
end
NDOSK=zeros(D,1);
for j=1:7
    NDOSK(j,1)=(DOSK(j,1)-DOSK(1,1))/j;
end
for j=8:D
    NDOSK(j,1)=(DOSK(j,1)-DOSK(j-7,1))/7;
end

% Hyogo
PHYG=5.425850;
XHYG=VN(:,29)/PHYG;
CHYG=VC(:,29)/PHYG;
DHYG=VD(:,29)/PHYG;
RHYG=VR(:,86)/PHYG;
AHYG=CHYG-RHYG-DHYG;
NHYG=zeros(D,1);
for j=1:7
    NHYG(j,1)=(XHYG(j+104,1)+XHYG(j+103,1)+XHYG(j+102,1)+XHYG(j+101,1)+XHYG(j+100,1)+XHYG(j+99,1)+XHYG(j+98,1))/7;
end
for j=8:D
    NHYG(j,1)=(CHYG(j,1)-CHYG(j-7,1))/7;
end
NDHYG=zeros(D,1);
for j=1:7
    NDHYG(j,1)=(DHYG(j,1)-DHYG(1,1))/j;
end
for j=8:D
    NDHYG(j,1)=(DHYG(j,1)-DHYG(j-7,1))/7;
end

% Hokkaido
PHKD=5.191355;
XHKD=VN(:,2)/PHKD;
CHKD=VC(:,2)/PHKD;
DHKD=VD(:,2)/PHKD;
RHKD=VR(:,5)/PHKD;
AHKD=CHKD-RHKD-DHKD;
NHKD=zeros(D,1);
for j=1:7
    NHKD(j,1)=(XHKD(j+104,1)+XHKD(j+103,1)+XHKD(j+102,1)+XHKD(j+101,1)+XHKD(j+100,1)+XHKD(j+99,1)+XHKD(j+98,1))/7;
end
for j=8:D
    NHKD(j,1)=(CHKD(j,1)-CHKD(j-7,1))/7;
end
NDHKD=zeros(D,1);
for j=1:7
    NDHKD(j,1)=(DHKD(j,1)-DHKD(1,1))/j;
end
for j=8:D
    NDHKD(j,1)=(DHKD(j,1)-DHKD(j-7,1))/7;
end

newcolors = [0 0 1; 
             1 131/255 0; 
             0 1 0; 
             1 0 0; 
             138/255 43/255 226/255; 
             169/255 80/255 69/255;
             1 0 1; 
             220/255 220/255 220/255; 
             1 191/255 17/255; 
             0 191/255 1];
colororder(newcolors) 
subplot(2,2,1)
plot([CJPN CTKY COKNW COSK CHKD],'LineWidth',2)
title('COVID-19: cases per 1M','data sourced by MoH of Japan')
xlabel('date')
ylabel('cases/1M')
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[lf]})
legend('Japan','Tokyo','Okinawa','Osaka','Hokkaido','Location','northwest');
%
subplot(2,2,2)
plot([DJPN DTKY DOKNW DOSK DHKD],'LineWidth',2)
title('COVID-19: deaths per 1M','data sourced by MoH of Japan')
xlabel('date')
ylabel('deaths/1M')
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[lf]})
legend('Japan','Tokyo','Okinawa','Osaka','Hokkaido','Location','northwest');
%
subplot(2,2,3)
plot([NJPN NTKY NOKNW NOSK NHKD],'LineWidth',2)
title('COVID-19: 7-day average of new cases per 1M','data sourced by MoH of Japan')
xlabel('date')
ylabel('cases/1M')
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[lf]})
legend('Japan','Tokyo','Okinawa','Osaka','Hokkaido','Location','northwest');
%
subplot(2,2,4)
plot([NDJPN NDTKY NDOKNW NDOSK NDHKD],'LineWidth',2)
title('COVID-19: 7-day average of daily deaths per 1M','data sourced by MoH of Japan')
xlabel('date')
ylabel('deaths/1M')
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[lf]})
legend('Japan','Tokyo','Okinawa','Osaka','Hokkaido','Location','northwest');
%
set(gcf,'Position',[600,200,1200,800]);
saveas(gcf,'mhlw.png');
