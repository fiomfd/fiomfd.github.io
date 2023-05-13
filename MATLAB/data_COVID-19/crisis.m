websave('./csv/jhu_cases.csv','https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv');
W=importdata('./csv/jhu_cases.csv');
[pw qw]=size(W.data);
websave('./csv/jhu_deaths.csv','https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv');
U=importdata('./csv/jhu_deaths.csv');
[pu qu]=size(U.data);

DD=qw-985;
dd0=datetime(2022,10,1);
dd1=dd0+days(DD-1);
ll0=datestr(dd0,'yyyy-mm-dd');
ll1=string('2022-11-01');
ll2=string('2022-12-01');
ll3=string('2023-01-01');
llf=datestr(dd1,'yyyy-mm-dd');

% NSW
PNSW=8.21;
WNSW=transpose(W.data(11,3:qw))/PNSW;
NSW=zeros(DD,1);
for j=1:DD
    NSW(j,1)=max(0,WNSW(j+983,1)-WNSW(j+976,1))/7;
end

% Victoria
PVIC=6.64;
WVIC=transpose(W.data(16,3:qw))/PVIC;
VIC=zeros(DD,1);
for j=1:DD
    VIC(j,1)=max(0,WVIC(j+983,1)-WVIC(j+976,1))/7;
end

% Australia
PAUS=25.981402;
WACT=transpose(W.data(10,3:qw));
WNSW=transpose(W.data(11,3:qw));
WNT=transpose(W.data(12,3:qw));
WQLD=transpose(W.data(13,3:qw));
WSA=transpose(W.data(14,3:qw));
WTAS=transpose(W.data(15,3:qw));
WVIC=transpose(W.data(16,3:qw));
WWA=transpose(W.data(17,3:qw));
WAUS=(WACT+WNSW+WNT+WQLD+WSA+WTAS+WVIC+WWA)/PAUS;
AUS=zeros(DD,1);
for j=1:DD
    AUS(j,1)=max(0,WAUS(j+983,1)-WAUS(j+976,1))/7;
end
UACT=transpose(U.data(10,3:qw));
UNSW=transpose(U.data(11,3:qw));
UNT=transpose(U.data(12,3:qw));
UQLD=transpose(U.data(13,3:qw));
USA=transpose(U.data(14,3:qw));
UTAS=transpose(U.data(15,3:qw));
UVIC=transpose(U.data(16,3:qw));
UWA=transpose(U.data(17,3:qw));
UAUS=(UACT+UNSW+UNT+UQLD+USA+UTAS+UVIC+UWA)/PAUS;
DEAUS=zeros(DD,1);
for j=1:DD
    DEAUS(j,1)=max(0,UAUS(j+983,1)-UAUS(j+976,1))/7;
end

% Western Australia
PWAU=2.69;
WWAU=transpose(W.data(17,3:qw))/PWAU;
WAU=zeros(DD,1);
for j=1:DD
    WAU(j,1)=max(0,WWAU(j+983,1)-WWAU(j+976,1))/7;
end

% Brazil
PBRA=215.019011;
UBRA=transpose(U.data(32,3:qu))/PBRA;
DEBRA=zeros(DD,1);
for j=1:DD
    DEBRA(j,1)=max(0,UBRA(j+983,1)-UBRA(j+976,1))/7;
end

% Brunei Darussalam
PBWN=0.4435;
WBWN=transpose(W.data(33,3:qw))/PBWN;
BWN=zeros(DD,1);
for j=1:DD
    BWN(j,1)=max(0,WBWN(j+983,1)-WBWN(j+976,1))/7;
end
UBWN=transpose(U.data(33,3:qu))/PBWN;
DEBWN=zeros(DD,1);
for j=1:DD
    DEBWN(j,1)=max(0,UBWN(j+983,1)-UBWN(j+976,1))/7;
end

% Hong Kong 
PHKG=7.600852;
WHKG=transpose(W.data(72,3:qw))/PHKG;
HKG=zeros(DD,1);
for j=1:DD
    HKG(j,1)=max(0,WHKG(j+983,1)-WHKG(j+976,1))/7;
end
UHKG=transpose(U.data(72,3:qu))/PHKG;
DEHKG=zeros(DD,1);
for j=1:DD
    DEHKG(j,1)=max(0,UHKG(j+983,1)-UHKG(j+976,1))/7;
end

% France Italy Spain 
PFRA=65.508662;
PITA=60.317073;
PESP=46.784213;
PFIE=PFRA+PITA+PESP;
WFRA=transpose(W.data(132,3:qw));
WITA=transpose(W.data(155,3:qw));
WESP=transpose(W.data(239,3:qw));
WFIE=(WFRA+WITA+WESP)/PFIE;
FIE=zeros(DD,1);
for j=1:DD
    FIE(j,1)=max(0,WFIE(j+983,1)-WFIE(j+976,1))/7;
end
UFRA=transpose(U.data(132,3:qw));
UITA=transpose(U.data(155,3:qw));
UESP=transpose(U.data(239,3:qw));
UFIE=(UFRA+UITA+UESP)/PFIE;
DEFIE=zeros(DD,1);
for j=1:DD
    DEFIE(j,1)=max(0,UFIE(j+983,1)-UFIE(j+976,1))/7;
end

% India
PIND=1402.124607;
WIND=transpose(W.data(149,3:qw))/PIND;
IND=zeros(DD,1);
for j=1:DD
    IND(j,1)=max(0,WIND(j+983,1)-WIND(j+976,1))/7;
end
UIND=transpose(U.data(149,3:qu))/PIND;
DEIND=zeros(DD,1);
for j=1:DD
    DEIND(j,1)=max(0,UIND(j+983,1)-UIND(j+976,1))/7;
end

% Indonesia
PIDN=278.239007;
WIDN=transpose(W.data(150,3:qw))/PIDN;
IDN=zeros(DD,1);
for j=1:DD
    IDN(j,1)=max(0,WIDN(j+983,1)-WIDN(j+976,1))/7;
end
UIDN=transpose(U.data(150,3:qu))/PIDN;
DEIDN=zeros(DD,1);
for j=1:DD
    DEIDN(j,1)=max(0,UIDN(j+983,1)-UIDN(j+976,1))/7;
end

% Japan
PJPN=125.845010;
WJPN=transpose(W.data(157,3:qw))/PJPN;
JPN=zeros(DD,1);
for j=1:DD
    JPN(j,1)=max(0,WJPN(j+983,1)-WJPN(j+976,1))/7;
end
UJPN=transpose(U.data(157,3:qu))/PJPN;
DEJPN=zeros(DD,1);
for j=1:DD
    DEJPN(j,1)=max(0,UJPN(j+983,1)-UJPN(j+976,1))/7;
end

% South Korea
PKOR=51.341022;
WKOR=transpose(W.data(163,3:qw))/PKOR;
KOR=zeros(DD,1);
for j=1:DD
    KOR(j,1)=max(0,WKOR(j+983,1)-WKOR(j+976,1))/7;
end
UKOR=transpose(U.data(163,3:qu))/PKOR;
DEKOR=zeros(DD,1);
for j=1:DD
    DEKOR(j,1)=max(0,UKOR(j+983,1)-UKOR(j+976,1))/7;
end

% Malaysia
PMYS=33.060794;
WMYS=transpose(W.data(179,3:qw))/PMYS;
MYS=zeros(DD,1);
for j=1:DD
    MYS(j,1)=max(0,WMYS(j+983,1)-WMYS(j+976,1))/7;
end
UMYS=transpose(U.data(179,3:qu))/PMYS;
DEMYS=zeros(DD,1);
for j=1:DD
    DEMYS(j,1)=max(0,UMYS(j+983,1)-UMYS(j+976,1))/7;
end

% Mexico
PMEX=131.137507;
UMEX=transpose(U.data(186,3:qu))/PMEX;
DEMEX=zeros(DD,1);
for j=1:DD
    DEMEX(j,1)=max(0,UMEX(j+983,1)-UMEX(j+976,1))/7;
end

% New Zealand
PNZ1=0.017593;
PNZ2=4.892946;
PNZL=PNZ1+PNZ2;
WNZ1=transpose(W.data(203,3:qw));
WNZ2=transpose(W.data(204,3:qw));
WNZL=(WNZ1+WNZ2)/PNZL;
NZL=zeros(DD,1);
for j=1:DD
    NZL(j,1)=max(0,WNZL(j+983,1)-WNZL(j+976,1))/7;
end
XNZ1=transpose(U.data(203,3:qw));
XNZ2=transpose(U.data(204,3:qw));
XNZL=(XNZ1+XNZ2)/PNZL;
DENZL=zeros(DD,1);
for j=1:DD
    DENZL(j,1)=max(0,XNZL(j+983,1)-XNZL(j+976,1))/7;
end

% Philippines
PPHI=112.027348;
WPHI=transpose(W.data(217,3:qw))/PPHI;
PHI=zeros(DD,1);
for j=1:DD
    PHI(j,1)=max(0,WPHI(j+983,1)-WPHI(j+976,1))/7;
end
UIDN=transpose(U.data(217,3:qu))/PIDN;
DEIDN=zeros(DD,1);
for j=1:DD
    DEIDN(j,1)=max(0,UIDN(j+983,1)-UIDN(j+976,1))/7;
end

% Russia
PRUS=146.036343;
WRUS=transpose(W.data(222,3:qw))/PRUS;
RUS=zeros(DD,1);
for j=1:DD
    RUS(j,1)=max(0,WRUS(j+983,1)-WRUS(j+976,1))/7;
end
URUS=transpose(U.data(222,3:qu))/PRUS;
DERUS=zeros(DD,1);
for j=1:DD
    DERUS(j,1)=max(0,URUS(j+983,1)-URUS(j+976,1))/7;
end

% Singapore
PSIN=5.925237;
WSIN=transpose(W.data(235,3:qw))/PSIN;
SIN=zeros(DD,1);
for j=1:DD
    SIN(j,1)=max(0,WSIN(j+983,1)-WSIN(j+976,1))/7;
end
USIN=transpose(U.data(235,3:qu))/PSIN;
DESIN=zeros(DD,1);
for j=1:DD
    DESIN(j,1)=max(0,USIN(j+983,1)-USIN(j+976,1))/7;
end

% Sri Lanka
PLKA=21.559415;
WLKA=transpose(W.data(243,3:qw))/PLKA;
LKA=zeros(DD,1);
for j=1:DD
    LKA(j,1)=max(0,WLKA(j+983,1)-WLKA(j+976,1))/7;
end
ULKA=transpose(U.data(243,3:qu))/PLKA;
DELKA=zeros(DD,1);
for j=1:DD
    DELKA(j,1)=max(0,ULKA(j+983,1)-ULKA(j+976,1))/7;
end

% Taiwan
PTWN=23.61;
WTWN=transpose(W.data(250,3:qw))/PLKA;
TWN=zeros(DD,1);
for j=1:DD
    TWN(j,1)=max(0,WTWN(j+983,1)-WTWN(j+976,1))/7;
end
UTWN=transpose(U.data(250,3:qu))/PTWN;
DETWN=zeros(DD,1);
for j=1:DD
    DETWN(j,1)=max(0,UTWN(j+983,1)-UTWN(j+976,1))/7;
end

% Thailand
PTHA=70.085127;
WTHA=transpose(W.data(253,3:qw))/PTHA;
THA=zeros(DD,1);
for j=1:DD
    THA(j,1)=max(0,WTHA(j+983,1)-WTHA(j+976,1))/7;
end
UTHA=transpose(U.data(253,3:qu))/PTHA;
DETHA=zeros(DD,1);
for j=1:DD
    DETHA(j,1)=max(0,UTHA(j+983,1)-UTHA(j+976,1))/7;
end

% United Staes
PUSA=334.207212;
WUSA=transpose(W.data(261,3:qw))/PUSA;
USA=zeros(DD,1);
for j=1:DD
    USA(j,1)=max(0,WUSA(j+983,1)-WUSA(j+976,1))/7;
end
UUSA=transpose(U.data(261,3:qu))/PUSA;
DEUSA=zeros(DD,1);
for j=1:DD
    DEUSA(j,1)=max(0,UUSA(j+983,1)-UUSA(j+976,1))/7;
end

% United Kingdom
PGBR=68.466544;
WGBR=transpose(W.data(279,3:qw))/PGBR;
GBR=zeros(DD,1);
for j=1:DD
    GBR(j,1)=max(0,WGBR(j+983,1)-WGBR(j+976,1))/7;
end
UGBR=transpose(U.data(279,3:qu))/PGBR;
DEGBR=zeros(DD,1);
for j=1:DD
    DEGBR(j,1)=max(0,UGBR(j+983,1)-UGBR(j+976,1))/7;
end

% Vietnam
PVNM=98.953541;
WVNM=transpose(W.data(284,3:qw))/PVNM;
VNM=zeros(DD,1);
for j=1:DD
    VNM(j,1)=max(0,WVNM(j+983,1)-WVNM(j+976,1))/7;
end
    
websave("./csv/mhlw_cases.csv","https://covid19.mhlw.go.jp/public/opendata/confirmed_cases_cumulative_daily.csv");
websave("./csv/mhlw_deaths.csv","https://covid19.mhlw.go.jp/public/opendata/deaths_cumulative_daily.csv");
XC=importdata("./csv/mhlw_cases.csv");
XD=importdata("./csv/mhlw_deaths.csv");
VC=XC.data(:,:);
VD=XD.data(:,:);

% Okinawa (1.46M): code 47
POKNW=1.469335;    
OKNW1=VC(:,48)/POKNW;
OKNW2=VD(:,48)/POKNW;
OKNW=zeros(DD,1);
for j=1:DD
    OKNW(j,1)=(OKNW1(j+875,1)-OKNW1(j+868,1))/7;
end
DEOKNW=zeros(DD,1);
for j=1:DD
    DEOKNW(j,1)=(OKNW2(j+875,1)-OKNW2(j+868,1))/7;
end

% Tokyo (14M): code 13, 
PTKY=13.988129;
TKY1=VC(:,14)/PTKY;
TKY2=VD(:,28)/PTKY;
TKY=zeros(DD,1);
for j=1:DD
    TKY(j,1)=(TKY1(j+875,1)-TKY1(j+868,1))/7;
end
DETKY=zeros(DD,1);
for j=1:DD
    DETKY(j,1)=(TKY2(j+875,1)-TKY2(j+868,1))/7;
end

% Osaka: code 28
POSK=8.797153;
OSK1=VC(:,28)/POSK;
OSK2=VD(:,28)/POSK;
OSK=zeros(DD,1);
for j=1:DD
    OSK(j,1)=(OSK1(j+875,1)-OSK1(j+868,1))/7;
end
DEOSK=zeros(DD,1);
for j=1:DD
    DEOSK(j,1)=(OSK2(j+875,1)-OSK2(j+868,1))/7;
end

% Hokkaido
PHKD=5.191355;
HKD1=VC(:,2)/PHKD;
HKD2=VD(:,2)/PHKD;
HKD=zeros(DD,1);
for j=1:DD
    HKD(j,1)=(HKD1(j+875,1)-HKD1(j+868,1))/7;
end
DEHKD=zeros(DD,1);
for j=1:DD
    DEHKD(j,1)=(HKD2(j+875,1)-HKD2(j+868,1))/7;
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
         
% plot
subplot(1,2,1)
plot([JPN,TKY,OKNW,OSK,HKD,AUS,KOR,HKG,TWN,SIN],'LineWidth',2)
title('COVID-19: 7-day average of new cases per 1M','data sourced by JHU and MoH of Japan')
xlabel('date');
ylabel('cases/1M');
xticks([1 32 62 DD])
xticklabels({[ll0],[ll1],[ll2],[llf]})
legend('Japan','Tokyo','Okinawa','Osaka','Hokkaido','Australia','South Korea','Hong Kong','Taiwan','Singapore','Location','northwest');
%
subplot(1,2,2)
plot([DEJPN,DETKY,DEOKNW,DEOSK,DEHKD,DEAUS,DEKOR,DEHKG,DETWN,DESIN],'LineWidth',2)
title('COVID-19: 7-day average of deaths per 1M','data sourced by JHU and MoH of Japan')
xlabel('date');
ylabel('deaths/1M');
xticks([1 32 62 DD])
xticklabels({[ll0],[ll1],[ll2],[llf]})
legend('Japan','Tokyo','Okinawa','Osaka','Hokkaido','Australia','South Korea','Hong Kong','Taiwan','Singapore','Location','northwest');
%
set(gcf,'Position',[600,200,1200,400]);
saveas(gcf,'crisis.png');
