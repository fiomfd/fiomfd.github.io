websave("./csv/mhlw_cases.csv","https://covid19.mhlw.go.jp/public/opendata/confirmed_cases_cumulative_daily.csv");
websave("./csv/mhlw_deaths.csv","https://covid19.mhlw.go.jp/public/opendata/deaths_cumulative_daily.csv");
XC=importdata("./csv/mhlw_cases.csv");
XD=importdata("./csv/mhlw_deaths.csv");
l=length(XC.data(:,1));
VC=XC.data(:,:);
VD=XD.data(:,:);

d0=datetime('2022-10-01');
d1=datetime('2022-11-01');
d2=datetime('2022-12-01');
d3=datetime('2023-01-01');

xdf=strrep(XC.textdata(l+1,1),'/','-');
df=datetime(xdf);
D=1+days(df-d0);
l0=datestr(d0,'yyyy-mm-dd');
l1=datestr(d1,'yyyy-mm-dd');
l2=datestr(d2,'yyyy-mm-dd');
l3=datestr(d3,'yyyy-mm-dd');
lf=datestr(df,'yyyy-mm-dd');

% Japan
PJPN=125.845010;
CJPN=VC(:,1)/PJPN;
DJPN=VD(:,1)/PJPN;
NJPN=zeros(D,1);
for j=1:D
    NJPN(j,1)=(CJPN(j+875,1)-CJPN(j+868,1))/7;
end
NDJPN=zeros(D,1);
for j=1:D
    NDJPN(j,1)=(DJPN(j+875,1)-DJPN(j+868,1))/7;
end

% Tokyo
PTKY=13.988129;
CTKY=VC(:,14)/PTKY;
DTKY=VD(:,14)/PTKY;
NTKY=zeros(D,1);
for j=1:D
    NTKY(j,1)=(CTKY(j+875,1)-CTKY(j+868,1))/7;
end
NDTKY=zeros(D,1);
for j=1:D
    NDTKY(j,1)=(DTKY(j+875,1)-DTKY(j+868,1))/7;
end

% Okinawa
POKNW=1.469335;    
COKNW=VC(:,48)/POKNW;
DOKNW=VD(:,48)/POKNW;
NOKNW=zeros(D,1);
for j=1:D
    NOKNW(j,1)=(COKNW(j+875,1)-COKNW(j+868,1))/7;
end
NDOKNW=zeros(D,1);
for j=1:D
    NDOKNW(j,1)=(DOKNW(j+875,1)-DOKNW(j+868,1))/7;
end

% Osaka
POSK=8.797153;
COSK=VC(:,28)/POSK;
DOSK=VD(:,28)/POSK;
NOSK=zeros(D,1);
for j=1:D
    NOSK(j,1)=(COSK(j+875,1)-COSK(j+868,1))/7;
end
NDOSK=zeros(D,1);
for j=1:D
    NDOSK(j,1)=(DOSK(j+875,1)-DOSK(j+868,1))/7;
end

% Hokkaido
PHKD=5.191355;
CHKD=VC(:,2)/PHKD;
DHKD=VD(:,2)/PHKD;
NHKD=zeros(D,1);
for j=1:D
    NHKD(j,1)=(CHKD(j+875,1)-CHKD(j+868,1))/7;
end
NDHKD=zeros(D,1);
for j=1:D
    NDHKD(j,1)=(DHKD(j+875,1)-DHKD(j+868,1))/7;
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
subplot(1,2,1)
plot([NJPN NTKY NOKNW NOSK NHKD],'LineWidth',1.5)
title('COVID-19: 百万人あたりの日毎感染者数の7日間平均','厚生労働省の公開データに基づく')
xlabel('日付')
ylabel('百万人あたりの人数')
xticks([1 32 62 D])
xticklabels({[l0],[l1],[l2],[lf]})
legend('日本国','東京都','沖縄県','大阪府','北海道','Location','northwest');
%
subplot(1,2,2)
plot([NDJPN NDTKY NDOKNW NDOSK NDHKD],'LineWidth',1.5)
title('COVID-19: 百万人あたりの死亡者数の7日間平均','厚生労働省の公開データに基づく')
xlabel('日付')
ylabel('百万人あたりの人数')
xticks([1 32 62 D])
xticklabels({[l0],[l1],[l2],[lf]})
legend('日本国','東京都','沖縄県','大阪府','北海道','Location','northwest');
%
set(gcf,'Position',[500,200,1400,500]);
saveas(gcf,'xmhlwj.png');