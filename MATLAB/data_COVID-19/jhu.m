websave('./csv/jhu_cases.csv','https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv');
websave('./csv/jhu_deaths.csv','https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv');
% websave('./csv/jhu_discharged.csv','https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv');

A=importdata('./csv/jhu_cases.csv');
B=importdata('./csv/jhu_deaths.csv');
% C=importdata('./csv/jhu_discharged.csv');
[pw qw]=size(A.data);

D=qw-2;
d0=datetime(2020,1,22);
d1=d0+days(D-1);
l0=datestr(d0,'yyyy-mm-dd');
l1=datestr(d0+days(floor(D/3)),'yyyy-mm-dd');
l2=datestr(d0+days(floor(2*D/3)),'yyyy-mm-dd');
l3=datestr(d1,'yyyy-mm-dd');

DD=qw-528;
dd0=datetime('2021-07-01');
dd1=dd0+days(DD-1);
ll0=datestr(dd0,'yyyy-mm-dd');
ll1=datestr(dd0+days(floor((DD-1)/2)),'yyyy-mm-dd');
ll2=datestr(dd1,'yyyy-mm-dd');

% Argentina
PARG=45.870820;
AARG=transpose(A.data(8,3:qw))/PARG;
BARG=transpose(B.data(8,3:qw))/PARG;
NARG=zeros(D,1);
for j=1:7
    NARG(j,1)=AARG(j,1)/j;
end
for j=8:D
    NARG(j,1)=max(0,AARG(j,1)-AARG(j-7,1))/7;
end
NDARG=zeros(D,1);
for j=1:7
    NDARG(j,1)=BARG(j,1)/7;
end
for j=8:D
    NDARG(j,1)=(BARG(j,1)-BARG(j-7,1))/7;
end

% Brazil
PBRA=215.019011;
ABRA=transpose(A.data(32,3:qw))/PBRA;
BBRA=transpose(B.data(32,3:qw))/PBRA;
NBRA=zeros(D,1);
for j=1:7
    NBRA(j,1)=ABRA(j,1)/j;
end
for j=8:D
    NBRA(j,1)=max(0,ABRA(j,1)-ABRA(j-7,1))/7;
end
NDBRA=zeros(D,1);
for j=1:7
    NDBRA(j,1)=BBRA(j,1)/7;
end
for j=8:D
    NDBRA(j,1)=(BBRA(j,1)-BBRA(j-7,1))/7;
end


% Colombia
PCOL=51.765589;
ACOL=transpose(A.data(94,3:qw))/PCOL;
BCOL=transpose(B.data(94,3:qw))/PCOL;
NCOL=zeros(D,1);
for j=1:7
    NCOL(j,1)=ACOL(j,1)/j;
end
for j=8:D
    NCOL(j,1)=max(0,ACOL(j,1)-ACOL(j-7,1))/7;
end
NDCOL=zeros(D,1);
for j=1:7
    NDCOL(j,1)=BCOL(j,1)/7;
end
for j=8:D
    NDCOL(j,1)=max(0,(BCOL(j,1)-BCOL(j-7,1)))/7;
end

% France Italy Spwin
PFRA=65.508662;
PITA=60.317073;
PESP=46.784213;
PFIE=PFRA+PITA+PESP;
AFRA=transpose(A.data(132,3:qw));
AITA=transpose(A.data(155,3:qw));
AESP=transpose(A.data(242,3:qw));
AFIE=(AFRA+AITA+AESP)/PFIE;
BFRA=transpose(B.data(132,3:qw));
BITA=transpose(B.data(155,3:qw));
BESP=transpose(B.data(242,3:qw));
BFIE=(BFRA+BITA+BESP)/PFIE;
NFIE=zeros(D,1);
for j=1:7
    NFIE(j,1)=AFIE(j,1)/j;
end
for j=8:D
    NFIE(j,1)=max(0,AFIE(j,1)-AFIE(j-7,1))/7;
end
NDFIE=zeros(D,1);
for j=1:7
    NDFIE(j,1)=BFIE(j,1)/7;
end
for j=8:D
    NDFIE(j,1)=max(0,(BFIE(j,1)-BFIE(j-7,1)))/7;
end

% Mexico
PMEX=131.137507;
AMEX=transpose(A.data(186,3:qw))/PMEX;
BMEX=transpose(B.data(186,3:qw))/PMEX;
NMEX=zeros(D,1);
for j=1:7
    NMEX(j,1)=AMEX(j,1)/j;
end
for j=8:D
    NMEX(j,1)=max(0,AMEX(j,1)-AMEX(j-7,1))/7;
end
NDMEX=zeros(D,1);
for j=1:7
    NDMEX(j,1)=BMEX(j,1)/7;
end
for j=8:D
    NDMEX(j,1)=max(0,(BMEX(j,1)-BMEX(j-7,1)))/7;
end

% United Staes
PUSA=334.207212;
AUSA=transpose(A.data(261,3:qw))/PUSA;
BUSA=transpose(B.data(261,3:qw))/PUSA;
NUSA=zeros(D,1);
for j=1:7
    NUSA(j,1)=AUSA(j,1)/j;
end
for j=8:D
    NUSA(j,1)=max(0,AUSA(j,1)-AUSA(j-7,1))/7;
end
NDUSA=zeros(D,1);
for j=1:7
    NDUSA(j,1)=BUSA(j,1)/7;
end
for j=8:D
    NDUSA(j,1)=max(0,(BUSA(j,1)-BUSA(j-7,1)))/7;
end

% United Kingdom
PGBR=68.466544;
AGBR=transpose(A.data(279,3:qw))/PGBR;
BGBR=transpose(B.data(279,3:qw))/PGBR;
NGBR=zeros(D,1);
for j=1:7
    NGBR(j,1)=AGBR(j,1)/j;
end
for j=8:D
    NGBR(j,1)=max(0,AGBR(j,1)-AGBR(j-7,1))/7;
end
NDGBR=zeros(D,1);
for j=1:7
    NDGBR(j,1)=BGBR(j,1)/7;
end
for j=8:D
    NDGBR(j,1)=max(0,(BGBR(j,1)-BGBR(j-7,1)))/7;
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
subplot(2,2,1)
plot([AUSA,AGBR,AARG,ABRA,ACOL,AMEX,AFIE],'LineWidth',2)
title('COVID-19: cases per 1M','data sourced by JHU Coronavirus Resource Center')
xlabel('date');
ylabel('cases/1M');
xticks([1 floor(D/3) floor(2*D/3) D]);
xticklabels({[l0],[l1],[l2],[l3]});
legend('United States','United Kingdom','Argentina','Brazil','Colombia','Mexico','FRA+ITA+ESP','Location','northwest');
% plot
subplot(2,2,2)
plot([BUSA,BGBR,BARG,BBRA,BCOL,BMEX,BFIE],'LineWidth',2)
title('COVID-19: deaths per 1M','data sourced by JHU Coronavirus Resource Center')
xlabel('date');
ylabel('deaths/1M');
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[l3]})
legend('United States','United Kingdom','Argentina','Brazil','Colombia','Mexico','FRA+ITA+ESP','Location','northwest');
% plot
subplot(2,2,3)
plot([NUSA,NGBR,NARG,NBRA,NCOL,NMEX,NFIE],'LineWidth',2)
title('COVID-19: 7-day average of new cases per 1M','data sourced by JHU Coronavirus Resource Center')
xlabel('date');
ylabel('cases/1M');
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[l3]})
legend('United States','United Kingdom','Argentina','Brazil','Colombia','Mexico','FRA+ITA+ESP','Location','northwest');
% plot
subplot(2,2,4)
plot([NDUSA,NDGBR,NDARG,NDBRA,NDCOL,NDMEX,NDFIE],'LineWidth',2)
title('COVID-19: 7-day average deaths per 1M','data sourced by JHU Coronavirus Resource Center')
xlabel('date');
ylabel('deaths/1M');
xticks([1 floor(D/3) floor(2*D/3) D])
xticklabels({[l0],[l1],[l2],[l3]})
legend('United States','United Kingdom','Argentina','Brazil','Colombia','Mexico','FRA+ITA+ESP','Location','northeast');
set(gcf,'Position',[600,200,1200,800]);
saveas(gcf,'jhu.png');
