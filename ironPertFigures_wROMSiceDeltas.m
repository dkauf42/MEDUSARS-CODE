% SCRIPT: ironPertFigures
%
% run after ironPertTimeSeries.m

[ ft, daysOfYr, ice, deltaice, feLvls ] = ironPertTimeSeries_wROMSiceDeltas;

dates = doy2date(daysOfYr,repmat(2012,length(daysOfYr),1));
XL = [min(dates) max(dates)]; % XL = [258 396];
xlab = '';

% Iron Input
figure;
plot(dates,ft.pertdfe(ft.k==1 & ft.t<399));
xlim(XL); xlabel(xlab); grid on;
ylabel('dfe (mmol/m3)')
datetick('x','mm/dd/yyyy','keeplimits')

% Ice Cover Time Series
figure;
plot(dates,ice);
xlim(XL); xlabel(xlab); grid on;
ylabel('ice %')
datetick('x','mm/dd/yyyy','keeplimits')

% Delta-ice time series
figure;
plot(dates,deltaice);
xlim(XL); xlabel(xlab); grid on;
ylabel('delta ice')
datetick('x','mm/dd/yyyy','keeplimits')

% lowest dFe level
figure;
plot(dates,feLvls);
yl=ylim; ylim([0 yl(2)]);
xlim(XL); xlabel(xlab); grid on;
ylabel('lowest dfe level')
datetick('x','mm/dd/yyyy','keeplimits')


% Clean Up
figFormat('all');