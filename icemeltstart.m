% SCRIPT:  icemeltstart
%
%  REMEMBER, to set date range in ice_timeseries.m and era_sol_20151026.m
%            to 2010 to '01 September 2010' : '2 February 2014'
%

[ outtable,newmmtable ] = limit_solrad_byice_wROMSiceDeltas( datenum('01 September 2012'):datenum('2 February 2013') );
clearvars('-except','outtable');

its=outtable.ice;

%%
% Try: Find Simple 80% Cross Points
crossPt = 80;
itsAbove=its>crossPt; itsBelow=its<crossPt; % points ABOVE and BELOW 80% ice cover
xts=itsBelow(2:end) - itsAbove(1:end-1);

xpts1 = xts==0; % points where ice level crosses 80% 
dtx1 = datetime(outtable.timeLocal(xpts1),'ConvertFrom','datenum'); % datetimes of cross points
dtx1( dtx1.Month < 9 )=[]; % remove cross points that occur before september each year

%%
% Try: Find Simple 80% Cross Points WITH a subsequent week(7 days) of below-80% ice cover
crossPt = 80;
itsAbove=its>crossPt; itsBelow=its<crossPt; % points ABOVE and BELOW 80% ice cover
xts=itsBelow(2:end) - itsAbove(1:end-1);

xpts2 = findPattern2( xts, [ 0 1 1 1 1 1 1 1 1 1 1 ] );
dtx2 = datetime(outtable.timeLocal(xpts2),'ConvertFrom','datenum'); % datetimes of cross points
dtx2( dtx2.Month < 9 )=[]; % remove cross points that occur before september each year

%%
% Try: Find Simple 80% Cross Points, and where at least 75% of the next 30 days are
     % still below 80% ice cover  ("sustains low ice for >= 75% of subsequent time window")
crossPtforRangePlot = 75; % 80 (Cross Point)
loopcount=1;
xts=[];
xpts3=[];
earlyDays = datetime([],[],[]);
crossPtRange = 70:1:95;
for crossPtforRangePlot=crossPtRange
itsAbove=its>crossPtforRangePlot; itsBelow=its<crossPtforRangePlot; % points ABOVE and BELOW 80% ice cover
xts=itsBelow(2:end) - itsAbove(1:end-1);

xpts3 = find(xts==0) + 1; % points where ice level crosses 80% 
windowL = 240; % time steps are 3 hours, so there are 8 steps in a day; e.g. for a 30 day window, enter 240 (240=30*8)
prctThrsh = 0.75; % 0.75 % in decimal form (at least this prct of the next days in window (e.g. 30)
     % Next remove points that are too close to the end to possibly match the pattern.
windowEnds = xpts3+windowL-1;
xpts3(windowEnds>length(its)) = [];
k=1;
for xii=1:numel(xpts3)
    twin = xts( xpts3(xii)+1 : xpts3(xii)+windowL);
    numBelow = sum(twin == 1);
    prctBelow(xii) = numBelow / windowL;
    if prctBelow(xii) >= prctThrsh
        newxpts3(k,1) = xpts3(xii);
        k=k+1;
    end
end
dtx3 = datetime(outtable.timeLocal(newxpts3),'ConvertFrom','datenum'); % datetimes of cross points
dtx3( dtx3.Month < 9 )=[]; % remove cross points that occur before september each year
dtx3 = datetime(unique([dtx3.Year dtx3.Month dtx3.Day],'rows'));
[c,i]=unique(dtx3.Year); % find first days (after september) for each year

earlyDays(:,loopcount)=dtx3(i);
loopcount=loopcount+1;

if crossPtforRangePlot==80
    fprintf('80 cross point\n');
    disp(earlyDays(:,loopcount-1))
end
end

%% Make Figures
FH = figure;
co=get(groot,'defaultaxescolororder');
% subplot(4,1,1);
ax = gca;
    yrstr = '2012';
    plot(crossPtRange,datenum(earlyDays(1,:)),'linestyle','none','marker','o','markerfacecolor',co(1,:))
    ylim([ datenum(['01-Sep-',yrstr]) datenum(['31-Dec-',yrstr]) ])
    yl=ylim; ax.YTick=yl(1):10:yl(end);
    ax.YMinorTick='on';
    ax.YMinorGrid='on';   grid on;
    datetick('y','mmm dd','keepticks','keeplimits')
    hold all
    hline(datenum(['10-Nov-',yrstr]),'r-');
    ax.XDir='reverse';
%     title(sprintf(['Dates When Ice Drops Below Threshold \n',...
%         'and Remains Below for ',num2str(prctThrsh*100),...
%         '% of the next ',num2str(windowL/8),' days']));
%     xlabel('Ice Threshold');
%     ylabel(['Dates - ',yrstr])
% % % subplot(4,1,2);
% % % ax = gca;
% % %     yrstr = '2011';
% % %     plot(crossPtRange,datenum(earlyDays(2,:)),'linestyle','none','marker','o','markerfacecolor',co(1,:))
% % %     ylim([ datenum(['01-Oct-',yrstr]) datenum(['31-Dec-',yrstr]) ])
% % %     yl=ylim; ax.YTick=yl(1):10:yl(end);
% % %     ax.YMinorTick='on';
% % %     ax.YMinorGrid='on';   grid on;
% % %     datetick('y','mmm dd','keepticks','keeplimits')
% % %     hold all
% % %     hline(datenum(['10-Nov-',yrstr]),'r-');
% % %     ax.XDir='reverse';
% % % subplot(4,1,3);
% % % ax = gca;
% % %     yrstr = '2012';
% % %     plot(crossPtRange,datenum(earlyDays(3,:)),'linestyle','none','marker','o','markerfacecolor',co(1,:))
% % %     ylim([ datenum(['01-Oct-',yrstr]) datenum(['31-Dec-',yrstr]) ])
% % %     yl=ylim; ax.YTick=yl(1):10:yl(end);
% % %     ax.YMinorTick='on';
% % %     ax.YMinorGrid='on';   grid on;
% % %     datetick('y','mmm dd','keepticks','keeplimits')
% % %     hold all
% % %     hline(datenum(['10-Nov-',yrstr]),'r-');
% % %     ax.XDir='reverse';
% % % subplot(4,1,4);
% % % ax = gca;
% % %     yrstr = '2013';
% % %     plot(crossPtRange,datenum(earlyDays(4,:)),'linestyle','none','marker','o','markerfacecolor',co(1,:))
% % %     ylim([ datenum(['01-Oct-',yrstr]) datenum(['31-Dec-',yrstr]) ])
% % %     yl=ylim; ax.YTick=yl(1):10:yl(end);
% % %     ax.YMinorTick='on';
% % %     ax.YMinorGrid='on';   grid on;
% % %     datetick('y','mmm dd','keepticks','keeplimits')
% % %     hold all
% % %     hline(datenum(['10-Nov-',yrstr]),'r-');
% % %     ax.XDir='reverse';
% % % %     xlabel('% Ice Cover Threshold')
    
figFormat(FH)