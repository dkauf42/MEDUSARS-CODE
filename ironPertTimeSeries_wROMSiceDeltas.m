function [ ft, daysOfYr, ice, deltaice, feLvls ] = ironPertTimeSeries_wROMSiceDeltas()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Specify the Ice Melt Day: %%%%%
day0 = date2doy(datenum('01-Nov-2012')); % this changes the day when iron is allowed to start entering the water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a template table from file for inserting the time series into
[ forcingTable ] = temperatureTimeSeriesFromGliderMLDwROMSdeltas;
ft = forcingTable;
clearvars('-except','ft');
% Reformat the Table Variables
ft(:,5)=[];
ft.Properties.VariableNames{4} = 'pertdfe';

% Get Ice Time Series
[ outtable, ~ ] = limit_solrad_byice_wROMSiceDeltas( );
daysOfYr = date2doy(outtable.timeLocal);
daysOfYr(daysOfYr < 40) = daysOfYr(daysOfYr < 40) + 366;
inDayRng = daysOfYr > min(ft.t) & daysOfYr < max(ft.t);
daysOfYr = daysOfYr(inDayRng);
ice = outtable.ice(inDayRng); 
[ud,di,~]=unique(floor(daysOfYr));
daysOfYr = ud;
ice = ice(di);
nd = numel(daysOfYr);
clear inDayRng ud di

% Assign Iron Levels and Release From Sea Ice
ice0 = 100; % start with 100% ice coverage regardless
fe0 = 0.00116;  % starting amount in mmol/m3, Sea Ice Fe Source of WESTERN Ross Sea (from McGillicuddy et al 2015, Table 1)

pertdfe = nan(nd,1); deltaice = nan(nd,1);
feLvls = nan(nd,1);  fel = fe0;
lowIceLvl = ice0;
for iii = 1:nd % loop over days
    if daysOfYr(iii) >= day0 % don't start adding until after day number 'day0'
        if ice(iii) < lowIceLvl % this day hits a new low, so iron is released
            deltaice(iii) = lowIceLvl - ice(iii);
            pertdfe(iii) = (deltaice(iii)/100) * fe0;
            
            lowIceLvl = ice(iii); % set the new lowest ice level
            feLvls(iii) = fel - pertdfe(iii); fel = feLvls(iii); % track felvl
        else
            deltaice(iii) = 0; pertdfe(iii) = 0; feLvls(iii) = fel;
        end
    else
            deltaice(iii) = 0; pertdfe(iii) = 0; feLvls(iii) = fel;
    end
    
    
    surfForDay = (ft.t==daysOfYr(iii) & ft.k==1); % get each day's surface row index
    ft.pertdfe(surfForDay) = roundsd(pertdfe(iii),4); % put iron release into output table
end

% Add 0.1 umol Fe m-2 yr-1 (5.4e-8 mmol m-3 d-1 over a 5m surface layer)
%     from discussion pg 8, McGillicuddy eet al 2015, referencing Edwards
%     and Sedwick, 2001.
% ft.pertdfe(ft.k==1) = ft.pertdfe(ft.k==1) + 0.000000054;

% % % % idxIceMelt = ft.t > 326  &  ft.t < 366; % 326 is Nov 21 for leapyear, 
% % % % ft.pertdfe( ft.k==1  &  idxIceMelt) = 0.0000032;


return % END FUNCTION