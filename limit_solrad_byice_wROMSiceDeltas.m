function [ outtable,newmmtable ] = limit_solrad_byice_wROMSiceDeltas( datesq )
%LIMIT_SOLRAD_BYICE Summary of this function goes here
%
% use output tables from:
% [ solradtable, dailytable ] = era_sol_20151026( filename, latq, lonq )
% and
% dateicetable = ice_timeseries()
%
%   MFILE:   limit_solrad_byice.m
%   MATLAB:  8.6.0.267246 (R2015b)
%   AUTHOR:  Daniel Edward Kaufman (USA)
%            @ The Virginia Institute of Marine Science
%   CONTACT: dkauf42@gmail.com
%   LATEST_UPDATE: October, 2015
%   REVISION HISTORY:
%   - Initial Generation (Oct. 2015)
%

% Before I was using:
% load('iceLimitedSolarRad_Sep2010.mat')

%% Get Data Files and Perform Other Housekeeping Operations
% Get ROMS Output Climate Deltas
ncsourcedir = '~/Desktop/CATEGORIES/CAREER_MANAGEMENT/VIMS/Dissertation/Chapter_2_ForwardModel/climate_scenarios_roms_outputs/';
ncA = [ncsourcedir, 'kaufman.data.future.delta.A.nc'];  iA = ncinfo(ncA);
ncB = [ncsourcedir, 'kaufman.data.future.delta.B.nc'];  iB = ncinfo(ncB);

%%%%%%%%%% COMMON TWEAKS %%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Apply deltas from the ROMS climate scenario outputs %%%%%
% deltasToApply = '2050';
% deltasToApply = '2100';
deltasToApply = 'none';
deltaSite = 'A';  % A: 77.2S  169.5E
% deltaSite = 'B';  % B: 76.5S  176E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Reduce Ice Coverage by a Percentage: %%%%%
% reduceIceByPrctMultiplier = 0.56; % 56% reduction
% reduceIceByPrctMultiplier = 0.78; % 78% reduction
% reduceIceByPrctMultiplier = 1; % 100% reduction
reduceIceByPrctMultiplier = 0; % No reduction
% reduceIceByPrctMultiplier = -0.25; % 25% amplification
% reduceIceByPrctMultiplier = -0.50; % 50% amplification
% reduceIceByPrctMultiplier = -0.56; % 56% amplification
% reduceIceByPrctMultiplier = -0.78; % 78% amplification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Specify the Ice Melt Day: %%%%%
iceMeltDay = '01-Nov-2012'; % this changes the day when light is allowed to start entering the water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Reduce Ice Coverage by fitting it to a new Curve %%%%%
% TransformIce = true;
TransformIce = false;
% whichReduction = '56';
% whichReduction = '78';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOAD DATA
cdir = pwd;
if ~exist('datesq','var')
    [ solradtable, dailytable ] = era_sol_20151026();
    dateicetable = ice_timeseries_areaAvg();
%     dateicetable = ice_timeseries_atPt();
else
    [ solradtable, dailytable ] = era_sol_20151026( datesq );
    dateicetable = ice_timeseries_areaAvg( datesq );
end

cd(cdir)
% Load Ice Data
timedays = dateicetable.date; % in matlabdates, by days
ice = dateicetable.ice;

% Load Solar Radiation Data
timehrs = solradtable.timeLocal; % in matlabdates, by 3 hours
sol = solradtable.ssrd;

%%%%%%%%%%
if TransformIce
    [ ice ] = ice_transformFnc( timedays, ice, whichReduction );
end

%% Get Ice Coverage Deltas for each day from the ROMS scenario outputs

% Time    % days since Sep. 15, 2012, in 5 day increments
timeA = double( ncread(ncA,'time') );
timeB = double( ncread(ncB,'time') );
% AICE: sea ice cover fraction difference
% 73 x 1 single
ice_A2050 = double( ncread(ncA,'aice_2050') );
ice_A2100 = double( ncread(ncA,'aice_2100') );
ice_B2050 = double( ncread(ncB,'aice_2050') );
ice_B2100 = double( ncread(ncB,'aice_2100') );

% Convert ROMS days (from Sep.15,2012) to days-of-year
doyA = date2doy( datenum('Sep 15 2012') + timeA );
doyB = date2doy( datenum('Sep 15 2012') + timeB );

% Determine where the year splits
firsttimestepofnewyear = 23;
lasttimestepoffirstyear = 22;

ice_dV_A2050 = [ ones(1,5) .* ice_A2050(1),...
                 interp1(doyA(1:lasttimestepoffirstyear), ice_A2050(1:lasttimestepoffirstyear), 262:365) ,...
                 interp1([-4.5; 0.5; doyA(firsttimestepofnewyear:end)], ice_A2050(21:end), 1:255),...
                 ice_A2050(end) ];
ice_dV_A2100 = [ ones(1,5) .* ice_A2100(1),...
                 interp1(doyA(1:lasttimestepoffirstyear), ice_A2100(1:lasttimestepoffirstyear), 262:365) ,...
                 interp1([-4.5; 0.5; doyA(firsttimestepofnewyear:end)], ice_A2100(21:end), 1:255),...
                 ice_A2100(end) ];
ice_dV_B2050 = [ ones(1,5) .* ice_B2050(1),...
                 interp1(doyB(1:lasttimestepoffirstyear), ice_B2050(1:lasttimestepoffirstyear), 262:365) ,...
                 interp1([-4.5; 0.5; doyB(firsttimestepofnewyear:end)], ice_B2050(21:end), 1:255),...
                 ice_B2050(end) ];
ice_dV_B2100 = [ ones(1,5) .* ice_B2100(1),...
                 interp1(doyB(1:lasttimestepoffirstyear), ice_B2100(1:lasttimestepoffirstyear), 262:365) ,...
                 interp1([-4.5; 0.5; doyB(firsttimestepofnewyear:end)], ice_B2100(21:end), 1:255),...
                 ice_B2100(end) ];
             
days_ts = 257:621;
days_ts_dates = doy2date(days_ts,repmat(2012,numel(days_ts),1));
days_ts( days_ts>365 ) = days_ts( days_ts>365 )-365;

switch deltasToApply
    case '2050'
        switch deltaSite
            case 'A'
                deltaIce_beforeSort = ice_dV_A2050(:);
                DaysAndIce_deltas = sortrows( [ days_ts(:) deltaIce_beforeSort] );
            case 'B'
                deltaIce_beforeSort = ice_dV_B2050(:);
                DaysAndIce_deltas = sortrows( [ days_ts(:) deltaIce_beforeSort] );
        end
    case '2100'
        switch deltaSite
            case 'A'
                deltaIce_beforeSort = ice_dV_A2100(:);
                DaysAndIce_deltas = sortrows( [ days_ts(:) deltaIce_beforeSort] );
            case 'B'
                deltaIce_beforeSort = ice_dV_B2100(:);
                DaysAndIce_deltas = sortrows( [ days_ts(:) deltaIce_beforeSort] );
        end
    case 'none'
        deltaIce_beforeSort = zeros(365,1);
        DaysAndIce_deltas = [ days_ts(:) deltaIce_beforeSort];
end
deltaDays = days_ts_dates
deltaIce = deltaIce_beforeSort;

days_ts_dates_sep1start = circshift(days_ts_dates,[0,12]);
    daysoffyear = datetime(days_ts_dates_sep1start(1:12),'ConvertFrom','datenum');
    daysoffyear.Year = 2012;
    days_ts_dates_sep1start(1:12) = datenum(daysoffyear);
deltaIce_beforeSort_sep1start = circshift(deltaIce_beforeSort,[12,0]);

% Apply the Deltas
ice = ice + ( deltaIce_beforeSort_sep1start(1:155) .* 100 );
ice(ice > 100) = 100;
ice(ice < 0)   = 0;
putvar(deltaIce_beforeSort_sep1start,days_ts_dates_sep1start)

%% Calculations
outtable = solradtable; % copy solar radiation table to place output in
outtable.ice = nan(length(sol),1);

% Match hours with days (expand hours timeseries to the length(days)), e.g.
% if hours== 1     and days== 1   ,then make days into: 1
%            1.5              2                         1
%            2                3                         2
%            2.5              4                         2
%            3                                          3
%            3.5                                        3
%            4                                          4
%            4.5                                        4
% and expand ice time series accordingly with the days
ntd = numel(timedays);
flooredhrs = floor(timehrs);
for dii = 1:ntd % loop over days in ice table
    thisdayidx = ( flooredhrs == timedays(dii) ); % get day indices in solradtable
    outtable.ice(thisdayidx) = ice(dii);
end


%%%%%% REDUCE ICE COVERAGE BY A PERCENTAGE %%%%%%
outtable.ice = outtable.ice - (outtable.ice .* reduceIceByPrctMultiplier);

% calculate limited solar radiation by the percentage of ice coverage
outtable.solarIced = ( (100-outtable.ice) ./ 100 ) .* abs(sol);


%% Make Figures
% figure;
% plot(outtable.timeLocal, outtable.ssrd);
% datetick('x',2,'keeplimits')
% 
% figure;
% plot(outtable.timeLocal, outtable.solarIced);
% datetick('x',2,'keeplimits')

%% Visualize both Timeseries for the model time period
% XL = [min(outtable.timeLocal)-1 max(outtable.timeLocal)];
% FH = figure('color','w','position',[651   498   888   830]);
% spm = 3;  spn = 1;
% h(1)=subplot(spm,spn,1);
%     plot(outtable.timeLocal,outtable.ice,'marker','o');
%     xlim(XL); datetick('x','mm/dd/yyyy','keeplimits'); grid on;
%     ylabel('% Ice Coverage')
% %     dateaxis('x',2)
% h(2)=subplot(spm,spn,2);
%     plot(outtable.timeLocal,outtable.ssrd,'marker','o');
%     xlim(XL); datetick('x','mm/dd/yyyy','keeplimits'); grid on;
%     ylabel(sprintf('Surface Solar Radiation\nDownwards (W m^{-2})'));
% %     dateaxis('x',2)
% h(3)=subplot(spm,spn,3);
%     plot(outtable.timeLocal,outtable.solarIced,'marker','o');
%     xlim(XL); datetick('x','mm/dd/yyyy','keeplimits'); grid on;
%     ylabel(sprintf('Ice Limited Solar\nRadiation (W m^{-2})'));
% %     dateaxis('x',2)
% figFormat(FH);
% axs=findall(gcf,'type','axes');
% linkaxes(axs,'x');

%% For MarMOT Table

% Get Days for MarMOT
% marmot_base_DoY = 244; % day before september 1, 2012 (in leap year)
marmot_base_DoY = 245; % september 1, 2012 (in leap year)
daysInYear = 366; % 2012 is a leap year
daydif = timehrs(1) - marmot_base_DoY;
marmotT = timehrs(:) - daydif;

% KEEP solrad 0 until Nov. 15, which we are defining as 'iceMeltDay'
imdIdx = find(doy2date(marmotT(:),ones(size(marmotT(:)))*2012)==datenum(iceMeltDay));
outtable.solarIced( 1 : imdIdx-1 ) = 0;

mmtable = readtable('ft_base2012_ncepIceLim_sol_rsA.dat', 'delimiter', ' ');
mmsite = repmat(mmtable.site(1),numel(marmotT),1);

newmmtable = table(mmsite,marmotT,round(outtable.solarIced,4),'VariableNames',{'site','t','sol'});

end

