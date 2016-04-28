function [forcingTable] = temperatureTimeSeriesFromGliderMLDwROMSdeltas(varargin)

%% Initial Function Housekeeping 
    p=inputParser;
    addParameter(p,'makefig',false,@islogical);
parse(p,varargin{:});
makefig = p.Results.makefig;

%% Get Data Files and Perform Housekeeping Operations
sourcedir = '~/Desktop/CATEGORIES/CAREER_MANAGEMENT/VIMS/Dissertation/Chapter_2_ForwardModel/Glider_Data_ALL/RMJ_Glider_2012-2013_Data_Structures/';
load([sourcedir,'gdivegcast_SG503_2012.mat'],'gcast'); % get 'gcast' variable (structure type)
forcingTable=marmot_tabread([sourcedir,'fkt_base2010_roms_vdcFromMLDandLlort_rsA.dat']);

gcast = gcast(cell2mat({gcast(:).used}) ~= 0,:); % remove 'unused' data rows

% Get ROMS Output Climate Deltas
ncsourcedir = '~/Desktop/CATEGORIES/CAREER_MANAGEMENT/VIMS/Dissertation/Chapter_2_ForwardModel/climate_scenarios_roms_outputs/';
ncA = [ncsourcedir, 'kaufman.data.future.delta.A.nc'];  iA = ncinfo(ncA);
ncB = [ncsourcedir, 'kaufman.data.future.delta.B.nc'];  iB = ncinfo(ncB);

%%%%%%%%%% COMMON TWEAKS %%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Apply deltas from the ROMS climate scenario outputs %%%%%
% deltasToApply = '2050';
deltasToApply = '2100';
% deltasToApply = 'none';
deltaSite = 'A';  % A: 77.2S  169.5E
% deltaSite = 'B';  % B: 76.5S  176E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Make Temperatures from the glider higher by degrees: %%%%%%%%%
% higherTempsByDeg = 1.66; % 1.66 deg C Increase
% higherTempsByDeg = 2.50; % 2.50 deg C Increase
higherTempsByDeg = 4.0; % 4 deg C Increase
% higherTempsByDeg = 0; % No Increase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Get Temperature Profiles for each day in the Glider Data Table
nc = numel(gcast);
posCounter = 1;
for cii = 1:nc
    thisDaten = gcast(cii).datetimeno;
    castLength = numel(thisDaten);
    nextPos = posCounter + castLength;
    
    daten(posCounter:nextPos-1) = gcast(cii).datetimeno;
    tempc(posCounter:nextPos-1) = gcast(cii).tempc;
    depthm(posCounter:nextPos-1) = gcast(cii).depthm;
    
    posCounter = nextPos;
end

% Get Indices for each day in glider data table
daysVector = floor( daten(1) ) : 1 : ceil( daten(end) );
[~,binEdge,binIdx]=histcounts(daten, daysVector);
binMid = (binEdge(1:end-1) + binEdge(2:end) ) ./ 2;

% Convert glider dates to DaysOfYear
GliderDoYs = date2doy( floor(binMid(:)) );

% remove leap year numbers (glider was from 2012, a leap year; so, dec 31
%   is day #366.  But we just want a year of days 1 to 365
GliderDoYs(1:40) = GliderDoYs(1:40) - 1;

clear posCounter thisDaten castLength cii gcast daysVector nextPos binMid

%% Go through each day in glider file...
%   and interpolate temperatures to model k levels
% K to Z Transformation
klev = (1:41).'; nk = numel(klev);
zmid = (2.5:5:202.5).';

nd = numel(GliderDoYs);
placecounter = 1;
tempMat = nan(nk,365);
for dii = 1:nd % loop over the glider days
    
    thisDayIdx = binIdx==dii;
    thisDoY = GliderDoYs(dii);
    thisDoY = mod(thisDoY,365);
    thisDoY( thisDoY == 0) = 365;
    
    thisDayZ = depthm(thisDayIdx);
    thisDayT = tempc(thisDayIdx);
    
    nnz = ~isnan(thisDayZ);
    nnt = ~isnan(thisDayT);
    thisDayTK = interp1(thisDayZ(nnz & nnt),thisDayT(nnz & nnt),zmid);
    
    if dii==1
        firstTK = thisDayTK;
    elseif dii==nd
        lastTK = thisDayTK;
    end
    
    tempMat(klev,thisDoY) = thisDayTK;
    
    placecounter = placecounter + nk;
end

clear nd thisDoY placecounter dii

% % % % %% Interpolate Glider Temps for rest of year (between day 39 and day 327)
% % % % winterDays = 40:325;
% % % % winterTemps = interp2([39 326], klev, [lastTK(:), firstTK(:)], [40:325], klev);
% % % % % nw = numel(winterDays);

%% Copy the first glider profile to the rest of the year
winterDays = 40:325;
winterTemps = interp2([39 326], klev, [firstTK(:), firstTK(:)], 40 : 325, klev);

tempMat(klev,winterDays) = winterTemps;
% fulltempvec = tempMat(:);

clear winterDays winterTemps

%% Get Temperatures for each day from the ROMS scenario outputs

% Time    % days since Sep. 15, 2012, in 5 day increments
timeA = double( ncread(ncA,'time') );
timeB = double( ncread(ncB,'time') );
zrA   = -ncread(ncA,'zr');  % (zr = depths of center of each 24 vertical layer
zrB   = -ncread(ncB,'zr');  % (zr = depths of center of each 24 vertical layer
%  difference values given are future - 20th century
% 24 x 73 single
temp_A2050 = double( ncread(ncA,'temp_2050') );
temp_A2100 = double( ncread(ncA,'temp_2100') );
temp_B2050 = double( ncread(ncB,'temp_2050') );
temp_B2100 = double( ncread(ncB,'temp_2100') );

% Convert ROMS days (from Sep.15,2012) to days-of-year
doyA = date2doy( datenum('Sep 15 2012') + timeA );
doyB = date2doy( datenum('Sep 15 2012') + timeB );

% Determine where the year splits
firsttimestepofnewyear = 23;
lasttimestepoffirstyear = 22;
             
size(temp_A2050)

% Interpolate to a full year
temp_dV_A2050_beforeZinterp = [ ones(24,5) .* repmat(temp_A2050(:,1),1,5),...
                 interp1(doyA(1:lasttimestepoffirstyear), temp_A2050(:,1:lasttimestepoffirstyear).', 262:365).' ,...
                 interp1([-4.5; 0.5; doyA(firsttimestepofnewyear:end)], temp_A2050(:,21:end).', 1:255).',...
                 temp_A2050(:,end) ];
temp_dV_A2100_beforeZinterp = [ ones(24,5) .* repmat(temp_A2100(:,1),1,5),...
                 interp1(doyA(1:lasttimestepoffirstyear), temp_A2100(:,1:lasttimestepoffirstyear).', 262:365).' ,...
                 interp1([-4.5; 0.5; doyA(firsttimestepofnewyear:end)], temp_A2100(:,21:end).', 1:255).',...
                 temp_A2100(:,end) ];
temp_dV_B2050_beforeZinterp = [ ones(24,5) .* repmat(temp_B2050(:,1),1,5),...
                 interp1(doyB(1:lasttimestepoffirstyear), temp_B2050(:,1:lasttimestepoffirstyear).', 262:365).' ,...
                 interp1([-4.5; 0.5; doyB(firsttimestepofnewyear:end)], temp_B2050(:,21:end).', 1:255).',...
                 temp_B2050(:,end) ];
temp_dV_B2100_beforeZinterp = [ ones(24,5) .* repmat(temp_B2100(:,1),1,5),...
                 interp1(doyB(1:lasttimestepoffirstyear), temp_B2100(:,1:lasttimestepoffirstyear).', 262:365).' ,...
                 interp1([-4.5; 0.5; doyB(firsttimestepofnewyear:end)], temp_B2100(:,21:end).', 1:255).',...
                 temp_B2100(:,end) ];
size(temp_dV_B2100_beforeZinterp)

% Interpolate to the k level depths
temp_dV_A2050 = interp1( zrA, temp_dV_A2050_beforeZinterp, zmid, 'linear', 'extrap');
temp_dV_A2100 = interp1( zrA, temp_dV_A2100_beforeZinterp, zmid, 'linear', 'extrap');
temp_dV_B2050 = interp1( zrA, temp_dV_B2050_beforeZinterp, zmid, 'linear', 'extrap');
temp_dV_B2100 = interp1( zrA, temp_dV_B2100_beforeZinterp, zmid, 'linear', 'extrap');
size(temp_dV_B2100)

days_ts = 257:621;
days_ts_dates = doy2date(days_ts,repmat(2012,numel(days_ts),1));
days_ts( days_ts>365 ) = days_ts( days_ts>365 )-365;

switch deltasToApply
    case '2050'
        switch deltaSite
            case 'A'
                deltaTemps_beforeSort = temp_dV_A2050.';
                DaysAndTemps_deltas = sortrows( [ days_ts(:) deltaTemps_beforeSort] );
            case 'B'
                deltaTemps_beforeSort = temp_dV_B2050.';
                DaysAndTemps_deltas = sortrows( [ days_ts(:) deltaTemps_beforeSort] );
        end
    case '2100'
        switch deltaSite
            case 'A'
                deltaTemps_beforeSort = temp_dV_A2100.';
                DaysAndTemps_deltas = sortrows( [ days_ts(:) deltaTemps_beforeSort] );
            case 'B'
                deltaTemps_beforeSort = temp_dV_B2100.';
                DaysAndTemps_deltas = sortrows( [ days_ts(:) deltaTemps_beforeSort] );
        end
    case 'none'
        deltaTemps_beforeSort = zeros(365,nk);
        DaysAndTemps_deltas = [ days_ts(:) deltaTemps_beforeSort];
end
deltaDays  = DaysAndTemps_deltas(:,1);
deltaTemps = DaysAndTemps_deltas(:,2:end).';

deltaTemps_sep15start = circshift(deltaTemps,[0,108]);

% TO PLOT THE DELTAS:
%   figure; pcolor(deltaTemps_sep15start)
%   set(gca,'Ydir','reverse');
%   colorbar; caxis([-1.5 1.5]); colormap(cbrewer('div','RdBu',64))
%   sh=get(gca,'children');
%   sh.EdgeColor = 'none';

% % Apply the Deltas
tempMat = tempMat + deltaTemps;
putvar(deltaDays,deltaTemps,deltaTemps_sep15start)


%% Put everything in the Forcing Table's order of dys
% Match the Glider days up with the days (And Depths) in the Forcing File
frc_daysOfYr = unique(forcingTable.t);
nd = numel( frc_daysOfYr );

fulltempvec  = nan(nd, 1);
fulltimevec = nan(nd, 1);
fullkvec    = nan(nd, 1);
placecounter = 1;
for dii = 1:nd % loop over the forcing days (target times)
    % Get Days of Year from Forcing Table (even day numbers > 365)
    thisDoY = frc_daysOfYr(dii);
    thisDoY = mod(thisDoY,365);
    thisDoY( thisDoY == 0) = 365;
    
    fullkvec(placecounter:placecounter + nk-1) = klev;
    try
    fulltempvec(placecounter:placecounter + nk-1) = tempMat(:,thisDoY);
    catch
        error('a')
    end
    fulltimevec(placecounter:placecounter + nk-1) = thisDoY;
    
    placecounter = placecounter + nk;
end


forcingTable.fullkvec = fullkvec;
forcingTable.fulltempvec = roundsd(fulltempvec,4) + higherTempsByDeg;
forcingTable.fulltimevec = fulltimevec;

clear klev zmid frc_daysOfYr nk
clear fullkvec fulltempvec fulltimevec
clear nd thisDoY tempMat placecounter dii

%% Reorder and reformat the forcing table
forcingTable.t(:) = forcingTable.t(:)+1; % For leap year 2012, add one to the days past February
forcingTable.fulltimevec = [];
forcingTable.fullkvec = [];
forcingTable.w = zeros( size(forcingTable,1), 1);
forcingTable = forcingTable(:,[1:3 7 4:6]);
forcingTable.vdc = [];
forcingTable = forcingTable(:,[1:4 6 5]);
forcingTable.Properties.VariableNames{5} = 'temp';



