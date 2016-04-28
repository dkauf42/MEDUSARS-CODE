function [newmmtable,forcingTable,gliderMLDs] = vdcTimeSeriesFromGliderMLDwROMSdeltas(varargin)


%% Initial Function Housekeeping 
    p=inputParser;
    addParameter(p,'makefig',false,@islogical);
parse(p,varargin{:});
makefig = p.Results.makefig;

%% Get Data Files and Perform Other Housekeeping Operations
sourcedir = '~/Desktop/CATEGORIES/CAREER_MANAGEMENT/VIMS/Dissertation/Chapter_2_ForwardModel/Glider_Data_ALL/RMJ_Glider_2012-2013_Data_Structures/';
load([sourcedir,'Glider_SG503_2012.mat']); % get 'SG503_2012' variable
gliderTable = SG503_2012;  clear SG503_2012
forcingTable=marmot_tabread([sourcedir,'fkt_base2010_roms_vdcFromMLDandLlort_rsA.dat']);

gliderTable = gliderTable(gliderTable.used(:) ~= 0,:); % remove 'unused' data rows

% Get ROMS Output Climate Deltas
ncsourcedir = '~/Desktop/CATEGORIES/CAREER_MANAGEMENT/VIMS/Dissertation/Chapter_2_ForwardModel/climate_scenarios_roms_outputs/';
ncA = [ncsourcedir, 'kaufman.data.future.delta.A.nc'];  iA = ncinfo(ncA);
ncB = [ncsourcedir, 'kaufman.data.future.delta.B.nc'];  iB = ncinfo(ncB);
ncSpatialAvg1 = [ncsourcedir, 'spatially_averaged/', 'kaufman.data.future.delta.area1.nc'];  iSA1 = ncinfo(ncSpatialAvg1);
ncSpatialAvg2 = [ncsourcedir, 'spatially_averaged/', 'kaufman.data.future.delta.area2.nc'];  iSA2 = ncinfo(ncSpatialAvg2);
ncSpatialAvg3 = [ncsourcedir, 'spatially_averaged/', 'kaufman.data.future.delta.area3.nc'];  iSA3 = ncinfo(ncSpatialAvg3);


%%%%%%%%%% COMMON TWEAKS %%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Apply deltas from the ROMS climate scenario outputs %%%%%

deltasToApply = '2050';
% deltasToApply = '2100';
% deltasToApply = 'none';

% deltaSite = 'A';  % A: 77.2S  169.5E
% % deltaSite = 'B';  % B: 76.5S  176E
% deltaSite = 'spatialAvg1';  % 1:  167 - 173 E  ;  77.5 - 76.6 S
deltaSite = 'spatialAvg2';  % 2:  166 - 174 E  ;  77.5 - 76 S
% deltaSite = 'spatialAvg3';  % 3:  164 - 176 E  ;  78 - 75 S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Make MLDs from the glider deeper by a percentage: %%%%%
% deeperMLDbyPrctMultiplier = 0.12; % 12% reduction
% deeperMLDbyPrctMultiplier = 0.44; % 44% reduction
% deeperMLDbyPrctMultiplier = 1; % 100% reduction
deeperMLDbyPrctMultiplier = 0; % No reduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Make MLDs from the glider start on an earlier date %%%%%%%%%%%%%%%
%%%%%        (by stretching the time-series):            %%%%%%%%%%%%%%%
stretchEarlierByDays = 0; % don't stretch the time-series any earlier
% stretchEarlierByDays = ceil(8.5); % stretch by 8.5 days (rounded up to 9)
% stretchEarlierByDays = ceil(19.2); % stretch by 19.2 days (rounded up to 20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% What to do w/period between IceMeltDay and 1st Glider MLD Datum %%%%
% winterToGliderMLDmethod = 'interpBetween';
winterToGliderMLDmethod = 'copyGliderMLDs';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Ice Melt Day %%%%%
IceMeltDay = datenum('13-Oct-2012');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffusion Values to use (comment out the Ksurface variable to make the
%           Kmixlayer cover everything up to the surface

Ksurface = 10^-3.5 * 86400; % 10-4 (86400 converts m2/s converted to m2/d)
Kmixlayer = 10^-0.5 * 86400;
Kbottom = 10^-3 * 86400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get MLDs for each day from the ROMS scenario outputs

switch deltaSite
    case 'A'
        ncROMSclimateFile = ncA;
    case 'B'
        ncROMSclimateFile = ncB;
    case 'spatialAvg1'
        ncROMSclimateFile = ncSpatialAvg1;
    case 'spatialAvg2'
        ncROMSclimateFile = ncSpatialAvg2;
    case 'spatialAvg3'
        ncROMSclimateFile = ncSpatialAvg3;
end

% Time    % days since Sep. 15, 2012, in 5 day increments
nctime = double( ncread(ncROMSclimateFile,'time') );
%  MLD: depth where sigmaT is 0.01 different from the surface value
% 73 x 1 single
mld_2050 = double( ncread(ncROMSclimateFile,'mld_2050') );
mld_2100 = double( ncread(ncROMSclimateFile,'mld_2100') );

% Convert ROMS days (from Sep.15,2012) to days-of-year
ncdoy = date2doy( datenum('Sep 15 2012') + nctime );

% Determine where the year splits
firsttimestepofnewyear = 23;
lasttimestepoffirstyear = 22;
             
mld_dV_2050 = [ ones(1,5) .* mld_2050(1),...
                 interp1(ncdoy(1:lasttimestepoffirstyear), mld_2050(1:lasttimestepoffirstyear), 262:365) ,...
                 interp1([-4.5; 0.5; ncdoy(firsttimestepofnewyear:end)], mld_2050(21:end), 1:255),...
                 mld_2050(end) ];
mld_dV_2100 = [ ones(1,5) .* mld_2100(1),...
                 interp1(ncdoy(1:lasttimestepoffirstyear), mld_2100(1:lasttimestepoffirstyear), 262:365) ,...
                 interp1([-4.5; 0.5; ncdoy(firsttimestepofnewyear:end)], mld_2100(21:end), 1:255),...
                 mld_2100(end) ];

days_ts = 257:621;
days_ts_dates = doy2date(days_ts,repmat(2012,numel(days_ts),1));
days_ts( days_ts>365 ) = days_ts( days_ts>365 )-365;

switch deltasToApply
    case '2050'
                deltaMLDs_beforeSort = mld_dV_2050(:);
                DaysAndMlds_deltas = sortrows( [ days_ts(:) deltaMLDs_beforeSort] );
    case '2100'
                deltaMLDs_beforeSort = mld_dV_2100(:);
                DaysAndMlds_deltas = sortrows( [ days_ts(:) deltaMLDs_beforeSort] );
    case 'none'
        deltaMLDs_beforeSort = zeros(365,1);
        DaysAndMlds_deltas = [ days_ts(:) deltaMLDs_beforeSort];
end
deltaDays = DaysAndMlds_deltas(:,1);
deltaMLDs = DaysAndMlds_deltas(:,2);















% % % % % Time    % days since Sep. 15, 2012, in 5 day increments
% % % % timeA = double( ncread(ncA,'time') );
% % % % timeB = double( ncread(ncB,'time') );
% % % % %  MLD: depth where sigmaT is 0.01 different from the surface value
% % % % % 73 x 1 single
% % % % mld_A2050 = double( ncread(ncA,'mld_2050') );
% % % % mld_A2100 = double( ncread(ncA,'mld_2100') );
% % % % mld_B2050 = double( ncread(ncB,'mld_2050') );
% % % % mld_B2100 = double( ncread(ncB,'mld_2100') );
% % % % 
% % % % % Convert ROMS days (from Sep.15,2012) to days-of-year
% % % % doyA = date2doy( datenum('Sep 15 2012') + timeA );
% % % % doyB = date2doy( datenum('Sep 15 2012') + timeB );
% % % % 
% % % % % Determine where the year splits
% % % % firsttimestepofnewyear = 23;
% % % % lasttimestepoffirstyear = 22;
% % % %              
% % % % mld_dV_A2050 = [ ones(1,5) .* mld_A2050(1),...
% % % %                  interp1(doyA(1:lasttimestepoffirstyear), mld_A2050(1:lasttimestepoffirstyear), 262:365) ,...
% % % %                  interp1([-4.5; 0.5; doyA(firsttimestepofnewyear:end)], mld_A2050(21:end), 1:255),...
% % % %                  mld_A2050(end) ];
% % % % mld_dV_A2100 = [ ones(1,5) .* mld_A2100(1),...
% % % %                  interp1(doyA(1:lasttimestepoffirstyear), mld_A2100(1:lasttimestepoffirstyear), 262:365) ,...
% % % %                  interp1([-4.5; 0.5; doyA(firsttimestepofnewyear:end)], mld_A2100(21:end), 1:255),...
% % % %                  mld_A2050(end) ];
% % % % mld_dV_B2050 = [ ones(1,5) .* mld_B2050(1),...
% % % %                  interp1(doyB(1:lasttimestepoffirstyear), mld_B2050(1:lasttimestepoffirstyear), 262:365) ,...
% % % %                  interp1([-4.5; 0.5; doyB(firsttimestepofnewyear:end)], mld_B2050(21:end), 1:255),...
% % % %                  mld_A2050(end) ];
% % % % mld_dV_B2100 = [ ones(1,5) .* mld_B2100(1),...
% % % %                  interp1(doyB(1:lasttimestepoffirstyear), mld_B2100(1:lasttimestepoffirstyear), 262:365) ,...
% % % %                  interp1([-4.5; 0.5; doyB(firsttimestepofnewyear:end)], mld_B2100(21:end), 1:255),...
% % % %                  mld_A2050(end) ];
% % % % 
% % % % days_ts = 257:621;
% % % % days_ts_dates = doy2date(days_ts,repmat(2012,numel(days_ts),1));
% % % % days_ts( days_ts>365 ) = days_ts( days_ts>365 )-365;
% % % % 
% % % % switch deltasToApply
% % % %     case '2050'
% % % %         switch deltaSite
% % % %             case 'A'
% % % %                 deltaMLDs_beforeSort = mld_dV_A2050(:);
% % % %                 DaysAndMlds_deltas = sortrows( [ days_ts(:) deltaMLDs_beforeSort] );
% % % %             case 'B'
% % % %                 deltaMLDs_beforeSort = mld_dV_B2050(:);
% % % %                 DaysAndMlds_deltas = sortrows( [ days_ts(:) deltaMLDs_beforeSort] );
% % % %         end
% % % %     case '2100'
% % % %         switch deltaSite
% % % %             case 'A'
% % % %                 deltaMLDs_beforeSort = mld_dV_A2100(:);
% % % %                 DaysAndMlds_deltas = sortrows( [ days_ts(:) deltaMLDs_beforeSort] );
% % % %             case 'B'
% % % %                 deltaMLDs_beforeSort = mld_dV_B2100(:);
% % % %                 DaysAndMlds_deltas = sortrows( [ days_ts(:) deltaMLDs_beforeSort] );
% % % %         end
% % % %     case 'none'
% % % %         deltaMLDs_beforeSort = zeros(365,1);
% % % %         DaysAndMlds_deltas = [ days_ts(:) deltaMLDs_beforeSort];
% % % % end
% % % % deltaDays = DaysAndMlds_deltas(:,1);
% % % % deltaMLDs = DaysAndMlds_deltas(:,2);


%% Get MLDs for each day in the Glider Data Table
t_sg = gliderTable.m_date;
mld_sg = gliderTable.mld;

% Get Indices for each day in glider data table
daysVector = floor( t_sg(1) ) : 1 : ceil( t_sg(end) );
[~,binEdge,binIdx]=histcounts(t_sg, daysVector);
binMid = (binEdge(1:end-1) + binEdge(2:end) ) ./ 2;

% Get Average MLD for each day in glider data table
DailyAvgMLD = accumarray(binIdx,mld_sg,[],@nanmean);

%%%%%% DECREASE MIXED LAYER DEPTHS BY A PERCENTAGE %%%%%%
DailyAvgMLD = DailyAvgMLD - (DailyAvgMLD .* deeperMLDbyPrctMultiplier);
DailyAvgMLD(DailyAvgMLD<2.5) = 2.6;

% Make a variable to hold the MLDs
DaysAndMlds = [floor(binMid(:)) DailyAvgMLD]; %%%### @floor changes datenums from noon to midnight of each day

%%%%%% STRETCH MIXED LAYER DEPTHS EARLIER BY NUMBER OF DAYS %%%%%%
stretchedDays = [ ((DaysAndMlds(1,1)-stretchEarlierByDays) : DaysAndMlds(1,1)).'; DaysAndMlds(2:end,1) ];
origDaysHigherResSpace = linspace( DaysAndMlds(1,1), DaysAndMlds(end,1), numel(stretchedDays) );
stretchedMLDs = interp1( DaysAndMlds(:,1), DaysAndMlds(:,2), origDaysHigherResSpace );
DaysAndMlds = [stretchedDays stretchedMLDs(:)];

% Convert glider dates to DaysOfYear
DaysAndMlds(:,3) = date2doy( DaysAndMlds(:,1) );

% remove leap year numbers (glider was from 2012, a leap year; so, dec 31
%   is day #366.  But we just want a year of days 1 to 365
lastDayOfYear = find(DaysAndMlds(:,3)==366);
DaysAndMlds(1:lastDayOfYear,3) = DaysAndMlds(1:lastDayOfYear,3) - 1;
% DaysAndMlds(1,:)=[];

clear AvgDuration binEdge binIdx binMid t_sg mld_sg daysVector DailyAvgMLD

% % % % %% Interpolate Glider MLDs for rest of year (between day 39 and day 327)
% % % % winterDays = 40:325;
% % % % winterMLDs = interp1([39 326], [DaysAndMlds(end,2), DaysAndMlds(1,2)], 40:325);
% % % % nw = numel(winterDays);

%% Copy the first glider MLDs to the rest of the year

wd = DaysAndMlds(end,3)+1:DaysAndMlds(1,3)-1; % winterDays (from last day of glider, back around to the first day of glider data)
winterMLDs = [];
% winterMLDs = interp1([39 326], [DaysAndMlds(end,2), DaysAndMlds(end,2)], 40 : 325);
% winterMLDs = interp1([39 326], [200, 200], 40 : 325);

    % Make MLD 200m deep during winter, and then linearly interpolate to
    %   first glider datum starting on 'iceMeltDay' or just copy the first
    %   glider day MLD back to the iceMeltDay.
                            % imd = 315; % 315 is day of year for Nov. 10 (during a leap year)
imd = date2doy(IceMeltDay); % 320 is day of year for Nov. 15 (during a leap year); 313 is day of year for Nov. 08 (during a leap year)
interpDays = wd(end)-imd;
%%%%%%%%%%%%%%% Put 200m MLDs in for all of the 'winter days' %%%%%%%%%%
winterMLDs(1:imd-wd(1)) = interp1([wd(1)-1 imd], [200, 200], wd(1) : imd-1);

%%%%%%%%%%%%%%% What to do with intermediate period between %%%%%%%%%%
%%%%%%%%%%%%%%%      winter 200m MLDs and glider MLDs       %%%%%%%%%%
switch winterToGliderMLDmethod
    case 'interpBetween'
        winterMLDs(end+1 : end+interpDays+1) = interp1([imd-1 wd(end)+1],...
            [200, DaysAndMlds(1,2)], imd : wd(end));  % linearly interpolate to first glider datum starting from 'iceMeltDay'
    case 'copyGliderMLDs'
        winterMLDs(end+1 : end+interpDays+1) = interp1([imd-1 wd(end)+1],...
            [ DaysAndMlds(1,2), DaysAndMlds(1,2)], imd : wd(end)); % copy the first glider day MLD back to the iceMeltDay.
end
    
nw = numel(wd);

% Make a complete year time-series (and sort it in order of Days of Year)
DaysAndMlds = [DaysAndMlds ; [ NaN(nw,1) , winterMLDs(:) , wd(:) ] ];
DaysAndMlds = sortrows( DaysAndMlds, 3);
gliderMLDs = array2table(DaysAndMlds,'VariableNames',{'matlabdates','mld','DoY'});

clear nw wd winterMLDs imd interpDays

%% Apply the climate scenario deltas from the ROMS outputs, if specified


% DaysAndMlds_A2050 = [ days_ts(:) mld_dV_A2050(:)];

% shallower future mixed layers should give positive numbers (remember that
% mixed layer depths are defined in the negative direction) ** FOR ROMS **

% For example, a delta of -0.7 means that the future case is 0.7m deeper,
% so we need to add +0.7 to the glider MLDs because they are in the
% positive direction
gliderMLDs.mld = gliderMLDs.mld - deltaMLDs;

% reset the MLDs that have become deeper than 200m with the climate deltas, back to 200m
gliderMLDs.mld( gliderMLDs.mld > 200 ) = 200;
% reset the MLDs that have become above the surface, back to 0m
gliderMLDs.mld( gliderMLDs.mld < 0 ) = 0;

%% Go through each day in the Forcing File..
%   and calculate the profile for each day based on the glider MLD for that day

% Match the Glider days up with the days (And Depths) in the Forcing File
frc_daysOfYr = unique(forcingTable.t);
nd = numel( frc_daysOfYr );

% K to Z Transformation
klev = (1:41).'; nk = numel(klev);
ztop = (0:5:200).';
% % KtoZtable=readtable('zlevel_rsA.dat','delimiter','\t');
% % nk = numel(KtoZtable.k);

fullvdcvec  = nan(nd, 1);
fulltimevec = nan(nd, 1);
fullkvec    = nan(nd, 1);
vdc = nan(nk,1);
placecounter = 1;
for dii = 1:nd % loop over the forcing days (target times)
    
    % Get Days of Year from Forcing Table (even day numbers > 365)
    thisDoY = frc_daysOfYr(dii);
    thisDoY = mod(thisDoY,365);
    thisDoY( thisDoY == 0) = 365;
    
    % Find this Day in the glider Table, and get the mixed layer depth for it
    mx = abs( gliderMLDs.mld( gliderMLDs.DoY == thisDoY ) );

    % Find the k-level & z that corresponds closest to the mixed layer depth
    [~,klevi] = min( abs( ztop - mx ) );
    %mxZ = ztop(klevi);

    if exist('Ksurface','var')
        vdc(1) = Ksurface;
    else
        vdc(1) = Kmixlayer;
    end
    vdc(2:klevi) = Kmixlayer;
    vdc(klevi+1:end) = Kbottom;
    
    % Smooth Diffusivities for levels just above and below mixed layer
    if klevi > 2
        if klevi<nk-1
            kx = [klevi-2 klevi+2];
            ky = [vdc(klevi-2) vdc(klevi+2)];
            kq = klevi-1:klevi+1;
            vdc(klevi-1:klevi+1) = interp1(kx,ky,kq);
        else
            kx = [klevi-2 klevi];
            ky = [vdc(klevi-2) vdc(klevi)];
            kq = klevi-1:klevi;
            vdc(klevi-1:klevi) = interp1(kx,ky,kq);
        end
    elseif klevi == 1
        kx = [1 klevi+2];
        ky = [vdc(1) vdc(klevi+2)];
        kq = klevi:klevi+1;
        vdc(klevi:klevi+1) = interp1(kx,ky,kq);
    else
        kx = [1 klevi+2];
        ky = [vdc(1) vdc(klevi+2)];
        kq = klevi-1:klevi+1;
        vdc(klevi-1:klevi+1) = interp1(kx,ky,kq);
    end
    
    fullkvec(placecounter:placecounter + nk-1) = klev;
    fullvdcvec(placecounter:placecounter + nk-1) = vdc;
    fulltimevec(placecounter:placecounter + nk-1) = thisDoY;
    
    placecounter = placecounter + nk;
end

forcingTable.fullkvec = fullkvec;
forcingTable.fullvdcvec = fullvdcvec;
forcingTable.fulltimevec = fulltimevec;

newmmtable = forcingTable;
newmmtable(:,'fulltimevec') = [];
newmmtable(:,'fullkvec') = [];
newmmtable(:,'vdc') = [];
newmmtable.Properties.VariableNames{5} = 'vdc';
newmmtable = newmmtable(:,[1:3 5 4]);

clear nd klev ztop frc_daysOfYr nk kx ky kq
clear fullkvec fullvdcvec fulltimevec
clear thisDoY mx klevi mxZ vdc placecounter dii
clear Ksurface Kmixlayer Kbottom


%% Make Figures
if makefig
    dates = doy2date(newmmtable.t,repmat(2012,length(newmmtable.t),1));
    min(dates)
    XL = [min(dates) 735266]; % XL = [258 396]; % until feb-1-2013
    xlab = '';
    
    FH = figure();
    depths = repmat((0:5:200).',length(newmmtable.k)./41,1);
    tkzPlot(dates,depths,newmmtable.vdc); % Plot Diffusion time-depth profiles
    hold all;
    lh=line(gliderMLDs.matlabdates,gliderMLDs.mld,'color','w','linewidth',2); % Plot MLD time-series
    lh.ZData=repmat(99999999,1,365); % increase the Z value, so that the line is on 'top' of the surface plot
    uistack(lh,'top')
    
    xlim(XL); xlabel(xlab); grid on;
    ylabel('depth (m)')
    datetick('x','mm/dd/yyyy','keeplimits')
    
    % % FH = figure(); h.ax1 = gca;
    % % plot(gliderMLDs.matlabdates,gliderMLDs.mld); % Plot MLD time-series
    % % ylim([0 200]);
    % %
    % % h.ax1.YDir = 'reverse';
    % % xlim(XL); xlabel(xlab); grid on;
    % % ylabel('depth (m)')
    % % datetick('x','mm/dd/yyyy','keeplimits')
    
    aaa=gliderMLDs.matlabdates;
    putvar(aaa)
    figure();
    lh=line(days_ts_dates(1:149),deltaMLDs_beforeSort(1:149),'color','k','linewidth',2); % Plot MLD time-series
    title('Delta MLDs');
    putvar(days_ts,days_ts_dates,DaysAndMlds,deltaDays,deltaMLDs,lh);
    xlim(XL); xlabel(xlab); grid on;
%     xlim([257,432])
    ylabel('depth (m)')
    datetick('x','mm/dd/yyyy','keeplimits')
end




putvar(newmmtable,forcingTable)

end


