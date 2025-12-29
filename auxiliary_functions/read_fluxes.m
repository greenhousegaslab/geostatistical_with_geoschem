%---------------------------------------------------------------------------------------%
% FUNCTION: read_fluxes.m								%
% PURPOSE:  Read in different CO2 flux components for the inverse model			%
% K. Morgan, 11/04/2024									%
%											%
%---------------------------------------------------------------------------------------%


%-----------------%
% NOTES:
% time dimension? %
%-----------------%


%-------------------%
% BEGIN FUNCTION    %
%-------------------%

function [ X, ocean ] = read_fluxes(yr);

	disp('Running the read_fluxes function');

	if yr == 2014; 
		start_day = datetime('20140901', 'InputFormat', 'yyyyMMdd');
	else;
                start_day = datetime(strcat(num2str(yr),'0101'), 'InputFormat', 'yyyyMMdd');
	end;
	
	if yr == 2024;
		end_day   = datetime(strcat('20240731'), 'InputFormat', 'yyyyMMdd');
	else;
                end_day   = datetime(strcat(num2str(yr),'1231'), 'InputFormat', 'yyyyMMdd');
	end; 

	nday      = days(end_day - start_day + 1);


%----------------------------%
% Loop over individual days  %
%----------------------------%

	X = [];

	for d = 1:nday;

	yyyyMMdd = datestr(start_day + d - 1, 'yyyymmdd');
	MMdd     = datestr(start_day + d - 1, 'mmdd');

%-----------------------------%
% READ IN BIOSPHERIC FLUXES   %
%-----------------------------%
        
        % CASA available through 2020
	if yr <= 2019;
        	NEE = ncread(strcat('/scratch4/smill191/OCO2_MIP_inputs/CASA_GFED4/regridded/avg_daily_CT2022_CASA_',yyyyMMdd,'.nc'), 'b4');
	else;
        	NEE = ncread(strcat('/scratch4/smill191/OCO2_MIP_inputs/CASA_GFED4/regridded/avg_daily_CT2022_CASA_2019',MMdd,'.nc'), 'b4');
	end;

%------------------------------%
% READ IN FOSSIL FUEL FLUXES   %
%------------------------------%

        % Read in ODIAC fluxes
        % ODIAC available through 2024
	% Leap years: 2016, 2020, 2024
        odiac = ncread(strcat('/scratch4/smill191/OCO2_MIP_inputs/ODIAC/regridded/avg_daily_fossil_fuel_4x5_',yyyyMMdd,'.nc'), 'emission');  % units micromol m-2 s-1

%----------------------------------%
% READ IN BIOMASS BURNING FLUXES   %
%----------------------------------%

        % Read in GFED fluxes
        % GFED available through 2020
        % Leap years: 2016, 2020
	if yr <= 2020;
		gfed = ncread(strcat('/scratch4/smill191/OCO2_MIP_inputs/GFED5/regridded/avg_daily_GFED5_',yyyyMMdd,'.nc'), 'CO2');
	else;
                gfed = ncread(strcat('/scratch4/smill191/OCO2_MIP_inputs/GFED5/regridded/avg_daily_GFED5_2020',MMdd,'.nc'), 'CO2');
	end;
	
%-------------------------------%
% ADD TOGETHER ALL THE FLUXES   %
%-------------------------------%

	totflux = odiac + NEE + gfed;  % no ocean

        % convert the daily emission matrix to a vector
        emis_day   = reshape(totflux, [], 1);

        % attach the next day to a vector
        X = [X; emis_day];

        end; % End loop over different days


%------------------------%
% READ IN OCEAN FLUXES   %
%------------------------%

        % Some flux products aren't available after 2018. Reset the flux year to 2018 where needed
        yr2 = yr;
        if yr > 2017; yr2 = 2017; end;
	if yr < 2015; yr2 = 2015; end;

        % Read in ocean fluxes 
        load(strcat('/scratch4/smill191/OCO2_MIP_inputs/takahashi/ocean_',num2str(yr2),'.mat'));
        ocean = double(ocean);

        % Deal with leap year for year 2020 (repeat Feb fluxes for an additional day)
        if (yr==2020 | yr==2024);

        % Feb 28 is the 59th day of the year
        % I want to repeat the 59th day in NEE and ocean (other flux types account for leap year already)
	% For year 2024, cut the ocean fluxes at July 31st (Aug - Dec 2024 not included in v11 MIP)

        m1   = length(ocean)./365;
        sel  = ((59-1).*m1 + 1):(59.*m1);
        sela = 1:(59.*m1);
	if yr==2020; selb = (59.*m1 +1):length(ocean); end;
	if yr==2024; selb = (59.*m1 +1):(212.*m1); end;
        ocean = [ ocean(sela); ocean(sel); ocean(selb) ];

        end; % End of yr if statement

	% For 2014, just take Sept 1 - Dec 31
	if yr==2014;
        	m1   = length(ocean)./365;
		sel = (1 + m1.*(31 + 28 + 31 + 30 + 31 + 30 + 31 + 31)):length(ocean);
		ocean = ocean(sel);
	end;

	disp('End of read_fluxes function');

%---------------------------------------------------------------------------------------%
% END OF SCRIPT
%---------------------------------------------------------------------------------------%
