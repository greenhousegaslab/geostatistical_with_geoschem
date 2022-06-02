%-------------------------------------------------------------------------------------------------------------------------------%
% Script process_obs.m														%
% Purpose: Cut the OCO-2 observations into daily data files. The inverse model requires the observations in daily files.	%
% Written by Zichong Chen. Commented by Scot Miller, Aug. 10, 2021.								%
%																%
%-------------------------------------------------------------------------------------------------------------------------------%

%------------%
% NOTES:     %
%------------%

	% For some weird reason, Matlab will only create new netcdf files in the current working directory.
	% Hence, you need to start matlab from wherever you want to create new netcdf files.

	% Path to the observations: /scratch/groups/smill191/smiller/data/OCO2_MIP/ 

%---------------------------%
% Required function inputs  %
%---------------------------%

	filename='~/work/smiller/data/OCO2_MIP/OCO2_obs/OCO2_b10c_10sec_GOOD_r5.nc4';


%---------------------------------------------%
% Read in the OCO-2 data and data attributes  %
%---------------------------------------------%

	disp('Read in the OCO-2 data and data attributes');

	sid=ncread(filename,'sounding_id');
	dat=ncread(filename,'date');
	lon=ncread(filename,'longitude');
	lat=ncread(filename,'latitude');
	qff=ncread(filename,'xco2_quality_flag');
	profile_ap=ncread(filename,'co2_profile_apriori');			
	ak=ncread(filename,'xco2_averaging_kernel');
	% pl=ncread(filename,'pressure_levels'); make sure return the same name
	% pl here. from D Baker, using psurf *sigma levels
	% unc=ncread(filename,'xco2_uncertainty');
	apr=ncread(filename,'xco2_apriori'); 				
	co2=ncread(filename,'xco2');
	pw=ncread(filename,'pressure_weight');
	%load p_level.mat pl %calcuate beforehand and now load it here%
	altitude = ncread(filename,'altitude');
	type = ncread(filename,'data_type');
	assimflag = ncread(filename,'assimilate_flag');

	% GET PRESSURE LEVELS %
	ps=ncread(filename,'psurf');
	sl=ncread(filename,'sigma_levels');
	pls=sl*ps'; 
	clear sl ps; %hPa%

	% Uncertainties
	% Instead of using the measurement uncertainty, I'll use the measurement uncertainty plus the model uncertainty (in the data file)
	% In order to add two uncertainties, we need to convert them to variances and then add them
	% We can then take the square root to convert back to standard deviation
	unc1 = ncread(filename,'xco2_uncertainty');
	unc2 = ncread(filename,'model_error');
	unc  = sqrt(unc1.^2 + unc2.^2); 
	
	% Scot: I'm going to multiply the uncertainties by 1.5 -- more realistic given actual model-data uncertainties
	% unc = unc .* 1.5;


%-------------------------------------%
% Screen out unsuitable observations  %
%-------------------------------------%

	disp('Screen out unsuitable observations');

	% Remove the following observations:
	% Observations at extreme high or low latitudes (due to low signal-to-noise ratio)
	% Observations collected above land >3000m elevation (indicates mountainous terrain which GEOS-Chem is unlikely to adaquately capture)
	% Also, only keep land nadir and land glint
	% Only keep observations with a good quality flag (qff=0)
	% Only use observations over land (for the current simulations)
	% comment: Data type derived from surface_type, operation_mode, and land_fraction: 1=land nadir, 2=land glint, 3=land target, 4=land transition, 
	% 5=water nadir, 6=water glint, 7=water target, 8=water transition, 9=other (mixed land/water scenes)
	% We're also supposed to withold some data for independent evaluation. Data to be included are marked with 'assimilation_flag' = 1;

%	sel = lat>-60 & lat<60 & altitude < 3000 & type < 3 & qff == 0 & assimflag == 1;
	
	% LNLG case:
	% sel = altitude < 3000 & type < 3 & qff == 0 & assimflag == 1;
	% sel = altitude < 3000 & type < 3 & assimflag == 1 & lat < 60 & lat > -60;
	% sel = altitude < 3000 & type < 3 & assimflag == 1;

	% OG case:
	% sel = altitude < 3000 & type == 6 & qff == 0 & assimflag == 1;
	% sel = altitude < 3000 & type == 6 & assimflag == 1;

	% LNLGOG case
        % sel = altitude < 3000 & (type < 3 | type == 6) & qff == 0 & assimflag == 1;
	sel = altitude < 3000 & (type < 3 | type == 6) & assimflag == 1;

	% Case using all OCO2 observations (i.e., for evaluating model outputs)
	% sel = 1:length(assimflag);

	% Remove unwanted observations
	sid        = sid(sel);
	dat        = dat(:,sel);
	lon        = lon(sel);
	lat        = lat(sel);
	qff        = qff(sel);
	profile_ap = profile_ap(:,sel);
	ak         = ak(:,sel);
	unc        = unc(sel);
	apr        = apr(sel);
	co2        = co2(sel);
	pw         = pw(:,sel);
	pls        = pls(:,sel);
	assimflag  = assimflag(sel);

	disp('Total number of observations that fit criteria');
	disp(num2str(length(sid)));


%----------------------------------------------------%
% Loop over each year and process the observations   %
%----------------------------------------------------%

	disp('Loop over each year and process the observations');

for year = 2014:2021,
    disp(num2str(year));
    for month = 1:12,
        for day = 1:31,

	% Find all OCO-2 observations that are from the specific day of interest.
	sel = find(dat(1,:)==year & dat(2,:)==month & dat(3,:)==day);


%-------------------------------------------------------------------------%
% Collect the observations (if there are more than 5 obs on a given day)  %
%-------------------------------------------------------------------------%
 
	if length(sel)>5, %make sure we have at least more than five data points for the specific day%


%--------------------------------------------%
% Initialize the netcdf file                 %
%--------------------------------------------%

	if month<10; month1=strcat('0',num2str(month)); else; month1=num2str(month); end;
        if day<10; day1=strcat('0',num2str(day)); else; day1=num2str(day); end;

        % Open the netCDF file.
        ncid = netcdf.create(strcat('oco2_LtCO2_',num2str(year),month1,day1,'.nc'),'CLOBBER');
        % Define the DIMENSIONS of the variable.
        dimid = netcdf.defDim(ncid,'sounding_id',length(sel));
        dimid1 = netcdf.defDim(ncid,'levels',20);
        dimid2 = netcdf.defDim(ncid,'epoch_dimension',7);
                
%-------------------------------------%
% Populate the netcdf file with data  %
%-------------------------------------%

	% Define a new VARIALBE in the file.
        my_varID = netcdf.defVar(ncid,'co2_profile_apriori','double',[dimid1 dimid]);
        my_varID1 = netcdf.defVar(ncid,'date','double',[dimid2 dimid]);
        my_varID2 = netcdf.defVar(ncid,'sounding_id','double',dimid);
        my_varID3 = netcdf.defVar(ncid,'xco2_averaging_kernel','double',[dimid1 dimid]);
        my_varID4 = netcdf.defVar(ncid,'pressure_levels','double',[dimid1 dimid]);
        my_varID5 = netcdf.defVar(ncid,'latitude','double',dimid);
        my_varID6 = netcdf.defVar(ncid,'longitude','double',dimid);
        my_varID7 = netcdf.defVar(ncid,'xco2_uncertainty','double',dimid);
        my_varID8 = netcdf.defVar(ncid,'xco2_apriori','double',dimid);
        my_varID9 = netcdf.defVar(ncid,'xco2','double',dimid);
        my_varID10 = netcdf.defVar(ncid,'pressure_weight','double',[dimid1 dimid]);
	my_varID11 = netcdf.defVar(ncid,'assimilate_flag','double',dimid);
                
        netcdf.endDef(ncid);
        % WRITE DATA to variable.
        netcdf.putVar(ncid,my_varID,profile_ap(:,sel));
        netcdf.putVar(ncid,my_varID1,dat(:,sel));
        netcdf.putVar(ncid,my_varID2,sid(sel));
        netcdf.putVar(ncid,my_varID3,ak(:,sel));
        netcdf.putVar(ncid,my_varID4,pls(:,sel));
        netcdf.putVar(ncid,my_varID5,lat(sel));
        netcdf.putVar(ncid,my_varID6,lon(sel));
        netcdf.putVar(ncid,my_varID7,unc(sel));
        netcdf.putVar(ncid,my_varID8,apr(sel));
        netcdf.putVar(ncid,my_varID9,co2(sel));
        netcdf.putVar(ncid,my_varID10,pw(:,sel));
	netcdf.putVar(ncid,my_varID11,assimflag(sel));

        % Verify that the variable was created.
        netcdf.close(ncid)

	% End if statements and loops                
	end; % End of length(sel) if statement
	end; % End of day loop
	end; % End of month loop
	end; % End of year loop               

%-------------------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT
%-------------------------------------------------------------------------------------------------------------------------------%


