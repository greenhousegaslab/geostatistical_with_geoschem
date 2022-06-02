%-------------------------------------------------------------------------------------------------------------------------------%
% SCRIPT: process_tccon.m													%
% PURPOSE: Reformat the TCCON data files to look the same as the OCO-2 data files.						%
% 	Then I can run the TCCON observations directly through the GOSAT observation operator in the GEOS-Chem adjoint.		%
% S. Miller, May 7, 2022													%
%																%
%-------------------------------------------------------------------------------------------------------------------------------%

%------------%
% NOTES:     %
%------------%

	% Scot: script is debugged and works. I need to actually run this script to generate the daily observations.

	% For some weird reason, Matlab will only create new netcdf files in the current working directory.
	% Hence, you need to start matlab from wherever you want to create new netcdf files.

	% Path to the observations: /scratch/groups/smill191/smiller/data/OCO2_MIP/ 

%---------------------------%
% Required function inputs  %
%---------------------------%

	tccondir = '/scratch/groups/smill191/smiller/data/OCO2_MIP/tccon_obs/';

%----------------------------------------------------%
% Loop over each year and process the observations   %
%----------------------------------------------------%

	disp('Loop over each year and process the observations');

	counter = 1;
	tcconobnall = [];

	for year = 2015:2020,
    	disp(num2str(year));
    	for month = 1:12,

        % Loop over days of the month
        ndays = eomday(year,month);
        for day=1:ndays;

        monlead = '';
        if(length(num2str(month))==1); monlead = '0'; end;
        dayslead = '';
        if(length(num2str(day))==1); dayslead = '0'; end;


%------------------------------------%
% Read in data from the TCCON file   %
%------------------------------------%

        % Read in TCCOn data descriptors
        obsfile    = strcat(tccondir,'tccon_timeavg_',num2str(year),monlead,num2str(month),dayslead,num2str(day),'.nc4');
	tccondate  = ncread(obsfile,'CO2/cdate');
        tcconid    = ncread(obsfile,'CO2/station_id');
        tcconobn   = ncread(obsfile,'CO2/obs_num');
	tcconobnall= [tcconobnall; tcconobn];
        tcconlat   = ncread(obsfile,'CO2/latitude');
        tcconlon   = ncread(obsfile,'CO2/longitude');

        % Read in TCCON data
        tcconco2   = ncread(obsfile,'CO2/column_mixing'); % TCCON dry air CO2 mole fraction
        tcconak    = ncread(obsfile,'CO2/avg_kernel'); % TCCON averaging kernel
        tcconp     = ncread(obsfile,'CO2/p_levels_ak'); % Pressure boundaries for avg kernel 
        tcconprior = ncread(obsfile,'CO2/prior_mixing'); % Dry air prior for TCCON
	psurf      = ncread(obsfile,'CO2/p_surf'); % Surface pressure

        % Convert TCCON pressure from pascals to hPA
        tcconp = tcconp ./ 100;
	psurf  = psurf ./ 100;

        % Calculate the TCCON pressure weighting function
        % Number of pressure levels is the same as the number of CO2 column levels --
        % means I must have pressure levels for grid midpoints (not edges)
        % OR maybe the last pressure edge is at 0.
	pwf    = zeros(71,length(tcconid));	
	xprior = zeros(length(tcconid),1);
	for(k = 1:length(tcconid));
	pedges      = [psurf(k); (tcconp(1:70,k) + tcconp(2:71,k))./2; 0];
	pwf(:,k)    = (pedges(1:71) - pedges(2:72)) ./ psurf(k);

	% disp('----------');
	% disp('pressure edges');
	% disp(num2str(pedges));
	% disp('pressure weighting');
	% disp(num2str(pwf(:,k)));
	% disp('pwf sum');
	% disp(num2str(sum(pwf)));
	% pause;

        % pwf(:,k)    = [ (tcconp(1:70,k) - tcconp(2:71,k)); tcconp(71,k)' ] ./ (tcconp(1,k)-tcconp(71,k));
	xprior(k) = sum(pwf(:,k) .* tcconprior(:,k)); 
	end; % End of k loop


%-----------------------------------%
% Reformat the TCCON data objects   %
%-----------------------------------%

	% Flip the presssure levels to be top to bottomo of the atmosphere
	tcconak    = flip(tcconak,1);
	tcconp     = flip(tcconp,1); 
	tcconprior = flip(tcconprior,1);
	pwf        = flip(pwf,1);

	% Reformat the date vector
	% I think tccondate is the right format but just needs an additional column	
	% Additional column is for millisecond
	tccondate = [tccondate; zeros(1,length(tcconid))];


%--------------------------------------------%
% Initialize the netcdf file                 %
%--------------------------------------------%

	if month<10; month1=strcat('0',num2str(month)); else; month1=num2str(month); end;
        if day<10; day1=strcat('0',num2str(day)); else; day1=num2str(day); end;

        % Open the netCDF file.
        ncid = netcdf.create(strcat('tccon_',num2str(year),month1,day1,'.nc'),'CLOBBER');

        % Define the DIMENSIONS of the variable.
        dimid = netcdf.defDim(ncid,'sounding_id',length(tcconid));
        dimid1 = netcdf.defDim(ncid,'levels',71);
        dimid2 = netcdf.defDim(ncid,'epoch_dimension',7);
        
        
%-------------------------------------%
% Populate the netcdf file with data  %
%-------------------------------------%

	% Define a new VARIALBE in the file.
        my_varID = netcdf.defVar(ncid,'co2_profile_apriori','double',[dimid1 dimid]);
        my_varID1 = netcdf.defVar(ncid,'date','double',[dimid2 dimid]);
        my_varID2 = netcdf.defVar(ncid,'sounding_id','NC_INT',dimid);
        my_varID3 = netcdf.defVar(ncid,'xco2_averaging_kernel','double',[dimid1 dimid]);
        my_varID4 = netcdf.defVar(ncid,'pressure_levels','double',[dimid1 dimid]);
        my_varID5 = netcdf.defVar(ncid,'latitude','double',dimid);
        my_varID6 = netcdf.defVar(ncid,'longitude','double',dimid);
        my_varID7 = netcdf.defVar(ncid,'xco2_uncertainty','double',dimid);
        my_varID9 = netcdf.defVar(ncid,'xco2','double',dimid);
        my_varID10 = netcdf.defVar(ncid,'pressure_weight','double',[dimid1 dimid]);
	my_varID11 = netcdf.defVar(ncid,'station_id','NC_INT',dimid);               
	my_varID12 = netcdf.defVar(ncid,'xco2_apriori','double',dimid);

        netcdf.endDef(ncid);
        % WRITE DATA to variable.
        netcdf.putVar(ncid,my_varID,tcconprior);
        netcdf.putVar(ncid,my_varID1,tccondate);
        netcdf.putVar(ncid,my_varID2,counter:(counter+length(tcconobn)-1));
        counter = counter+length(tcconobn);
        netcdf.putVar(ncid,my_varID3,tcconak);
        netcdf.putVar(ncid,my_varID4,tcconp);
        netcdf.putVar(ncid,my_varID5,tcconlat);
        netcdf.putVar(ncid,my_varID6,tcconlon);
        netcdf.putVar(ncid,my_varID7,ones(length(tcconid),1));
        netcdf.putVar(ncid,my_varID9,tcconco2);
        netcdf.putVar(ncid,my_varID10,pwf);
	netcdf.putVar(ncid,my_varID11,tcconid);
	netcdf.putVar(ncid,my_varID12,xprior);
	
        % Verify that the variable was created.
        netcdf.close(ncid)

	% End if statements and loops                
	end; % End of day loop
	end; % End of month loop
	end; % End of year loop               

	disp('Netcdf files written.');
	writematrix([(1:(counter-1))' tcconobnall],'tccon_id_information.csv');
	disp('ID written text file');

%-------------------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT
%-------------------------------------------------------------------------------------------------------------------------------%


