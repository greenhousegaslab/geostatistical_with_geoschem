%---------------------------------------------------------------------------------------%
% SCRIPT: make_restart.m								%
% PURPOSE: Create a netcdf restart file from a binary punch restart file.		%
% S. Miller, Dec. 26, 2025								%
%											%
%---------------------------------------------------------------------------------------%

%------------%
% NOTES:     %
%------------%


%---------------%
% Begin script  %
%---------------%

	% Add path to the binary punch matlab files
	addpath('/home/smill191/scripts/OCO2_MIP/rundir_toedit/auxiliary_files/');


	% Read in the existing restart file for 20140901 generated from CarbonTracker
	geosdir = '/scratch4/smill191/smiller/geoschem_adjoint_co2/v11/IS/IS_2015/';
	[ restart ] = readAllBPCHData( strcat(geosdir,'CT.20140901'), strcat(geosdir,'tracerinfo.dat'),strcat(geosdir,'diaginfo.dat'),true,true,false);


%-----------------------------------------------------------------%
% Pull out variables from the existing binary punch restart file  %
%-----------------------------------------------------------------%

	CO2   = restart.C_IJ_AVG.T_CO2.data;
	CO2   = CO2 ./ 1e9; % Convert from ppb to v/v
	lon   = restart.modelData.longC; 
	lat   = restart.modelData.latC;
	time  = restart.C_IJ_AVG.T_CO2.datenum; 
	dims  = restart.C_IJ_AVG.T_CO2.dims;
	layer = 1:dims(3);
	layer = layer';

%-------------------------------%
% Create the new restart file   %
%-------------------------------%

	outname = 'CT.20140901.nc';

        % Open the netcdf file
        ncid = netcdf.create(strcat(geosdir,outname),'NETCDF4');

        % Define the dimensions of the netcdf file
        dimid1 = netcdf.defDim(ncid,'lon',   length(lon));
        dimid2 = netcdf.defDim(ncid,'lat',   length(lat));
        dimid3 = netcdf.defDim(ncid,'layer', length(layer));

	% Add tracer details 
	my_varID   = netcdf.defVar(ncid,'tracer','NC_FLOAT',[dimid1 dimid2,dimid3]);
	my_varID1  = netcdf.defVar(ncid,'lon'     ,'NC_FLOAT',[dimid1]);
        my_varID2  = netcdf.defVar(ncid,'lat'     ,'NC_FLOAT',[dimid2]);
        my_varID3  = netcdf.defVar(ncid,'layer'   ,'NC_FLOAT',[dimid3]);
	netcdf.putAtt(ncid, my_varID, 'long_name', 'GEOS-CHEM Restart File: Instantaneous Tracer Concentrations (v/v)')
        netcdf.putAtt(ncid, my_varID, 'units', 'V/V')
	netcdf.putVar(ncid,my_varID  , CO2);
        netcdf.putVar(ncid,my_varID1 ,single(lon));
        netcdf.putVar(ncid,my_varID2 ,single(lat));
        netcdf.putVar(ncid,my_varID3 ,single(layer));

	% Add global attributes
	varid = netcdf.getConstant("NC_GLOBAL");
	netcdf.putAtt(ncid,varid,'MODELNAME','MERRA2_47L');
	netcdf.putAtt(ncid,varid,'LONRES',restart.modelData.hRes(2));
	netcdf.putAtt(ncid,varid,'LATRES',restart.modelData.hRes(1));
	netcdf.putAtt(ncid,varid,'HALPOLAR',1);
	netcdf.putAtt(ncid,varid,'CENTER180',1);
	netcdf.putAtt(ncid,varid,'CATEGORY','IJ-AVG-$');
	netcdf.putAtt(ncid,varid,'NTRACER',1);
	netcdf.putAtt(ncid,varid,'UNIT','v/v');
        netcdf.putAtt(ncid,varid,'ZTAU0',time(1));
        netcdf.putAtt(ncid,varid,'ZTAU1',time(2));
        netcdf.putAtt(ncid,varid,'RESERVED','');
        netcdf.putAtt(ncid,varid,'NI',dims(1));
        netcdf.putAtt(ncid,varid,'NJ',dims(2));
        netcdf.putAtt(ncid,varid,'NL',dims(3));
        netcdf.putAtt(ncid,varid,'IFIRST',1);
        netcdf.putAtt(ncid,varid,'JFIRST',1);
        netcdf.putAtt(ncid,varid,'LFIRST',1);
        netcdf.putAtt(ncid,varid,'NSKIP',0);

        % Close the netcdf file
        netcdf.close(ncid);


%---------------------------------------------------------------------------------------%
% END OF SCRIPT										%
%---------------------------------------------------------------------------------------%

