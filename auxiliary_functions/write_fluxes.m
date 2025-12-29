%-----------------------------------------------------------------------------------------------------------------------------------------------%
% SCRIPT: write_fluxes.m															%
% PURPOSE: Write a vector of fluxes from the inverse model for GEOS-Chem									%
% S. Miller, Dec 24, 2025															%
%																		%
%-----------------------------------------------------------------------------------------------------------------------------------------------%

	%-----------%
	% NOTES:    %
	%-----------%

		% This script is written for a global 4x5 degreee run.
		% This script needs to be modified if you're running a different spatial resolution or domain.

%-------------------%
% BEGIN FUNCTION:   %
%-------------------%

function [ ] = write_fluxes(fluxpath, landall, ocean, ntimes, shat_backtransform, year);


        disp('Write fluxes to netcdf file for GEOS-Chem model');

        % Add ocean fluxes back into the flux vector
        temp = ocean;
        temp(landall) = temp(landall) + shat_backtransform;

        % Convert the fluxes to netcdf. The GEOS-Chem adjoint will read in the fluxes in netcdf format.
        temp = temp .* 6.02e+13;
        % temp= shat_backtransform .* 6.02e+13; % convert from micromol m-2 s-1 to  molec/cm2/s, the unit needed for GEOS-Chem
        %! NOTE!!! This line needs to be changed if you're running any simulation other than 4x5 degree global
        temp = reshape(temp,72,46,ntimes); % 72 lon and 46 lat

        % One additional, important note: the function "netcdf.create" only works if you are currently in the folder where you intend to write
        % the netcdf file. Hence, I've included a command to "cd" into that folder.
        currentdir = pwd;
        cd(fluxpath);

        % Define latitudes and longitudes for the netcdf file
        %! NOTE!!! This line needs to be changed if you're running any simulation other than 4x5 degree global
        lons = -180:5:175;
        lats = [-89 -86:4:86 89];

        for j=1:ntimes,

        % Convert day of year to month and day of month
        if year > 2015; [yy month day HH MM] = datevec(datenum(year,1,j)); end;

        % If year=2015, then the inverse model starts Sept 1, 2014. Adjust the dates accordingly
        if year==2015;
                if j <= 122;  [yy month day HH MM] = datevec(datenum(2014,9,j));   end;
                if j >  122;  [yy month day HH MM] = datevec(datenum(year,1,j-122)); end;
        end; % End of if statement

        if month<10; month1=strcat('0',num2str(month)); else; month1=num2str(month); end;
        if day<10; day1=strcat('0',num2str(day)); else; day1=num2str(day); end;

        % Open the netcdf file
        ncid = netcdf.create(strcat('CO2.daily.geos.4x5.',num2str(yy),".",month1,".",day1,'.nc'),'NETCDF4');

        % Define the dimensions of the netcdf file
        dimid1 = netcdf.defDim(ncid,'lon', length(lons));
        dimid2= netcdf.defDim(ncid,'lat', length(lats));

        % Define new variables for the netcdf file
        my_varID = netcdf.defVar(ncid,'CO2_flux','double',[dimid1 dimid2]);
        my_varID1= netcdf.defVar(ncid,'lon','double',[dimid1]);
        my_varID2= netcdf.defVar(ncid,'lat','double',[dimid2]);
        netcdf.endDef(ncid);

        % Write data to the netcdf variable
        netcdf.putVar(ncid,my_varID,temp(:,:,j));
        netcdf.putVar(ncid,my_varID1,lons);
        netcdf.putVar(ncid,my_varID2,lats);

        % Close the netcdf file
        netcdf.close(ncid);

        end; % End of j loop
        clear temp;

        % Move from the flux directory back to the previous directory
        cd(currentdir);


%-----------------------------------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT
%-----------------------------------------------------------------------------------------------------------------------------------------------%

