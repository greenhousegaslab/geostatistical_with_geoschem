%-----------------------------------------------------------------------------------------------------------------------------------------------%
% FUNCTION; cost_gradient_fun															%
% PURPOSE:  Compute the geostatistical inversion cost function. Use the Kronecker product to make these calculations more effecient.		%
%	     This script will calculate both the cost function and the gradient (when the latter is requested).					%
% S. Miller, Sept. 21, 2018															%
% Modified for the GEOS-Chem adjoint on-line modeling runnning											%
%																		%
%-----------------------------------------------------------------------------------------------------------------------------------------------%

%----------------%
% NOTES:         %
%----------------%

	% Note: I assume here that we're running the inverse model on the 4x5 GEOS-Chem grid. If using a different grid, makes sure to modify the 
	% variables "lats" and "lons" below.

%-------------------%
% BEGIN FUNCTION:   %
%-------------------%

function [ f,g ] = cost_gradient_fun( X, A, B, Einv, Dinv, Dinv1, CD, CE, shat, landmap, year, geosdir, fluxpath, ocean);


%-------------------------------------------------------------------%
% Advance the counter and save the flux estimate at each iteration  %
%-------------------------------------------------------------------%

	disp('Calculate the cost and (possibly) gradient functions');

	% Read in the counter and advance it by 1
	load(strcat(fluxpath,'count.mat'),'count');
	count = count +  1;
	save(strcat(fluxpath,'count.mat'),'count');

	% Set a timer to determine how long each iteration takes
	tic;
	
	% Set the dimensions of the flux vector for subsequent calculations
	ntimes = size(Dinv,1);
	m = size(X,1);
	m1 = m ./ ntimes; % Divide m by the total number of time periods in the inversion

%-----------------------------------------------------% 
% Back-transform the current estimate of the fluxes.  %
%-----------------------------------------------------% 

	% shat is not in flux space -- it has been transformed via Eq. 13 in Miller et all (2020, GMD).
	% We need to back-transform into flux space (Eq. 14 in Miller et al) before calculating drift coefficients
	% and before passing the fluxes through the forward model to calculate the cost function and gradient.

	% Calculate chol(Q)*shat
	% I've included a transpose here.
	% CDt = CD';
	% CDt1=CD1';
	Qx = [];

	for j = 1:ntimes;
    		Qx1   = zeros(m1,1);
    		% QQx1  = zeros(m1,1);
    		for i = 1:ntimes;
        	sel = (m1.*(i-1)+1):(i.*m1);
        	Qx1 = Qx1 + shat(sel) .* CD(j,i);
        	% QQx1= QQx1 + shat(sel).* CD1(j,i);
    		end; % End of i loop
    		% Qx1(landmap==0)=QQx1(landmap==0); clear QQx1;
    		temp =  CE' * Qx1; clear Qx1;
    		Qx   =  [Qx; temp];
	end; % End of i loop
	clear Qx1 temp;
	
	shat_backtransform = Qx;
	clear Qx; 

	%--------------------------------------------%
        % Save the flux estimate at each iteration   %
	%--------------------------------------------%

	% Scot: I've decided to save the back-transformed fluxes at each iteration.
	% 	In Zichong's original code, he saves the fluxes in transformed space.
	% Saved fluxes have units of micromol m-2 s-1

        save(strcat(fluxpath,'shat_',num2str(count),'.mat'),'shat_backtransform');

%-------------------------------------------------------------------------------------%
% Estimate the drift coefficients (beta) based on the current estimate of the fluxes. %
%-------------------------------------------------------------------------------------%

	% Note that the estimate for beta will change as the flux estimate is updated at each iteration
	% of the inverse modeling algorithm. The estimate for beta should converge as the flux estimate
	% converges on a solution.

	disp('Estimate the drift coefficients (beta)');
	B2 = [];
	for j = 1:ntimes;
    		A1 = zeros(m1,1);
    		% AA1 = zeros(m1,1);
    		for i = 1:ntimes;
        		sel = (m1.*(i-1)+1):(i.*m1);
        		A1 = A1 + shat_backtransform(sel) .* Dinv(j,i);
        		% AA1= AA1 + shat_backtransform(sel).* Dinv1(j,i);
    		end; % End of i loop
		% A1(landmap==0)=AA1(landmap==0); clear AA1;
    		temp = Einv * A1;
    		B2 = [B2; temp];
	end; % End of i loop
	clear A1 temp;

	% Estimate the drift coefficients using the current estimate for the fluxes
	beta1 = (X' * B) \ (X' * B2);
	clear B2;
	disp('Estimated drift coefficients');
	disp(num2str(beta1));

	%------------------------------------------------------------%
	% Save the estimated drift coefficients from each iteration  %
	%------------------------------------------------------------%

	save(strcat(fluxpath,'beta_',num2str(count),'.mat'),'beta1');


%----------------------------------------------------------------------------------%
% Prepare to pass the back-transformed fluxes through the forward GEOS-Chem model  %
%----------------------------------------------------------------------------------%

	disp('Write fluxes to netcdf file for GEOS-Chem model');

	% Add ocean fluxes back into the flux vector
	landall = repmat(landmap,ntimes,1);
	landall = landall==1;
	temp = ocean;
	temp(landall) = temp(landall) + shat_backtransform;

	% Convert the fluxes to netcdf. The GEOS-Chem adjoint will read in the fluxes in netcdf format.
	temp = temp .* 6.02e+13;
	% temp= shat_backtransform .* 6.02e+13; % convert from micromol m-2 s-1 to  molec/cm2/s, the unit needed for GEOS-Chem
	temp = reshape(temp,72,46,ntimes); % 72 lon and 46 lat

	% One additional, important note: the function "netcdf.create" only works if you are currently in the folder where you intend to write
	% the netcdf file. Hence, I've included a command to "cd" into that folder.
	currentdir = pwd;
	cd(fluxpath);

	% Define latitudes and longitudes for the netcdf file
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


%--------------------------------------------------------%
% Run the fluxes through the forward and adjoint models. %
%--------------------------------------------------------%

	%% Step 2. run the GC forward and GC adjoint model, to get Hs and the gradient

	% Delete the file residuals.csv from the GEOS-Chem run directory.
	% GEOS-Chem will write a new file with residuals for the current flux estimate
	% (instead of appending to the existing residuals.csv file)
	if count > 1; 
		unix(['mv ',geosdir,'residuals_oco2.csv ',geosdir,'residuals_oco2',num2str(count-1),'.csv']); 
                unix(['mv ',geosdir,'residuals_insitu.csv ',geosdir,'residuals_insitu',num2str(count-1),'.csv']); 
                unix(['mv ',geosdir,'residuals_situ_identifier.csv ',geosdir,'residuals_situ_identifier',num2str(count-1),'.csv']); 
	end;
	unix(['rm ',geosdir,'residuals_oco2.csv']);
        unix(['rm ',geosdir,'residuals_insitu.csv']);
        unix(['rm ',geosdir,'residuals_situ_identifier.csv']);


%-----------------------------------------------------%
% OPTION 1: Launch GEOS-Chem from Matlab	      %
%-----------------------------------------------------%

	disp(['GEOS-Chem launched. Current iteration: ',num2str(count-1)]);
	disp(ctime(time ()))
	unix(['./run &> ADJOINT.out']);
	disp('GEOS-Chem model finished running. Continue calculating cost function.')
    	disp(ctime(time ()))


%-----------------------------------------------------%
% OPTION 2: Launch GEOS-Chem as a separate slurm job. %
%-----------------------------------------------------%

%	% Submit the model run -- forward and adjoint
%	% Save the job ID so that we can then check on the job status later
%	[status,cmdout] = unix(['sbatch --parsable ',geosdir,'run']);
%	jobid = str2num(cmdout);
%	disp('GEOS-Chem launched. Job ID:');
%	disp(cmdout);
%
%	% Below is the command to monitor if the GC model is still runing or has completed
%	% The required code may differ on different computers.
%
%	% Check if the model is still running or completed%
%	% This code will check for the status based on the job ID of the GEOS-Chem job submitted above.
%	[status,cmdout]=unix(['sacct --jobs=',num2str(jobid)]); %check the status, please adjust to whatever alias that you use in your supercomputer%
%
%	% Note that "19" refers to 19 characters from the end of "cmdout". You may need to adjust this number
%	% depending on the printed output on your computer system.
%	% (I'd love to find a more elegant way to check on the run status :P.)
%	while cmdout(length(cmdout)-19)~='C'; %C means complete%
%    		disp('GEOS-Chem forward and adjoint model is still running');
%                % disp(cmdout);
%                % disp(cmdout(length(cmdout)-19));
%		% disp(num2str(jobid));
%	 	% disp(['sacct --jobs=',num2str(jobid)]);
%    		pause(180); %pause for three minutes%
%    		[status,cmdout]=unix(['sacct --jobs=',num2str(jobid)]);
%	end; % 
%
%	disp('GEOS-Chem model finished running. Continue calculating cost function.')


%---------------------------------------------------------------------------%
% Read in the cost function and gradient values (after GEOS-Chem finishes)  %
%---------------------------------------------------------------------------%

	%------------------------%
	% Read in the gradient   %
	%------------------------%

	[ adj ] = readBPCHSingle(strcat(geosdir,'OptData/gctm.gdt.01'),'C_IJ_GDE','CO2bal',strcat(geosdir,'tracerinfo.dat'),strcat(geosdir,'diaginfo.dat'),true,true,false);

        adj=adj.*6.02.*(10^13);  % Convert units on the gradient to 1/(umol/m2/s)

	adj=reshape(adj,[],1);

	% Only keep land cells
	adj = adj(landall);


	%-----------------------------%
	% Read in the cost function   %
	%-----------------------------%

	% Read it the cost function
	temp = dlmread(strcat(geosdir,'OptData/cfn.01'));
	L1 = temp(2); 
	clear temp;
	disp('Cost function value of L1');
	disp(num2str(L1));

%----------------------------------------------------%
% Compute the prior component of the cost function   %
%----------------------------------------------------%

        %! In the adjoint code, the observation component of the cost function is multiplied by a factor of (1/2)
        %! (e.g., as in Eq. 4 of Michalak et al. 2004).
	%! E.g., see line 758 of obs_operators/gosat_co2_mod.f
        %! Hence, I will multiply L2 by a factor of 0.5.

	Gshat = shat - A * ((X' * B) \ (A' * shat));
	L2 = 0.5 .* shat' * Gshat;
	clear Gshat;

	disp('Cost function value of L2');
	disp(num2str(L2));

%-------------------------------------------------%
% Add up the two components of the cost function  %
%-------------------------------------------------%

	f = L1+ L2;
	disp('Time elapsed for GEOS-Chem and cost function calculations');
	disp(toc);

	disp('Cost function value');
	disp(f);

	% Save out the cost function information for each iteration
	save(strcat(fluxpath,'cfun_',num2str(count),'.mat'),'L1','L2','f');

%---------------------------------------%
%***** CALCULATE THE GRADIENT      *****%
%---------------------------------------%

if ( nargout > 1 );
    
    	disp('Calculating the gradient');
    
%------------------------------------------------------%
% Calculate the observation component of the gradient  %
%------------------------------------------------------%
    
	L1 = -adj; %to match how Scot builds L1 
    	% (as in https://github.com/greenhousegaslab/geostatistical_inverse_modeling/blob/master/cost_gradient_fun_transform.m)
    
    	% Multiply the variable "L1" by chol(Q)
    	Qx = [];
   
    	for j = 1:ntimes;
        	Qx1   = zeros(m1,1);
        	for i = 1:ntimes;
            		sel = (m1.*(i-1)+1):(i.*m1);
            		Qx1 = Qx1 + L1(sel,:) .* CD(j,i);
        	end; % End of i loop
        	temp =  CE * Qx1;
        	Qx   =  [Qx; temp];
    	end; % End of i loop
    	clear Qx1 temp;
    	L1 = Qx;
    
%------------------------------------------------%
% Calculate the prior component of the gradient  %
%------------------------------------------------%
    
    	% Calculate G * shat
    	% L2 = [I + (sQinv * X) * ((X' * Qinv * X) \ (sQinv * X)') ] * shat
    
    	L2 = shat - A * ((X' * B) \ (A' * shat));
    
%--------------------------------------------%
% Add up the two components of the gradient  %
%--------------------------------------------%
    
    	g = L2 - L1;
   
 	disp('Time elapsed for gradient calculations');
    	disp(toc);
    
    	disp('Gradient mean of L1');
    	disp(num2str(mean(L1)));
    
    	disp('Gradient mean of L2');
    	disp(num2str(mean(L2)));
    
end; % End of nargout if statement


%-----------------------------------------------------------------------------------------------------------------------------------------------%
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------------------------------------------------------%
