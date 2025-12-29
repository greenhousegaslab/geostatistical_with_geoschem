%-----------------------------------------------------------------------------------------------------------------------%
% SCRIPT:  inversion_run												%
% PURPOSE: Launch a geostatistical inverse model (GIM), implemented with Lagrange multipliers for non-negativity	%
% S. Miller, Jan. 8, 2016												%	
% Modified to be used to couple GEOS-Chem forward and adjoint model							%
% Z. Chen, Dec. 18, 2018												%
% Re-organized the code and added the KronMat object.									%
% S. Miller, Dec. 23, 2025												%
%															%
%-----------------------------------------------------------------------------------------------------------------------%


%-----------%
% NOTES:    %
%-----------%

	% Note that the format of the fluxes in this script is a little bit different from the format of fluxes used in
	% GEOS-Chem.
	% Notably, this script is formatted after the inverse modeling scripts that we use with STILT. The Matlab scripts
	% process fluxes in units of micromol m-2 s-1 and assume that the rows of the flux matrix are latitudes while the
	% columns are longitude.
	% GEOS-Chem, by contrast, processes fluxes in units of molec/cm2/s and assumes that the rows of the flux matrix
	% are longitudes (while the columns are latitudes). 
   
	% There are numerous different options that you can set or change in this file, and I would read through this file
	% in detail before launching the inverse modeling run.

	% Note that this version of the code will automatically create the "run" and "input.geos" files requied by GOES-Chem.
	% Also, nothing needs to be changed in the input.gcadj file.


%-------------------------------%
% Set required script options   %
%-------------------------------%
    
	% Set the name of the run
	runname = 'LNLG_2015';

	% Set the year of the inverse modeling simulations associated with this specific run.
	% Note that 2015 simulations include a few months in 2014
	yr = 2015;

	%---------------------------%		
	% Set required file paths   %
	%---------------------------%

	% Path to the GEOS-Chem run directory (This directory should include the file 'run' which launches GEOS-Chem.)
	geosdir = '/scratch4/smill191/kmorgan/OCO2_MIP_v11/IS/IS_2016/';

	% Path where the fluxes should be saved at each iteration (so the fluxes can be read in by GEOS-Chem).
	% Note that this path does not contain the final inversion outputs but rather is a temporary file path
	% for storing fluxes at each iteration.
	fluxpath = strcat(geosdir,'fluxes/');

	% Folder with all the auxiliary scripts and files required for the runs
	% These files include the land mask ('landmap.mat'), the Matlab scripts that read in binary punch files, and the 
	% files that are used to create the input.geos and run files necessary for running GEOS0-Chem.
	% fminlbfgs.m must also be in this folder.
	% Functions that read binary punch files are avialable for download here: 
	% http://wiki.seas.harvard.edu/geos-chem/index.php/Matlab_software_tools_for_use_with_GEOS-Chem
	auxdir = '~/scripts/utilities/co2_inversion_functions/';
	addpath(auxdir);

	% Set other directories required for GEOS-Chem
	extdir      = '/scratch16/smill191/GEOSChem/ExtData/';
	datadir     = strcat(extdir,'GEOS_4x5/');
	onebyonedir = strcat(extdir,'GEOS_1x1/');
        landfile    = strcat(datadir,'MERRA2/2015/01/MERRA2.20150101.CN.4x5.nc4');      
        
	% Path to the OCO2 observations:
	OCO2dir   = '/scratch4/smill191/OCO2_MIP_inputs/v11_obs/obs_v11_LNLG/';

	% Path to the in situ observations:
	flaskdir = '/scratch4/smill191/OCO2_MIP_inputs/v11_obs/insitu_obs/';


	%---------------------------------------------%
	% **Set the covariance matrix parameters**    %
        % Please note the unit for elements in theta  % 
        % are ppm, umol/m2/s, km, and days.           %                             
	%---------------------------------------------%

	% Note that the model-data mismatch error is defined in the observation netcdf files, not in the Matlab scripts here.     
     	% theta(1): defines the diagonal elements of Q
	% theta(2): decorrelation length
	% theta(3): decorrelation time
	theta = [ sqrt(0.31) 500 5.1 ];  % for land%

	% Display the covariance matrix parameters on screen
	disp('Covariance matrix parameters for land');
	disp(theta);


%----------------------------------------------------------------%
% Check to see if this run is a continuation of an existing run  %
%----------------------------------------------------------------%

        % If this run is the continuation of an existing run, there will be a 
        % file in fluxdir called "count.mat" that contains the iteration counter
        % from the previous run. 

        % If that file exists, we'll assume that this inversion run is the continuation
        % of a previous run and will use the last available iteration as the initial guess
        % in this launch of L-BFGS.
        restartrun = 0;
        if exist(strcat(fluxpath,'count.mat')) == 2;
                restartrun = 1;
	end;
	
%---------------------------------%
% **Read in the observations**    %
%---------------------------------%

	% Note: in an adjoint-based inverse model, the observations are read in by the GEOS-Chem adjoint fortran code.
	% Hence, we will not read in the observations here.
 
%----------------------------%
% **Create the X matrix**    %
%----------------------------%

	disp('Create the X matrix');

	% Note from Scot: for the MIP studies, I'll create a single column of X that includes
	% biospheric fluxes, ODIAC fossil fuel emissions, GFED biomass burning emissions,
	% and ocean fluxes from NASA

	[ X, ocean ] = read_fluxes(yr);

	% If yr is 2015, then we'll read in the last few months of 2014 and append to 2015
	if yr==2015;
		[ X2014, ocean2014 ] = read_fluxes(2014);
		X                    = [X2014; X];
		ocean                = [ocean2014; ocean];
	end;

        disp('X has been successfully processed');

	%--------------------------------%
	% Set dimensions and land mask   %
	%--------------------------------%

	% The current landmap file included with this code is for a global 4x5 degree run.
	% You'll need to change the landmap file if you want to run a different domain or spatial resolution.
	% If you do run a different domain, you'll also need to modify the function write_fluxes.m to match your new domain.

        disp('Set land map and inverse modeling dimensions');
	load('landmap.mat');
        landmap=reshape(landmap,[],1); % Reshape the land map to be a vector
	ntimes = size(X,1) ./ length(landmap);
	m = length(find(landmap==1)) .* ntimes; 

	landsel = landmap == 1;

        landall = repmat(landmap,ntimes,1);
        landall = landall==1;

        % Only keep land grid cells in X
        X = X(landall,:);


%----------------------------------------------%
% **Read in the geographic distance matrix**   %
%----------------------------------------------%

        % The deltamat file here is for a global 4x5 degree run (land only). You'll need to create a different deltamat file if you want to run
        % the inverse model at a different spatial resolution or a different domain.

        % This matrix (dimensions r x r) should define the distance between each emissions grid box
        % This matrix is symmetric with zeros on the diagonals. 
        % This matrix should typically has units of km. Whatever units you choose should match the units of theta(2)
        % Refer to Eq. 7 in Gourdji et al., 2010 for detail.
        load('deltamat.mat','deltamat');


%-------------------------%
% Run the inverse model   %
%-------------------------%

	inverse_model;


%------------------------------%
% Write the outputs to file    %
%------------------------------%

        disp('Writing outputs to file');
        outname = strcat(fluxpath,'fluxes_LBFGS.csv');
        dlmwrite(outname,full(shat),',');

        disp('Outputs written to file');
        disp(outname);


%-----------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT														%
%-----------------------------------------------------------------------------------------------------------------------%
