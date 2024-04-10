%-----------------------------------------------------------------------------------------------------------------------%
% SCRIPT:  inversion_run												%
% PURPOSE: Launch a geostatistical inverse model (GIM), implemented with Lagrange multipliers for non-negativity	%
% S. Miller, Jan. 8, 2016	
% Modified to be used to couple GEOS-Chem forward and adjoint model
% Z. Chen, Dec. 18, 2018
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
	% In the future, we could reformat this inverse modeling code to match GEOS-Chem
    
%-------------------------------%
% Set required script options   %
%-------------------------------%
    
     	% Set the maximum number of iterations that are allowed for the L-BFGS algorithm
	maxit = 50;

	% Set the year of the inverse modeling simulations associated with this specific run.
	% Note that 2015 simulations include a few months in 2014
	yr = 2021;

	%---------------------------%		
	% Set required file paths   %
	%---------------------------%

	% Create a path that contains other temporary inputs and outputs.
	% These include deltamat, which lists the distance in km among different GEOS-Chem model grid boxes.
	% This folder should contain the prior fluxes (in this case from ORCHIDEE, GFED, Ocean fluxes, and ODIAC.)
	% Lastly, the land mask is saved there.
        inpath = '/scratch4/smill191/smiller/data/OCO2_MIP/misc_inputs/';

	% Path to the GEOS-Chem run directory (This directory should include the file 'run' which launches GEOS-Chem.)
	geosdir = '/scratch4/smill191/smiller/geoschem_adjoint_co2/LNLG_2021/';

	% Path where the fluxes should be saved at each iteration (so the fluxes can be read in by GEOS-Chem).
	% Note that this path does not contain the final inversion outputs but rather is a temporary file path
	% for storing fluxes at each iteration.
	fluxpath = '/scratch4/smill191/smiller/outputs/OCO2_MIP/fluxes_LNLG_2021/';

	% Path to the Matlab functions that will read in binary punch files
	% Functions are avialable for download here: http://wiki.seas.harvard.edu/geos-chem/index.php/Matlab_software_tools_for_use_with_GEOS-Chem
	addpath('~/scripts/utilities/BPCH_functions/');

	% Note that the function fminlbfgs.m must either be in the same folder as this script, or you must add the path to that function.
	
	%---------------------------------------------%
	% **Set the covariance matrix parameters**    %
        % Please note the unit for elements in theta  % 
        % are ppm, umol/m2/s, km, and days.           %                             
	%---------------------------------------------%

	% Note that the model-data mismatch error is defined in the observation netcdf files, not in the Matlab scripts here.     
     	% theta(1): defines the diagonal elements of Q
	% theta(2): decorrelation length
	% theta(3): decorrelation time
	theta = [ sqrt(0.31) 1460 5.1 ];  % for land%

	% Display the covariance matrix parameters on screen
	disp('Covariance matrix parameters for land');
	disp(theta);


%----------------------------------------------------------------%
% Check to see if this run is a continuation of an existing run  %
%----------------------------------------------------------------%

        % If this run is the continuation of an existing run, there will be a 
        % file in fluxdir called "data.mat" that contains the output of the L-BFGS
	% algorithm.

        % If that file exists, we'll assume that this inversion run is the continuation
        % of a previous run and will use the last available iteration as the initial guess
        % in this launch of L-BFGS.
        restartrun = 0;
        if exist(strcat(fluxpath,'data.mat')) == 2;
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

	% Note: I've hard-coded the file paths to the fluxes below. Change these paths as needed for
	% your inverse modeling setup.

	% Note from Scot: for the MIP studies, I'll create a single column of X that includes
	% biospheric fluxes, ODIAC fossil fuel emissions, GFED biomass burning emissions,
	% and ocean fluxes from NASA

	% Some flux products aren't available after 2018. Reset the flux year to 2018 where needed
	yr1 = yr;
	yr2 = yr;
	if yr > 2018; yr1 = 2018; end;
	if yr > 2017; yr2 = 2017; end;

	% Read in CASA fluxes (in micromol m-2 s-1)
	% CASA is only available through 2018
	load(strcat('/scratch4/smill191/smiller/data/OCO2_MIP/CASA_GFED4/regridded/','CASA_',num2str(yr1),'.mat'));

	% Read in ODIAC fossil fuel emissions
	load(strcat('/scratch4/smill191/smiller/data/OCO2_MIP/ODIAC/regridded/','ODIAC_',num2str(yr),'.mat'));

	% Read in GFED biomass burning emissions
	load(strcat('/scratch4/smill191/smiller/data/OCO2_MIP/GFED/regridded/','gfed_',num2str(yr),'.mat'));

	% Read in ocean fluxes 	
	load(strcat('/scratch4/smill191/smiller/data/OCO2_MIP/takahashi/','ocean_',num2str(yr2),'.mat'));
	% Convert ocean fluxes from single precision to double precision
	ocean = double(ocean);

	%-------------------------------------------------------%
	% If year=2015, add fluxes for last few months of 2014  %
	%-------------------------------------------------------%

        % The inverse model for year 2015 also includes fluxes for Sept - Dec 2014.
        % If year=2015, I need to modify the flux inputs to include fluxes for the last few months of 2014

	if yr == 2015; 

	% CASA: Existing flux files include a full year of fluxes for 2014
	NEE2014 = load('/scratch4/smill191/smiller/data/OCO2_MIP/CASA_GFED4/regridded/CASA_2014.mat');
	NEE2014 = NEE2014.NEE;

	% Select out elements corresponding to Sept 1 - Dec 31, 2014 and add to existing vector of fluxes	
	sel = (1 + 3312.*(31 + 28 + 31 + 30 + 31 + 30 + 31 + 31)):length(NEE2014);
	NEE = [NEE2014(sel); NEE]; 

	% ECO-Darwin: Existing flux files do not include 2014 (only 2015)
	ocean = [ocean(sel); ocean];

	% ODIAC: Existing flux files only include Sept - Dec 2014
	odiac2014 = load(strcat('/scratch4/smill191/smiller/data/OCO2_MIP/ODIAC/regridded/','ODIAC_2014.mat'));
	odiac = [odiac2014.odiac; odiac]; 

	% GFED: Existing flux files include a full year of fluxes for 2014
	gfed2014 = load(strcat('/scratch4/smill191/smiller/data/OCO2_MIP/GFED/regridded/','gfed_2014.mat'));
	gfed2014 = gfed2014.gfed;
	gfed = [gfed2014(sel); gfed]; 

	end; % End of yr if statement

        %-----------------------------------------------------------%
        % If year=2000, add an extra day to account for leap year   %
        %-----------------------------------------------------------%

        if yr==2020;

        % Feb 28 is the 59th day of the year
        % I want to repeat the 59th day in NEE and ocean (other flux types account for leap year already)

        m1   = length(NEE)./365;
        sel  = ((59-1).*m1 + 1):(59.*m1);
        sela = 1:(59.*m1);
        selb = (59.*m1 +1):length(NEE);

        NEE   = [ NEE(sela); NEE(sel); NEE(selb) ];
        ocean = [ ocean(sela); ocean(sel); ocean(selb) ];

        end; % End of yr if statement

	%------------------------%
	% Assemble the X matrix  %
	%------------------------%

	% Add the differnet flux types together
	% X = NEE + odiac + gfed + ocean;
	X = NEE + odiac + gfed;

        disp('X has been successfully processed');


	%--------------------------------%
	% Set dimensions and land mask   %
	%--------------------------------%

        disp('Set land map and inverse modeling dimensions');
	load(strcat(inpath,'landmap.mat'));
        landmap=reshape(landmap,[],1); % Reshape the land map to be a vector
        % landmap(landmap==2)=nan; %ice is not considered% 
        % In Zichong's landmap.mat matrix, there are no model grid boxes with a designation of "2" (i.e., no boxes labelled as ice). 
	ntimes = size(X,1) ./ length(landmap);
    	% m=3312 * 365; %46 (lat) \times 72 (lon) \times 365 (days)%
	m = length(find(landmap==1)) .* ntimes; 

	landsel = landmap == 1;

        landall = repmat(landmap,ntimes,1);
        landall = landall==1;

        % Only keep land grid cells in X
        X = X(landall,:);


%-------------------------%
% Create the E matrix     %
%-------------------------%

	disp('Create the E and D matrices');

	% The E and D matrices are used to construct the Q covariance matrix. 
	% The E matrix describes how covariances in Q decay in space, and the D matrix describes how covariances in Q decay in time.
	% See the following paper by Vineet Yadav for more details on the E and D matrices: http://www.geosci-model-dev.net/6/583/2013/gmd-6-583-2013.html
	% The E matrix has dimensions (r x r)
	% The D matrix has dimensions (q x q)
	
	%----------------------------------------------%
	% **Read in the geographic distance matrix**   %
	%----------------------------------------------%
	
	% This matrix (dimensions r x r) should define the distance between each emissions grid box
	% This matrix is symmetric with zeros on the diagonals. 
	% This matrix should typically has units of km. Whatever units you choose should match the units of theta(2)
	% Refer to Eq. 7 in Gourdji et al., 2010 for detail.
	load(strcat(inpath,'deltamat.mat'),'deltamat');

	%------------%
	% Create E   %
	%------------%

%----------------------------%
% Create the sigmaQ vector   %
%----------------------------%

    	sigmaQ = repmat(theta(1),length(landmap),1);	
%    sigmaQ=zeros(size(landmap));
%    sigmaQ(landmap==1)=theta(1); %% for land%%
%    sigmaQ(landmap==0)=theta1(1); %% for ocean%%

	E = exp(-1 .* deltamat ./ theta(2));
	
%    E=zeros(3312,3312);
%    for i=1:3312,
%        for j=1:3312,
%            if landmap(i)+ landmap(j)==2, % when land meets land%
%		% Spherical model
%        	%       E(i,j)                      = 1 - 1.5 .* (deltamat(i,j) ./theta(2))  + 0.5 .* (deltamat(i,j).^3 ./ theta(2).^3);
%        	%        temp=E(i,j); temp(deltamat(i,j) > theta(2)) = 0; E(i,j)=temp; 
%
%		% Exponential model
%               	E(i,j) = exp(-1.*(deltamat(i,j)./theta(2)));  %exponential%
%              
%            elseif landmap(i)+ landmap(j)==0, %when ocean meets ocean%
%
%		% Spherical mdoel
%        	%      E(i,j)                      = 1 - 1.5 .* (deltamat(i,j) ./theta1(2))  + 0.5 .* (deltamat(i,j).^3 ./ theta1(2).^3);
%        	%      temp=E(i,j); temp(deltamat(i,j) > theta1(2)) = 0; E(i,j)=temp;
%              
%		% Exponential model
%              	E(i,j) = exp(-1.*(deltamat(i,j)./theta1(2)));
%            end
%        end
%    end
%    %%note currently we do not consiser the spatial covariance between land and ocean %%  
%     clear temp;
%      disp('E has been successfully processed'); 
      
    E = (sigmaQ*sigmaQ') .* E;  
 
    % Only keep land grid cells
    landsel = landmap == 1;
    E = E(landsel,landsel);
     
    % Take the inverse of the E matrix
    % Note: this step will produce an error or will produce values of "Inf" if E is not a positive definite matrix.
    Einv  = inv(E);

	
%----------------------------%
% Create the D matrix        %
%----------------------------%

	days = 1:ntimes;
	days = days';
	days = days * ones(1,length(days));
	days = abs(days - days');

	%------------%
	% Create D   %
	%------------%
	
    	% Option I
	% Spherical covariance model
	% Advantage: it tapers off to zero and is therefore faster to compute with large
	% datasets, e.g., Miller 2018 & 2020;
    	%D = 1 - 1.5 .* (days   ./theta(3))  + 0.5 .* (days.^3   ./ theta(3).^3);
    	%	D(days   > theta(3)) = 0;
    
    	% Option II
    	% Exponential, but pay particular attention to the fact that how we
    	% define decorrelation length and time  differently between exponential and
    	% spherical covariance models. 
    	D = exp(-1 .* days ./ theta(3));  Dinv =inv(D);
    	D1 = exp(-1 .* days ./ theta(3));  Dinv1 =inv(D1);

%------------------------%
% Create the R matrix    %
%------------------------%
	
	% Note that the R matrix is calculated as part of the GEOS-Chem adjoint code.
	% The model-data mismatch errors (i.e., the diagonal elements of R) need to be specified in the
	% observation files that are read into GEOS-Chem. These values are not needed by the Matlab code.
	

%------------------------------------------------------%
% Create the initial guess for the L-BFGS-B algorithm  %
%------------------------------------------------------%

        disp('Create an initial guess for the L-BFGS-B algorithm');

        % The L-BFGS-B algorithm requires an initial guess of the fluxes/emissions. It will then iterate toward the solution.
        % The initial guess can be very important, and a more accurate initial guess is always better than a poor initial guess.
        % The L-BFGS-B will converge quickly on the solution for flux elements that have a large impact on the cost function.
        % However, L-BFGS-B will converge slowly for flux elements that do not have a large impact on the cost function.
        % Here, I will use the deterministic model as the initial guess of the fluxes.

        % Use the prior flux estimate as the initial guess.
        % Note: in this case, I won't worry about estimating beta 
        % (then I would need to run all of the prior fluxes through the GEOS-Chem forward model).

        if restartrun == 0;

        disp('This run is a new run. Setting initial flux estimate equal to prior fluxes.');

        shat0 = NEE + odiac + gfed; %  + ocean;
	shat0 = shat0(landall);

        end;

        %-------------------------------------------------%
        % Read in existing flux estimate (if available)   %
        %-------------------------------------------------%

        % If this run is the restart of an existing run, we'll use the existing flux estimate
        % from the last run as the initial guess of the fluxes.
        % We'll also re-start the counter where we left off (instead of starting it at 0).

        if restartrun == 1;

        disp('This run is the restart of an existing run.');
        disp('Reading in existing flux estimate as the initial flux estimate.');

	% Read in the data.mat file and pull out the outputs from the previous iteration of L-BFGS

        % load(strcat(fluxpath,'count.mat'),'count');
        % load(strcat(fluxpath,'shat_',num2str(count),'.mat'),'shat_backtransform');
        % shat0 = shat_backtransform;

        disp('Existing value of the iteration counter:');
        disp(num2str(count));

        % cost_gradient_fun.m will +1 to count at the beginning of the script. 
        % Here, I subtract 1 from count so that we start cost_gradient_fun.m
        % on the same counter where we left off at the end of the last run.
        % count = count -  1;
        % save(strcat(fluxpath,'count.mat'),'count');

        end;


%-------------------------------------------------%
% Pre-calculate matrix products where possible    %
%-------------------------------------------------%

	disp('Calculate inv(Q) * X');
	% We'll refer to this matrix product as varible "B" from now on
	% B = inv(Q) * X;
    
	B = [];
	m1  = size(E,1);
        % landmap=repmat(landmap,1,size(X,2));        

	for j = 1:size(Dinv,1);
	B1 = zeros(m1,size(X,2));
        % B2 = zeros(m1,size(X,2));
		for i = 1:size(Dinv,1);
		sel = (m1.*(i-1)+1):(i.*m1);
		B1 = B1 + X(sel,:) .* Dinv(j,i);
                % B2 = B2 + X(sel,:) .* Dinv1(j,i);
		end; % End of i loop
        % B1(landmap==0) = B2(landmap==0); clear B2;
	temp = Einv * B1; clear B1;
	B = [B; temp];
	end; % End of j loop
	clear B1 temp;

    	disp('Calculate several matrix products that will be used in the preconditioner');
	CD = sqrtm(D); % CD1=sqrtm(D1);
	CE = sqrtm(E);
	invCD = inv(CD); % invCD1=inv(CD1);
	invCE = inv(CE);

        disp('Calculate inv(sqrtm(Q)) * X');   
        A = [];

        for j = 1:ntimes;
        A1 = zeros(m1,size(X,2));
        % AA1=zeros(m1,size(X,2));
                for i = 1:size(Dinv,1);
                sel = (m1.*(i-1)+1):(i.*m1);
                A1 = A1 + X(sel,:) .* invCD(j,i);
                % AA1 = AA1 + X(sel,:) .* invCD1(j,i); 
                end; % End of i loop
        % A1(landmap==0)=AA1(landmap==0); clear AA1;  
        temp = invCE * A1;  clear A1;
        A = [A; temp];
        end; % End of i loop
        clear A1 temp;	
   

        % landmap=landmap(:,1);
 	%transform shat0 to shat0*%
        Qx=[];
        for j = 1:ntimes;
        Qx1   = zeros(m1,1);
        % QQx1  =zeros(m1,1);
               for i = 1:ntimes;
               sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + shat0(sel) .* invCD(j,i);
                % QQx1 = QQx1 + shat0(sel) .* invCD1(j,i);
                end; % End of i loop
        % Qx1(landmap==0)=QQx1(landmap==0); clear QQx1;
        temp = invCE * Qx1;   clear QX1;
        Qx   =  [Qx; temp];
        end; % End of i loop
        clear Qx1 temp;
        shat0=Qx; clear Qx;

	
%-------------------------%
% Estimate the fluxes     %
%-------------------------%

	disp('Run inversion (optimized with the L-BFGS algorithm)');

	% Define the script that will calculate the cost function and gradient. 
	f1 = @(shat) cost_gradient_fun(X, A, B, Einv, Dinv, Dinv1, CD, CE, shat, landmap, yr, geosdir, fluxpath, ocean); 

        % Create an empty flux estimate
        shat = [];

        if restartrun == 0; 
		options = struct('HessUpdate','lbfgs','GradObj','on','Display','iter','MaxIter',maxit,'GradConstr',false,'restart',false);
                [shat,costfun,exitflag,gradient] = fminlbfgs(f1,shat0,options);
	end;
 	if restartrun == 1;
                % If this run is a restart run, then the intial guess doesn't matter. L-BFGS will grab all needed inputs from data.mat
		options = struct('HessUpdate','lbfgs','GradObj','on','Display','iter','MaxIter',maxit,'GradConstr',false,'restart',true);
                [shat,costfun,exitflag,gradient] = fminlbfgs(f1,ones(m,1),options);
  	end;


%-----------------------------------------------%
% Transform fluxes back to normal space         %
%-----------------------------------------------%

	stemp = [];
 
 	ntimes = size(D,1);
         for j = 1:ntimes;
         Qx1   = zeros(m1,1);
         % QQx1  = zeros(m1,1);
                 for i = 1:ntimes;
                 sel = (m1.*(i-1)+1):(i.*m1);
                Qx1 = Qx1 + shat(sel) .* CD(j,i);
               % QQx1 = QQx1 + shat(sel) .* CD1(j,i);
                end; % End of i loop
        % Qx1(landmap==0)=QQx1(landmap==0); clear QQx1;      
        temp =  CE * Qx1; clear Qx1;
        stemp   =  [stemp; temp];
        end; % End of i loop
        clear Qx1 temp;

	shat = stemp;
	
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
