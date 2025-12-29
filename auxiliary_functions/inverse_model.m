%-----------------------------------------------------------------------------------------------------------------------%
% SCRIPT: inverse_model													% 
% PURPOSE: Creates the inputs required for the inverse model and launches the L-BFGS optimization algorithm.		%
% S. Miller, Dec. 24, 2025												%
%															%
%-----------------------------------------------------------------------------------------------------------------------%


%-------------------------%
% Create the E matrix     %
%-------------------------%

	disp('Create the E and D matrices');

	% The E and D matrices are used to construct the Q covariance matrix. 
	% The E matrix describes how covariances in Q decay in space, and the D matrix describes how covariances in Q decay in time.
	% See the following paper by Vineet Yadav for more details on the E and D matrices: http://www.geosci-model-dev.net/6/583/2013/gmd-6-583-2013.html
	% The E matrix has dimensions (r x r)
	% The D matrix has dimensions (q x q)
	
	%------------%
	% Create E   %
	%------------%

    	sigmaQ = repmat(theta(1),length(landmap),1);	

	E = exp(-1 .* deltamat ./ theta(2));
	
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

	D = exp(-1 .* days ./ theta(3));  Dinv =inv(D);


%----------------------------------------------%
% Create Q and Q^(-1) using the kronMat class  %
%----------------------------------------------%

	Q    = kronMat(D,E);
	Qinv = kronMat(Dinv,Einv);

	%----------------------------------------%
	% Take the symmetric sq root of D and E  %
	%----------------------------------------%

	CD = sqrtm(D);
	CE = sqrtm(E);
	invCD = inv(CD);
	invCE = inv(CE);

	Qsqrt    = kronMat(CD,CE);
	Qsqrtinv = kronMat(invCD,invCE);


%------------------------%
% Create the R matrix    %
%------------------------%
	
	% Note that the R matrix is calculated as part of the GEOS-Chem adjoint code.
	% The model-data mismatch errors (i.e., the diagonal elements of R) need to be specified in the
	% observation files that are read into GEOS-Chem. These values are not needed by the Matlab code.


%--------------------------------------------------------------------%
% Create files for GEOS-Chem with run options (input.geos and run)   %
%--------------------------------------------------------------------%

	disp('Create the GEOS-Chem run and input.geos files');

        % By default, we'll assume that this run covers a full year (except for 2015, which starts in Sept 2014)
        if(yr > 2015);
                startday = juliandate(yr,1,1);
                endday   = juliandate(yr+1,1);
        end;
        if(yr == 2015);
                startday = juliandate(2014,9,1);
                endday   = juliandate(2016,1,1);
        end

        % Create the run files required by GEOS-Chem
        create_run_files(geosdir, datadir, auxdir, landfile, onebyonedir, OCO2dir, flaskdir, runname, startday, endday);


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

        shat0 = X;

        % Create an iteration counter
        count=0;
        save(strcat(fluxpath,'count.mat'),'count');
        clear count;
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

        load(strcat(fluxpath,'count.mat'),'count');
        load(strcat(fluxpath,'shat_',num2str(count),'.mat'),'shat_backtransform');
        shat0 = shat_backtransform;

        disp('Existing value of the iteration counter:');
        disp(num2str(count));

        % cost_gradient_fun.m will +1 to count at the beginning of the script. 
        % Here, I subtract 1 from count so that we start cost_gradient_fun.m
        % on the same counter where we left off at the end of the last run.
        count = count -  1;
        save(strcat(fluxpath,'count.mat'),'count');

        end;


%-------------------------------------------------%
% Pre-calculate matrix products where possible    %
%-------------------------------------------------%

	disp('Calculate inv(Q) * X');
	% We'll refer to this matrix product as varible "B" from now on
	B = Qinv * X;

        disp('Calculate inv(sqrtm(Q)) * X');   
  	A = Qsqrtinv * X; 

	% Transform shat0: shat0* = Q^(-1/2)*shat0
	shat0 = Qsqrtinv * shat0;


%-------------------------%
% Estimate the fluxes     %
%-------------------------%

	disp('Run inversion (optimized with the L-BFGS algorithm)');

	% Set maximum number of iterations
	maxit = 50;

	% Define the script that will calculate the cost function and gradient. 
	f1 = @(shat) cost_gradient_fun(X, A, B, Qinv, Qsqrt, shat, landall, yr, geosdir, fluxpath, ocean, ntimes); 

        % Create an empty flux estimate
        shat = [];

        options = struct('HessUpdate','lbfgs','GradObj','on','Display','iter','MaxIter',maxit,'GradConstr',false);

        [shat,costfun,exitflag,gradient] = fminlbfgs(f1,shat0,options);


%-----------------------------------------------%
% Transform fluxes back to normal space         %
%-----------------------------------------------%

	shat = Qsqrt * shat;	


%-----------------------------------------------------------------------------------------------------------------------%
% END OF SCRIPT														%
%-----------------------------------------------------------------------------------------------------------------------%
