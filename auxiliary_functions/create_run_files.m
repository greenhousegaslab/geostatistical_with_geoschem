%-----------------------------------------------------------------------------------------------%
% SCRIPT: create_run_files.m									%
% PURPOSE: Create the "run" and "input.geos" files necessary fo running GEOS-Chem. 		%
% 	The inputs for this function are set in the file inversion_run.m			%
% S. Miller, Oct 24, 2025									%
%-----------------------------------------------------------------------------------------------%

%-------------------%
% FUNCTION NOTES:   %
%-------------------%

	% Note: This script does not work for older versions of Matlab (e.g., only works for R2022a and forward)

	% strcat() will remove trailing spaces. Instead, the number 32 refers to the ascii character for a space.

	% I should set the start and end dates in inversion_run.m
	% In inversion_run.m, I can set the root ExtData directory, and then set a bunch of other directories based on that.
	% I should mention in inversion_run.m that we make some assumptions about how ExtData is organized, so I would double-check
	% the input.geos file the first time you run to make sure it's populated correctly. Otherwise, edit create_run_files.m.

	% A note on the time inputs:
	%	* The startday and endday should be objects created with the Matlab 'juliandate' function.
	%	* rundir: GEOS-Chem run directory
	%	* datadir: The ExtData directory. On Rockfish: /scratch16/smill191/GEOSChem/ExtData/GEOS_4x5/
	%	* landfile: The landmap.mat file (Gives the land and ocean mask.)
	%	* onebyonedir: The folder with 1x1 inputs. On Rockfish: /scratch16/smill191/GEOSChem/ExtData/GEOS_1x1/
	%	* OCO2dir: The folder containing OCO-2 data. On Rockfish: /scratch4/smill191/OCO2_MIP_inputs/v11_obs/obs_v11_LNLG/
	%	* flaskdir: The folder containing flask observations. On Rockfish: /scratch4/smill191/OCO2_MIP_inputs/v11_obs/insitu_obs/
	%	* runname: Name of the current inverse modeling run.
	%	* startday: First day of simulations, created with the 'juliandate' function.
	%	* endday: Last day of simulations, created with the 'juliandate' function.

%------------------%
% Begin function   %
%------------------%

function [ ] = create_run_files(rundir,  datadir, auxdir, landfile, onebyonedir, OCO2dir, flaskdir, runname, startday, endday);


%-------------------------%
% Create the run file     %
%-------------------------%

	outfile = strcat(rundir,'run');

	% Use the run_template files to concatenate and create the run file
	% Add custom inputs to that file
	writelines('#!/bin/bash',outfile,WriteMode='overwrite');
	writelines(strcat('#SBATCH --job-name=',runname),outfile,WriteMode='append');	

	unix(strcat('cat',32,auxdir,'run_template >>',32,outfile));

	writelines(strcat('DRUNDIR=',rundir),outfile,WriteMode='append')

	unix(strcat('cat',32,auxdir,'run_template2 >>',32,outfile));


%-------------------------------------------%
% Determine the number of days in the run   %
%-------------------------------------------%

	numdays = endday - startday;

	dt = datetime(startday, 'ConvertFrom', 'juliandate');

	% 2. Convert the datetime object to a numeric value in YYYYMMDD format
	yyyymmdd_start = convertTo(dt, 'yyyymmdd');

	dt = datetime(endday, 'ConvertFrom', 'juliandate');

        % 2. Convert the datetime object to a numeric value in YYYYMMDD format
        yyyymmdd_end = convertTo(dt, 'yyyymmdd');

%-----------------------------%
% Create the input.geos file  %
%-----------------------------%

	outfile = strcat(rundir,'/input.geos');

	writelines('GEOS-CHEM SIMULATION: Fill in options here for your run',				outfile,WriteMode='overwrite');
	writelines('------------------------+------------------------------------------------------',	outfile,WriteMode='append');
	writelines('%%% SIMULATION MENU %%% :',								outfile,WriteMode='append');
	writelines(strcat('Start YYYYMMDD, HHMMSS  :',32,num2str(yyyymmdd_start),' 000000'),		outfile,WriteMode='append');
	writelines(strcat('End   YYYYMMDD, HHMMSS  :',32,num2str(yyyymmdd_end),' 000000'),		outfile,WriteMode='append');
	writelines('Run directory           : ./',							outfile,WriteMode='append');
	writelines(strcat('Input restart file      : ./restart.',num2str(yyyymmdd_start),'.nc'),        outfile,WriteMode='append');
	writelines('Make new restart file?  : T',							outfile,WriteMode='append');
	writelines('Output restart file(s)  : Restart/restart.YYYYMMDD',				outfile,WriteMode='append');
	writelines(strcat('Root data directory     :',32,datadir),					outfile,WriteMode='append');
	writelines(' => GEOS-FP    subdir   : GEOS_FP/YYYY/MM/',					outfile,WriteMode='append');
	writelines(' => MERRA2     subdir   : MERRA2/YYYY/MM/',						outfile,WriteMode='append');
	writelines(strcat('Land mask file          :',32,landfile),					outfile,WriteMode='append');
	writelines(strcat('Dir w/ 1x1 emissions etc:',32,onebyonedir),					outfile,WriteMode='append');
	writelines('Temporary directory     : ./temp/',							outfile,WriteMode='append');
	writelines(strcat('OCO2 obs directory      :',32,OCO2dir),					outfile,WriteMode='append');
	writelines('Obs file prefix (ncdf)  : oco2_LtCO2_',						outfile,WriteMode='append');
	writelines(strcat('Flask obs directory     :',32,flaskdir),					outfile,WriteMode='append');
	writelines('Obs file prefix (ncdf)  : obspack_co2_1_OCO2MIP_v5.0_2024-10-03.',			outfile,WriteMode='append');
	writelines('In situ flag (0 or 1)   : 1',							outfile,WriteMode='append');
	writelines('Global offsets I0, J0   : 0 0',							outfile,WriteMode='append');
	writelines('Emissions file path     : ./fluxes/',						outfile,WriteMode='append');
	writelines(strcat('Number of time periods  :',32,num2str(numdays)),				outfile,WriteMode='append');


%-----------------------------------------------------------------------------------------------%
% END OF SCRIPT
%-----------------------------------------------------------------------------------------------%

