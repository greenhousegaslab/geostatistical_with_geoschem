#---------------------------------------------------------------------------------------------------------------#
# FUNCTION: cat_tccon.r												#
# PURPOSE: Concatinate CSV files that contain the model-data residuals for TCCON.				#
#	This script will concatinate residuals with identifier and will concatinate years.			#
# S. Miller, Jan 17, 2022											#
#														#
#---------------------------------------------------------------------------------------------------------------#


#-----------#
# NOTES:    #
#-----------#


#------------------------#
# User-required inputs   #
#------------------------#

        # Path to inverse modeling folders
        inpath <- "/scratch/groups/smill191/smiller/geoschem_adjoint_co2/"

        # Path where the OCO-2 files should be written
        outpath <- "/scratch/groups/smill191/smiller/outputs/OCO2_MIP/for_submission/"


#-------------------------------------#
# Read in TCCON obs ID information    #
#-------------------------------------#

        # The TCCON observation ID is a 16-digit integer, but Fortran will only write 
        # integers to 15-digit precision. Hence, I created my own observation identifier
        # that's shorter than 16 digits and will need to convert to the original observation
        # identifier.

        idinfo <- read.csv("/scratch/groups/smill191/smiller/data/OCO2_MIP/tccon/tccon_id_information.csv",header=F);
	colnames(idinfo) <- c("myid","tcconid")


#------------------------------------#
# Loop over each run type and year   #
#------------------------------------#

        # List of years to process
        yearlist <- 2015:2020

        # Names of the different MIP runs (contained in "inpath")
        runnames <- c("prior","IS","LNLG","LNLGIS","LNLGOGIS","OG")

        for(runtype in runnames){
        print("------------")
        print(runtype)
        print("------------")

        residall <- NULL

        for(yr in yearlist){
        print(yr)

        # Set the path to the fluxes
        if(yr==2015)         runpath <- paste(inpath,runtype,"/",sep="")
        if(yr >2015)         runpath <- paste(inpath,runtype,"_",yr,"/",sep="")


#------------------------#
# Concatinate residuals  #
#------------------------#

	# Read in obspack ID
	resid1 <- read.csv(paste(runpath,"tccon_id.csv",sep=""),header=F,strip.white=T)

	# Read in the model-data residuals
	resid2 <- read.csv(paste(runpath,"residuals_tccon.csv",sep=""),header=F)

	# Concatinate the two together
	resid1 <- resid1[,1:2]
	resid2 <- resid2[,1:4]
	resid <- cbind(resid1,resid2)

	print("Mean diff for given year")
	print(mean(resid[,5]-resid[,6]))
	
	print("RMSE of residuals")
	print(sqrt(mean( (resid[,5] - resid[,6])^2) ))

	print("SD of residuals")
	print(sqrt(sd(resid[,5] - resid[,6])))

	# Concatinate to existing residual file
	residall <- rbind(residall,resid)

	} # End of yr loop


#---------------------------------------------------------------#
# Replace my observation ID with the real TCCON observation ID  #
#---------------------------------------------------------------#

	# Set column names on the TCCON file
        colnames(residall) <- c("identifier","site","lat","lon","model","obs")

        # print("SCOT TEST")
	# print(residall[10000:10010,])

	# Match the residual file identifier to the IDs in the identifier file
	sel <- match(x=residall[,"identifier"],table=idinfo[,"myid"])	

	residall[,"identifier"] <- idinfo[sel,"tcconid"]

	# options(digits=16)
	# print(residall[10000:10010,])
	# print(idinfo[sel[10000:10010],])
        # print("------")


#----------------------------#
# Write the final csv file   #
#----------------------------#

	outname <- paste(outpath,"OCO2_v10mip.tccon.",runtype,".JHU.csv",sep="")

	print("Mean model")
	print(mean(residall[,"model"]))
	
	print("Mean obs")
	print(mean(residall[,"obs"]))

	print("RMSE")
	print(sqrt(mean( (residall[,"model"] - residall[,"obs"])^2) ))

	print("Standard deviation of residuals")
	print(sd(residall[,"model"] - residall[,"obs"]))

	# sel <- residall[,"site"]==9
	# print("Bias at a specific site")
	# print(mean( (residall[sel,"model"] - residall[sel,"obs"])))

	# print("RMSE at a specific site")
        # print(sqrt(mean( (residall[sel,"model"] - residall[sel,"obs"])^2) ))

	# print("SD at a specific site")
	# print(sd(residall[sel,"model"] - residall[sel,"obs"]))

	# Write to file	
	write.csv(residall,file=outname,row.names=F)

	} # End of runtype loop
	

#---------------------------------------------------------------------------------------------------------------#
# END OF SCRIPT
#---------------------------------------------------------------------------------------------------------------#

