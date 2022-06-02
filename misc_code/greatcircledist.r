#-------------------------------------------------------------------------------------------#
# FUNCTION: greatcircledist()                                                               #
# PURPOSE: Comute the great circle distance between two lat/lon coordinates. 	            #
# Written by Dan Matross                                                                    #
#   										     	    #
#-------------------------------------------------------------------------------------------#

	# Function: greatcircledist
	# usage: dist<-GreatCicleDist(lat1,lat2,lon1,lon2)
	# Purpose: Determines great cirle distance between two points on the suface of the earth
	# given a latitude and longitude for each.

	# FUNCTION INPUTS: lat/lon in degrees
	# FUNCTION OUTPUTS: distance in meters

	# EXPLANATION OF THE FUNCTION MATHEMATICS
	# Most easily considered as the arc-length of the circle between two points (that 
	# circle, of course has radius R.earth). We can use spherical polar coordinates
	# and the dot-product to determine this.

	# x=rcos(lat)cos(lon)
	# y=rcos(lat)sin(lon)
	# z=rsin(lat)

	# dot product p,q = length.p * length.q * cos(t)
	# where t is the angle between p and q

	# length(p)=length(q) = R

	# Take the dot product for (x1,y1,z1) (x2,y2,z2)
	# Rcos(lat1)cos(lon1)*R cos(lat2)cos(lon2) + Rcos(lat1)sin(lon1)*Rcos(lat2)sin(lat2) +
	# Rsin(lat1)*Rsin(lat2) = R*R*cos(t)

	# R^2s cancel
	# so t=arccos(cos(lat1)cos(lon1)cos(lat2)cos(lon2) + cos(lat1)sin(lon1)cos(lat2)sin(lat2) +sin(lat1)sin(lat2))

	# t represents the angle, distance is arc length
	# dist= t/2pi * 2pi * R = t*R

greatcircledist<-function(lat1,lat2,lon1,lon2)
{
	r.earth<-6.3710E6 #m

	theta1<-pi*lat1/180
	theta2<-pi*lat2/180
	lambda1<-pi*lon1/180
	lambda2<-pi*lon2/180

	trigout <- cos(theta1)*cos(lambda1)*cos(theta2)*cos(lambda2) +
                  cos(theta1)*sin(lambda1)*cos(theta2)*sin(lambda2) +
                  sin(theta1)*sin(theta2)
	trigout[trigout < -1] <- -1
	trigout[trigout >  1] <- 1

	ang<-acos(trigout)
	dist<-r.earth*ang

	return(dist)

} # End of the function
#---------------------------------------------------------------------------------
