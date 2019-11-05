

####  RealPlots is the code to create plots in real data application
#### Input
#### Argument		Type			            Description
#### mcmcinfer_object  	output of "inference.R"	      Results for "inference.R"
#### p			Vector of length P  in (0,1)	Probability grids
#### opt			1 or 0 	                  1 gives local significant regions while 0 does not give it.
#### edit	            scalar                        Adjust margin in legend
#### Output
#### Figure sumarizes estimation and inference for each covariate in the model.
#### For each covariate, there is one panel representing functional estimator b(p) and other panel contains density
#### estimates for different leves of each covariate


RealPlots = function( mcmcinfer_object ,p , edit=10, opt = 1 ){

par( mfrow=c(3,4), mar=c(2.5, 2.5, 3,2 ))
#################sex ###########################################
  plot(   0, type="n",    ylim= c( -20,  12), xlim=c(0,1)   , main="" )
 lines( p,  mcmcinfer_object$DataEst[,2], col="red"  , lty=1  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIu[,2], col="gray"  , lty=2  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIl[,2], col="gray"  , lty=2  , lwd=2 )  
 lines(  p , mcmcinfer_object$jointCI[ ,2,1 ], col="black" , lty=3  , lwd=2 )
 lines(  p , mcmcinfer_object$jointCI[ ,2,2 ], col="black" , lty=3  , lwd=2 )

mcmcinfer_object$local_p

title( paste0( paste0( "Sex (p=", round( mcmcinfer_object$global_p[1] ,3) ), sep="", ")", sep="" ), cex=1.5 )
points(p[mcmcinfer_object$local_p[,1]],  -20*mcmcinfer_object$local_p[,1][mcmcinfer_object$local_p[,1]==TRUE] , col="orange"   )
##points(p[flag2_R],  -20*flag2_R[flag2_R==TRUE] , col="orange"   )

if(opt==1){
legend( 0.1, 14,  pch=c(NA,NA,NA,19 ),
    c("Male effect", "95% CI (point)","95% CI (joint)","Flag") , 
	        lty= c(1,2,3,1)  ,
   col =c("red","gray","black" ,"orange") ,
   cex = 1.5 , bty = "n", ncol=1)
           }

if(opt==0){
legend( 0.1, 14,  pch=c(NA,NA,NA,NA ),
    c("Male effect", "95% CI (point)","95% CI (joint)",NA) , 
	        lty= c(1,2,3,NA)  ,
   col =c("red","gray","black" ,NA) ,
   cex = 1.5 , bty = "n", ncol=1)
           }



  plot(   0, type="n",    ylim=c(0,0.18), xlim=c(-10,60)  )
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,2] , col="black", lty=1 , lwd=1)
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,1] , col="red", lty=1 , lwd=2)  
title( "Male vs Female" , cex=1.5)

legend( -10+edit, 0.18,
    c("Male dist", "Female dist",
      paste0( paste0( "Mean shift (p=", round( mcmcinfer_object$mu_diff[1, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Var shift (p=", round( mcmcinfer_object$sigma_diff[1, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Skewed shift (p=", round(mcmcinfer_object$mu3_diff[1, 8] ,3) ), sep="", ")", sep="" )  ),
	        lty= c(1,1,NA,NA,NA)  ,
   col =c("red","black",NA ,NA, NA) ,
   cex = 1.5 , bty = "n", ncol=1)

####################age #######################################################

  plot(   0, type="n",    ylim= c( -25,  45), xlim=c(0,1)   , main="" )
 lines( p,  mcmcinfer_object$DataEst[,6], col="blue"  , lty=1  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIu[,6], col="gray"  , lty=2  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIl[,6], col="gray"  , lty=2  , lwd=2 )  
 lines(  p , mcmcinfer_object$jointCI[ ,6,1 ], col="black" , lty=3  , lwd=2 )
 lines(  p , mcmcinfer_object$jointCI[ ,6,2 ], col="black" , lty=3  , lwd=2 )


title( paste0( paste0( "Age (p=", round( mcmcinfer_object$global_p[5] ,3) ), sep="", ")", sep="" ) , cex=1.5)
points(p[mcmcinfer_object$local_p[,5]],  -20*mcmcinfer_object$local_p[,5][mcmcinfer_object$local_p[,5]==TRUE] , col="orange"   )

legend( 0.1, 47,  pch=c(NA,NA,NA,19 ),
    c("Age effect", "95% CI (point)","95% CI (joint)","") , 
	        lty= c(1,2,3,NA)  ,
   col =c("blue","gray","black" ,NA) ,
   cex = 1.5 , bty = "n", ncol=1)

  plot(   0, type="n",    ylim=c(0,0.18), xlim=c(-10,60)  )
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,4] , col="black", lty=1 , lwd=1)  
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,3] , col="blue", lty=1 , lwd=2)

title( "Age (84 years) vs Age (60 years)", cex=1.5 )
 ##  cbind(  AGE - median(AGE),  DATA$Age )

legend( -10+edit, 0.18,
    c("Age (84 years) dist", "Age (60 years) dist",
      paste0( paste0( "Mean shift (p=", round( mcmcinfer_object$mu_diff[2, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Var shift (p=", round( mcmcinfer_object$sigma_diff[2, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Skewed shift (p=", round(mcmcinfer_object$mu3_diff[2, 8] ,3) ), sep="", ")", sep="" )  ),
	        lty= c(1,1,NA,NA,NA)  ,
   col =c("blue","black",NA ,NA, NA) ,
   cex = 1.5 , bty = "n", ncol=1)

#################### gene 1 #######################################################

  plot(   0, type="n",    ylim= c( -15,  45), xlim=c(0,1)   , main="" )
 lines( p,  mcmcinfer_object$DataEst[,4], col="hotpink"  , lty=1  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIu[,4], col="gray"  , lty=2  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIl[,4], col="gray"  , lty=2  , lwd=2 )  
 lines(  p , mcmcinfer_object$jointCI[ ,4,1 ], col="black" , lty=3  , lwd=2 )
 lines(  p , mcmcinfer_object$jointCI[ ,4,2 ], col="black" , lty=3  , lwd=2 )

title( paste0( paste0( "DDIT3 (p=", round( mcmcinfer_object$global_p[3] ,3) ), sep="", ")", sep="" ) , cex=1.5)
points(p[mcmcinfer_object$local_p[,3]],  -15*mcmcinfer_object$local_p[,3][mcmcinfer_object$local_p[,3]==TRUE] , col="orange"   )

if(opt==1){
legend( 0.1, 47,  pch=c(NA,NA,NA,19 ),
    c("DDIT3 effect", "95% CI (point)","95% CI (joint)","Flag") , 
	        lty= c(1,2,3,1)  ,
   col =c("hotpink","gray","black" ,"orange") ,
   cex = 1.5 , bty = "n", ncol=1)
		}

if(opt==0){
legend( 0.1, 47,  pch=c(NA,NA,NA,NA ),
    c("DDIT3 effect", "95% CI (point)","95% CI (joint)",NA) , 
	        lty= c(1,2,3,NA)  ,
   col =c("hotpink","gray","black" ,  NA) ,
   cex = 1.5 , bty = "n", ncol=1)
		}



  plot(   0, type="n",    ylim=c(0,0.18), xlim=c(-10,60)  )
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,6] , col="black", lty=1 , lwd=1)  
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,5] , col="hotpink", lty=1 , lwd=2)

title( "DDIT3 vs None", cex=1.5 )


legend( -10+edit, 0.18,
    c("DDIT3 dist", "None dist",
      paste0( paste0( "Mean shift (p=", round(mcmcinfer_object$mu_diff[3, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Var shift (p=", round( mcmcinfer_object$sigma_diff[3, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Skewed shift (p=", round(mcmcinfer_object$mu3_diff[3, 8] ,3) ), sep="", ")", sep="" )  ),
	        lty= c(1,1,NA,NA,NA)  ,
   col =c("hotpink","black",NA ,NA, NA) ,
   cex = 1.5 , bty = "n", ncol=1)
########################gene 2###############################################

  plot(   0, type="n",    ylim= c( -15,  25), xlim=c(0,1)   , main="" )
 lines( p,  mcmcinfer_object$DataEst[,5], col="orange"  , lty=1  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIu[,5], col="gray"  , lty=2  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIl[,5], col="gray"  , lty=2  , lwd=2 )  
 lines(  p , mcmcinfer_object$jointCI[ ,5,1 ], col="black" , lty=3  , lwd=2 )
 lines(  p , mcmcinfer_object$jointCI[ ,5,2 ], col="black" , lty=3  , lwd=2 )

title( paste0( paste0( "EGFR (p=", round( mcmcinfer_object$global_p[4] ,3) ), sep="", ")", sep="" ) , cex=1.5)
points(p[mcmcinfer_object$local_p[,4]],  -25*mcmcinfer_object$local_p[,4][mcmcinfer_object$local_p[,4]==TRUE] , col="orange"   )

legend( 0.1, 27,  pch=c(NA,NA,NA,19 ),
    c("EGFR effect", "95% CI (point)","95% CI (joint)",NA) , 
	        lty= c(1,2,3,NA)  ,
   col =c("orange","gray","black" ,NA) ,
   cex = 1.5 , bty = "n", ncol=1)


  plot(   0, type="n",    ylim=c(0,0.18), xlim=c(-10,60)  )
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,8] , col="black", lty=1 , lwd=1)  
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,7] , col="orange", lty=1 , lwd=2)

title( "EGFR vs None" , cex=1.5)

legend( -10+edit, 0.18,
    c("EGFR dist", "None dist",
      paste0( paste0( "Mean shift (p=", round( mcmcinfer_object$mu_diff[4, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Var shift (p=", round(mcmcinfer_object$sigma_diff[4, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Skewed shift (p=", round(mcmcinfer_object$mu3_diff[4, 8] ,3) ), sep="", ")", sep="" )  ),
	        lty= c(1,1,NA,NA,NA)  ,
   col =c("orange","black",NA ,NA, NA) ,
   cex = 1.5 , bty = "n", ncol=1)
 
######################## tumor###############################################

  plot(   0, type="n",    ylim= c( -10,  35), xlim=c(0,1)   , main="" )
 lines( p,  mcmcinfer_object$DataEst[,3], col="green"  , lty=1  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIu[,3], col="gray"  , lty=2  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIl[,3], col="gray"  , lty=2  , lwd=2 )  
 lines(  p , mcmcinfer_object$jointCI[ ,3,1 ], col="black" , lty=3  , lwd=2 )
 lines(  p , mcmcinfer_object$jointCI[ ,3,2 ], col="black" , lty=3  , lwd=2 )

title( paste0( paste0( "Mesenchymal (p=", round( mcmcinfer_object$global_p[2] ,3) ), sep="", ")", sep="" ) , cex=1.5)
points(p[mcmcinfer_object$local_p[,2]],  -25*mcmcinfer_object$local_p[,2][mcmcinfer_object$local_p[,2]==TRUE] , col="orange"   )

legend( 0.1, 37,  pch=c(NA,NA,NA,19 ),
    c("Mesenchymal effect", "95% CI (point)","95% CI (joint)",NA) , 
	        lty= c(1,2,3,NA)  ,
   col =c("green","gray","black" ,NA) ,
   cex = 1.5 , bty = "n", ncol=1)


  plot(   0, type="n",    ylim=c(0,0.18), xlim=c(-10,60)  )
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,10] , col="black", lty=1 , lwd=1)  
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,9] , col="green", lty=1 , lwd=2)

title( "Mesenchymal vs None" , cex=1.5)

legend( -10+edit, 0.18,
    c("Mesenchymal dist", "None dist",
      paste0( paste0( "Mean shift (p=", round( mcmcinfer_object$mu_diff[5, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Var shift (p=", round( mcmcinfer_object$sigma_diff[5, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Skewed shift (p=", round( mcmcinfer_object$mu3_diff[5, 8] ,3) ), sep="", ")", sep="" )  ),
	        lty= c(1,1,NA,NA,NA)  ,
   col =c("green","black",NA ,NA, NA) ,
   cex = 1.5 , bty = "n", ncol=1)

######################## survival ###############################################

  plot(   0, type="n",    ylim= c( -20,  25), xlim=c(0,1)   , main="" )
 lines( p,  mcmcinfer_object$DataEst[,7], col="purple"  , lty=1  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIu[,7], col="gray"  , lty=2  , lwd=2 )  
 lines( p,  mcmcinfer_object$estCIl[,7], col="gray"  , lty=2  , lwd=2 )  
 lines(  p , mcmcinfer_object$jointCI[ ,7,1 ], col="black" , lty=3  , lwd=2 )
 lines(  p , mcmcinfer_object$jointCI[ ,7,2 ], col="black" , lty=3  , lwd=2 )

title( paste0( paste0( "Survival event (p=", round( mcmcinfer_object$global_p[6] ,3) ), sep="", ")", sep="" ), cex=1.5 )
points(p[mcmcinfer_object$local_p[,6]],  -20*mcmcinfer_object$local_p[,6][mcmcinfer_object$local_p[,6]==TRUE] , col="orange"   )

legend( 0.1, 27,  pch=c(NA,NA,NA,19 ),
    c("Survival event effect", "95% CI (point)","95% CI (joint)",NA) , 
	        lty= c(1,2,3,NA)  ,
   col =c("purple","gray","black" ,NA) ,
   cex = 1.5 , bty = "n", ncol=1)

  plot(   0, type="n",    ylim=c(0,0.18), xlim=c(-10,60)  )
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,12] , col="black", lty=1 , lwd=1) 
 lines(  xdomain[-1] ,  mcmcinfer_object$den_G[,11] , col="purple", lty=1 , lwd=2) 
title( "Event (< 12 months) vs Event (> 12 months)", cex=1.5 )

legend( -10+edit, 0.18,
    c("Event (< 12 months) dist", "Event (> 12 months) dist",
      paste0( paste0( "Mean shift (p=", round( mcmcinfer_object$mu_diff[6, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Var shift (p=", round( mcmcinfer_object$sigma_diff[6, 8] ,3) ), sep="", ")", sep="" ) ,
      paste0( paste0( "Skewed shift (p=", round(mcmcinfer_object$mu3_diff[6, 8] ,3) ), sep="", ")", sep="" )  ),
	        lty= c(1,1,NA,NA,NA)  ,
   col =c("purple","black",NA ,NA, NA) ,
   cex = 1.5 , bty = "n", ncol=1)

}
################################################################################
################################################################################
################################################################################

