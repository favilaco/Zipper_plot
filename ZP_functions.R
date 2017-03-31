##############################################################################################
##   Zipper plot: visualizing transcriptional activity of genomic regions				    ##
##                                                                                          ##
##   Author: Francisco Avila Cobos (Francisco.AvilaCobos@UGent.be)                          ##
##                                                                                          ##
##############################################################################################

##############################################################################################
##########################################  FANTOM5  #########################################
##############################################################################################

#Retrieves the closest CAGE peak on the same DNA strand for "CAGE_1tissue"
closest_CAGE <- function(User,CAGE){

	User <- User[,c(1:3),with=FALSE]
	setkey(User,TSS)

	w = CAGE[,start,]
	dt = CAGE[, val:=w]
	setkey(dt, val) 

	closest <- dt[User,roll="nearest",mult="first"] 

	if("+" %in% User$strand){

		closest[,`:=`(dist_start=start-val,dist_end=end-val)]

	} else {

		closest[,`:=`(dist_start=val-start,dist_end=val-end)]

	}

	return(closest[,list(dist_start,dist_end)])

}




#Retrieves the closest CAGE peak on the same DNA strand for "CAGE_ALL"
closest_CAGE_alltissues <- function(User,CAGE,detailed){

	User <- User[,c(1:3),with=FALSE]
	setkey(User,TSS)

	w=CAGE[,start,]
	dt = CAGE[, val:=w]
	setkey(dt, val) 
	closest <- dt[User,roll="nearest",mult="first"]

	if("+" %in% User$strand){

		closest[,`:=`(dist_start=start-val,dist_end=end-val)]

	} else {

		closest[,`:=`(dist_start=val-start,dist_end=val-end)]

	}

	if(detailed==TRUE){

		info <- closest[,4:(ncol(closest)-5),with=FALSE]
		closest$cases <- colnames(info)[apply(info, 1, which.max)]
		closest$signal <- apply(info, 1, function(x) max(x,na.rm=TRUE))
		return(closest[,list(chr,val,strand,dist_start,dist_end,signal,cases),])

	} else {

		return(closest[,list(dist_start,dist_end)])

	}

}




# Function employed to compute TSS p-val for "CAGE_1tissue"
closest_CAGE_signal <- function(User,CAGE){

	User <- User[,c(1:3),with=FALSE]
	setkey(User,TSS)

	w=CAGE[,start,] 
	dt = CAGE[, val:=w]
	setkey(dt, val) 
	closest <- dt[User,roll="nearest",mult="first"]

	if( "+" %in% User$strand){

		closest[,`:=`(dist_start=start-val,dist_end=end-val)]

	} else {

		closest[,`:=`(dist_start=val-start,dist_end=val-end)]

	}

	return(closest[,c(6,5,3,8,9,4),with=FALSE]) #Has to be with numbers because CNhs... ID changes each time
}




ORDERING_CAGE <- function(x,Plot=FALSE){

	colnames(x) <- c("strand","CAGE_start","CAGE_end")
	setkey(x,CAGE_start,CAGE_end)
	x[,':='(Min=pmin(CAGE_start,CAGE_end),Max=pmax(CAGE_start,CAGE_end))]

	max_width_right <- max(x$Max)
	max_width_left <- abs(min(x$Min))

	x_head <- x[(x[,CAGE_start]<=0 & x[,CAGE_end]>=0) | (x[,CAGE_start]>=0 & x[,CAGE_end]<=0),]
	setorder(x_head,Min)
	x_head[,c("Min","Max"):=NULL,]

	x_tail <- x[ (x[,CAGE_start]>0 & x[,CAGE_end]>0) | (x[,CAGE_start]<0 & x[,CAGE_end]<0) ,]
	x_tail[,Min:=pmin(abs(CAGE_start),abs(CAGE_end))]
	setorder(x_tail,Min)
	x_tail[,c("Min","Max"):=NULL,]

	x <- rbind(x_head,x_tail)

	Rel_zipper_height = round(nrow(x_head)/(nrow(x)),4)
	zipper_height = nrow(x_head)

	if(nrow(x_tail)>0){

		x_right <- x_tail[as.vector(x_tail[,2,with=FALSE]>0)]
		x_left <- x_tail[as.vector(x_tail[,2,with=FALSE]<0)]
		R <- nrow(x_head)+nrow(x_right)
		R_grid <- as.numeric(R)*as.numeric(max_width_right)
		L <- nrow(x_head)+nrow(x_left)
		L_grid <- as.numeric(L)*as.numeric(max_width_left)

		if(nrow(x_head)>0){
			
			ifelse(nrow(x_right)!=0, AUZ_right <- sum(x_right[,2,with=FALSE]/R_grid), AUZ_right <- 0)
			ifelse(nrow(x_left)!=0, AUZ_left <- sum(abs(x_left[,2,with=FALSE])/L_grid), AUZ_left <- 0)
			
		} else {

			ifelse(nrow(x_right)!=0, AUZ_right <- sum(x_right[,2,with=FALSE]/R_grid), AUZ_right <- NA)
			ifelse(nrow(x_left)!=0, AUZ_left <- sum(abs(x_left[,2,with=FALSE])/L_grid), AUZ_left <- NA)
			
		}

	} else { 

		R_grid <- nrow(x_head)*as.numeric(max_width_right)
		L_grid <- nrow(x_head)*as.numeric(max_width_left)
		AUZ_right <- 0
		AUZ_left <- 0
		
	}

		
	if(Plot){

		if(nrow(x)>1){
			#Visualization purposes: Modifying widths (to uniform) and ensuring CAGE peaks within the window
			width=Zipper_width/25 + 1
			x$CAGE_end <- x$CAGE_start + width

			x$CAGE_start[x$CAGE_start > Zipper_width] <- Zipper_width 
			x$CAGE_start[x$CAGE_start < -Zipper_width] <- -Zipper_width 

			x$CAGE_end[x$CAGE_end > Zipper_width] <- Zipper_width
			x$CAGE_end[x$CAGE_end < -Zipper_width] <- -Zipper_width
		}

		if(nrow(x_tail)>0){

			x_right$CAGE_start[x_right$CAGE_start > Zipper_width] <- Zipper_width
			x_left$CAGE_start[x_left$CAGE_start < -Zipper_width] <- -Zipper_width
			
			if(nrow(x_head)>0){
				
					ifelse(nrow(x_right)!=0, AUZr_window <- sum(x_right[,2,with=FALSE]/(R*Zipper_width)), AUZr_window <- 0)
					ifelse(nrow(x_left)!=0, AUZl_window <- sum(abs(x_left[,2,with=FALSE])/(L*Zipper_width)), AUZl_window <- 0)
					
				} else {

					ifelse(nrow(x_right)!=0, AUZr_window <- sum(x_right[,2,with=FALSE]/(R*Zipper_width)), AUZr_window <- NA)
					ifelse(nrow(x_left)!=0, AUZl_window <- sum(abs(x_left[,2,with=FALSE])/(L*Zipper_width)), AUZl_window <- NA)
					
				}

		} else {

			AUZr_window <- 0
			AUZl_window <- 0
		
		}
		
		return(list(x=x,Rel_zipper_height=Rel_zipper_height,zipper_height=zipper_height,
					R_grid=R_grid,L_grid=L_grid,AUZ_right=AUZ_right,AUZ_left=AUZ_left,
					AUZr_window=AUZr_window, AUZl_window=AUZl_window))
		

	} else {

		return(list(x=x,Rel_zipper_height=Rel_zipper_height,zipper_height=zipper_height,
				R_grid=R_grid,L_grid=L_grid,AUZ_right=AUZ_right,AUZ_left=AUZ_left))

		} 
		
}




CAGE_ZP <- function(Z,Zipper_width,COLOUR){

	if(nrow(Z$x)<100){
		PLOT <- ggplot(Z$x, aes(x=1:nrow(Z$x), ymin=CAGE_start,ymax=CAGE_end ))+
		geom_linerange(colour=COLOUR,size= (100 - nrow(Z$x) + 1)/nrow(Z$x)) + theme(legend.position="none") +
		coord_flip()+ 
		scale_x_reverse() +  
		theme_bw()+
		ylim(-Zipper_width,Zipper_width)+
		geom_hline(yintercept=0, linetype="dotted",size=0.3)+
		labs(title = "CAGE_FANTOM5_ZipperPlot",y="Distance_TSS_to_CAGE_peak (nt)",x="Num. Genomic Features") + theme( panel.grid= element_blank()) 

	} else {

		PLOT <- ggplot(Z$x, aes(x=1:nrow(Z$x), ymin=CAGE_start,ymax=CAGE_end ))+
		geom_linerange(colour=COLOUR) + theme(legend.position="none") +
		coord_flip()+ 
		scale_x_reverse() +  
		theme_bw()+
		ylim(-Zipper_width,Zipper_width)+
		geom_hline(yintercept=0, linetype="dotted",size=0.3)+
		labs(title = "CAGE_FANTOM5_ZipperPlot",y="Distance_TSS_to_CAGE_peak (nt)",x="Num. Genomic Features") + theme( panel.grid= element_blank()) 
	}

	return(PLOT)

}


##############################################################################################
##########################################  ROADMAP  #########################################
##############################################################################################

closest_ROADMAP_signal <- function(User,Roadmap,detailed){

	colnames(User) <- c("chr","TSS","strand")
	setkey(Roadmap,start)
	setkey(User,TSS)
	w=Roadmap[,start,]
	dt = Roadmap[, val:=w]
	setkey(dt, val) 
	closest <- dt[User,roll="nearest",mult="first"]

	if( "+" %in% User$strand){

		closest[,`:=`(dist_start=start-val,dist_end=end-val)]

	} else {

		closest[,`:=`(dist_start=val-start,dist_end=val-end)]

	}

	if(detailed==TRUE){

		info <- closest[,3:(ncol(closest)-5),with=FALSE]
		closest$cases <- colnames(info)[apply(info, 1, which.max)]
		closest$signal <- apply(info, 1, function(x) max(x,na.rm=TRUE))
		return(closest[,list(chr,val,strand,dist_start,dist_end,signal,cases),])

	} else {

		return(closest[,list(dist_start,dist_end)])

	}
}




ORDERING_ROADMAP <- function(x,Plot=FALSE){

	if(ncol(x)>3) x <- x[,list(strand,start,end),]
	
	colnames(x) <- c("strand","Roadmap_start","Roadmap_end")

	setkey(x,Roadmap_start,Roadmap_end)
	x[,':='(Min=pmin(Roadmap_start,Roadmap_end),Max=pmax(Roadmap_start,Roadmap_end))]

	max_width_right <- max(x$Max)
	max_width_left <- abs(min(x$Min))

	x_head <- x[(x[,Roadmap_start]<=0 & x[,Roadmap_end]>=0) | (x[,Roadmap_start]>=0 & x[,Roadmap_end]<=0),]
	setorder(x_head,Min)
	x_head[,c("Min","Max"):=NULL,]

	x_tail <- x[ (x[,Roadmap_start]>0 & x[,Roadmap_end]>0) | (x[,Roadmap_start]<0 & x[,Roadmap_end]<0) ,]
	x_tail[,Min:=pmin(abs(Roadmap_start),abs(Roadmap_end))]
	setorder(x_tail,Min)
	x_tail[,c("Min","Max"):=NULL,]

	x <- rbind(x_head,x_tail)

	Rel_zipper_height = round(nrow(x_head)/(nrow(x)),4)
	zipper_height = nrow(x_head)

	if(nrow(x_tail)>0){

		x_right <- x_tail[as.vector(x_tail[,2,with=FALSE]>0)]
		x_left <- x_tail[as.vector(x_tail[,2,with=FALSE]<0)]
		R <- nrow(x_head)+nrow(x_right)
		R_grid <- as.numeric(R)*as.numeric(max_width_right)
		L <- nrow(x_head)+nrow(x_left)
		L_grid <- as.numeric(L)*as.numeric(max_width_left)

		if(nrow(x_head)>0){
			
			ifelse(nrow(x_right)!=0, AUZ_right <- sum(x_right[,2,with=FALSE]/R_grid), AUZ_right <- 0)
			ifelse(nrow(x_left)!=0, AUZ_left <- sum(abs(x_left[,2,with=FALSE])/L_grid), AUZ_left <- 0)
			
		} else {

			ifelse(nrow(x_right)!=0, AUZ_right <- sum(x_right[,2,with=FALSE]/R_grid), AUZ_right <- NA)
			ifelse(nrow(x_left)!=0, AUZ_left <- sum(abs(x_left[,2,with=FALSE])/L_grid), AUZ_left <- NA)
			
		}

	} else {

		R_grid <- nrow(x_head)*as.numeric(max_width_right)
		L_grid <- nrow(x_head)*as.numeric(max_width_left)
		AUZ_right <- 0
		AUZ_left <- 0
		
	}

	if(Plot){

		x$Roadmap_start[x$Roadmap_start > Zipper_width] <- Zipper_width
		x$Roadmap_start[x$Roadmap_start < -Zipper_width] <- -Zipper_width

		x$Roadmap_end[x$Roadmap_end > Zipper_width] <- Zipper_width
		x$Roadmap_end[x$Roadmap_end < -Zipper_width] <- -Zipper_width

		if(nrow(x_tail)>0){

			# Calculating AUZ_window: Grid width = Zipper_width
			### Adjusting x_right & x_left to window [on unmodified data!]
			x_right$Roadmap_start[x_right$Roadmap_start > Zipper_width] <- Zipper_width
			x_left$Roadmap_start[x_left$Roadmap_start < -Zipper_width] <- -Zipper_width
			
			if(nrow(x_head)>0){
			
				ifelse(nrow(x_right)!=0, AUZr_window <- sum(x_right[,2,with=FALSE]/(R*Zipper_width)), AUZr_window <- 0)
				ifelse(nrow(x_left)!=0, AUZl_window <- sum(abs(x_left[,2,with=FALSE])/(L*Zipper_width)), AUZl_window <- 0)
				
			} else {

				ifelse(nrow(x_right)!=0, AUZr_window <- sum(x_right[,2,with=FALSE]/(R*Zipper_width)), AUZr_window <- NA)
				ifelse(nrow(x_left)!=0, AUZl_window <- sum(abs(x_left[,2,with=FALSE])/(L*Zipper_width)), AUZl_window <- NA)
				
			}

		} else {

			AUZr_window <- 0
			AUZl_window <- 0
		
		}

	return(list(x=x,Rel_zipper_height=Rel_zipper_height,zipper_height=zipper_height,
            	R_grid=R_grid,L_grid=L_grid,AUZ_right=AUZ_right,AUZ_left=AUZ_left,
				AUZr_window=AUZr_window, AUZl_window=AUZl_window))

	} else {

		return(list(x=x,Rel_zipper_height=Rel_zipper_height,zipper_height=zipper_height,
	            	R_grid=R_grid,L_grid=L_grid,AUZ_right=AUZ_right,AUZ_left=AUZ_left))

	}
}




ROADMAP_ZP <- function(Z,Zipper_width,COLOUR){

	if(nrow(Z$x)<100){

		PLOT <- ggplot(Z$x, aes(x=1:nrow(Z$x), ymin=Roadmap_start,ymax=Roadmap_end ))+
		geom_linerange(colour=COLOUR,size= (100 - nrow(Z$x) + 1)/nrow(Z$x)) + theme(legend.position="none") +
		coord_flip()+ 
		scale_x_reverse() +  
		theme_bw()+
		ylim(-Zipper_width,Zipper_width)+
		geom_hline(yintercept=0, linetype="dotted",size=0.3)+
		labs(title = paste(MARK,"ZipperPlot",sep="_"),y="Distance_TSS_to_hist_peak (nt)",x="Num. Genomic Features") + theme( panel.grid= element_blank()) 

	} else {

		PLOT <- ggplot(Z$x, aes(x=1:nrow(Z$x), ymin=Roadmap_start,ymax=Roadmap_end ))+
		geom_linerange(colour=COLOUR) + theme(legend.position="none") +
		coord_flip()+ 
		scale_x_reverse() +  
		theme_bw()+
		ylim(-Zipper_width,Zipper_width)+
		geom_hline(yintercept=0, linetype="dotted",size=0.3)+
		labs(title = paste(MARK,"ZipperPlot",sep="_"),y="Distance_TSS_to_hist_peak (nt)",x="Num. Genomic Features") + theme( panel.grid= element_blank()) 

	}

	return(PLOT)

}
