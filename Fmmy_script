#'Setting of runs ==============================================================================
#' 1) Read in data 
  #Japan_dat <- read.csv("Japan_data_1sp.csv") 
  Japan_dat <- readr::read_csv(system.file("extdata","Japan_data_1sp.csv",package="calcmmy"))
  #Japan_dat <- read.csv("Japan_data_extend.csv") 
		   
#' 2) Specify which quantile to use for steepness (h)
  set_hsd <- 1 # if 1 then sd, if 2 then 2sd		   

#' 3) Specify scenario
  set_analysis <- "species" # either "clark" | "group" | "species"
  
#' 4) for scenario "group" and "species", specify to use only BH or both (BH and RI). if 1= only BH, 2= both
  set_group_sr <- 1 # default is 1

#' 5) set name of output file for Tables
  output_folder <- "output_mmy"
if(!file.exists(output_folder)) dir.create(output_folder, recursive=TRUE)

#'End of settings ==============================================================================
#' Below is code (no need to change) 

#1: Read in stock name
Stock_name <- Japan_dat %>% distinct(Stock,.keep_all=TRUE)
stocks <- Stock_name$Stock

#2: Set range for steepness h
if(set_hsd==2){
set_low_quantile <- 0.023 #平均-2sdに対応  (-2sdなら0.023)(sdなら0.16)
set_high_quantile <- 0.977 #平均+2sdに対応 (+2sdなら0.977)(sdなら0.84)
} else {
set_low_quantile <- 0.16 #平均-1sdに対応  (-2sdなら0.023)(sdなら0.16)
set_high_quantile <- 0.84 #平均+1sdに対応 (+2sdなら0.977)(sdなら0.84)
}

 #3: Main code ==============================================================================

if(set_analysis == "clark"){
   ##1. Clark scenario ==============================================================================
   ##SR.coef = c(0.750,0.875,0.938,1.386,2.079,2.773),# coefficient for S-R relation
   ##is.ricker = c(F,F,F,T,T,T)
    pars <- c()
    for(i in 1:length(stocks)){
      temp<-Japan_dat %>% filter(Stock==stocks[[i]])
      bpara <- data.frame(waa=temp$waa,
                      maa=temp$maa,
                      M=temp$M,
                      F=temp$F,
                      Age=temp$Age,
                      Stock=temp$Stock)
  
      bpara$waa.mid <- bpara$waa
      bpara$maa.mid <- bpara$maa
  
      temp2 <- get.clark.data2(bpara,steepness=c(0.5,0.7,0.8,0.5,0.7,0.8), is.graph=T)
      pars <- rbind(pars,temp2)
    }
	
	## out put Table
    result_clark_h <- as.data.frame(pars, col.names = NULL)
    result_clark_h <- round_df(result_clark_h, 2)
    rownames(result_clark_h) <- stocks
    result_clark_h <- cbind(result_clark_h,as.data.frame(stocks))
	colnames(result_clark_h) <- c("Fmax","SPRmmy","Fmmy","Fmsy_BH05","Fmsy_BH07","Fmsy_BH08","Fmsy_RI05","Fmsy_RI07","Fmsy_RI08","SPRmsy_BH05","SPRmsy_BH07","SPRmsy_BH08","SPRmsy_RI05","SPRmsy_RI07","SPRmsy_RI08","Ymmy_Ymsy_BH05","Ymmy_Ymsy_BH07","Ymmy_Ymsy_BH08","Ymmy_Ymsy_RI05","Ymmy_Ymsy_RI07","Ymmy_Ymsy_RI08","Y35SPR_Ymsy_BH05","Y35SPR_Ymsy_BH07","Y35SPR_Ymsy_BH08","Y35SPR_Ymsy_RI05","Y35SPR_Ymsy_RI07","Y35SPR_Ymsy_RI08","Y40SPR_Ymsy_BH05","Y40SPR_Ymsy_BH07","Y40SPR_Ymsy_BH08","Y40SPR_Ymsy_RI05","Y40SPR_Ymsy_RI07","Y40SPR_Ymsy_RI08","h_low","h_med","h_hi","h_low2","h_med2","h_hi2","stock") 
	readr::write_csv(result_clark_h,path=str_c(output_folder,"/result_clark.csv"))
	##End of Clark scenario ==============================================================================
} 

if(set_analysis == "group"){
    ##2. Use h of family from FishLife ==============================================================================
      # Install and load package
      # devtools::install_github("james-thorson/FishLife")
    library( FishLife )
    vignette("tutorial","FishLife")
    
	# get h pars
    h_pars <-sd_pars_low <-sd_pars_high <-c()
  
    for(i in 1 : dim(Japan_dat)[[1]]){
      if(Japan_dat$Order[[i]]=="Perciformes"){
         Predict = Plot_taxa( Search_species(Family=Japan_dat$Family[[i]])$match_taxonomy, mfrow=c(3,2) )
        } else {
	     Predict = Plot_taxa( Search_species(Order=Japan_dat$Order[[i]])$match_taxonomy, mfrow=c(3,2) ) 
	    }
 
      temp <- Predict[[1]]$Mean_pred[[14]]
      temp_low <- qnorm(set_low_quantile,mean=Predict[[1]]$Mean_pred[[14]],sd=sqrt(Predict[[1]]$Cov_pred[14,14]))
      temp_high <- qnorm(set_high_quantile,mean=Predict[[1]]$Mean_pred[[14]],sd=sqrt(Predict[[1]]$Cov_pred[14,14]))
  
      h_pars <- rbind(h_pars,temp)
      sd_pars_low <- rbind(sd_pars_low,temp_low)
      sd_pars_high <- rbind(sd_pars_high,temp_high)
    }

    h_pars <- as.data.frame(h_pars) 
    sd_pars_low <- as.data.frame(sd_pars_low)
    sd_pars_high <- as.data.frame(sd_pars_high)

    h_pars <- 0.2 + 0.8*plogis(h_pars$V1) # back transform from logit
    hs_low <- 0.2 + 0.8*plogis(sd_pars_low$V1) # back transform from logit
    hs_high <- 0.2 + 0.8*plogis(sd_pars_high$V1) # back transform from logit

    hs<-h_pars

    Abs <- (hs-0.2)/(0.8*hs) # translate to BH shape parameter
    Abs_low <- (hs_low-0.2)/(0.8*hs_low) # translate to BH shape parameter
    Abs_high <- (hs_high-0.2)/(0.8*hs_high) # translate to BH shape parameter

    Ars <- log(-4*hs/(hs-1)) # translate to RI shape parameter (not used)
    Ars_low <- log(-4*hs_low/(hs_low-1))# translate to RI shape parameter (not used)
    Ars_high <- log(-4*hs_high/(hs_high-1))# translate to RI shape parameter (not used)
    #Ars_high[is.na(Ars_high)]<-5.981414

    Japan_dat_hs_Abs <- cbind(Japan_dat,hs_low,hs,hs_high,Abs_low,Abs,Abs_high,Ars_low,Ars,Ars_high)
    
	# apply h pars
    Stock_name <- Japan_dat_hs_Abs %>% distinct(Stock,.keep_all=TRUE)
    stocks <- Stock_name$Stock
    
    pars <- pars_h <- c()
    for(i in 1:length(stocks)){
      temp <- Japan_dat_hs_Abs %>% filter(Stock==stocks[[i]])
  
      bpara <- data.frame(waa=temp$waa,
                      maa=temp$maa,
                      M=temp$M,
                      F=temp$F,
                      Age=temp$Age,
                      Stock=temp$Stock)
  
      bpara$waa.mid <- bpara$waa
      bpara$maa.mid <- bpara$maa
  
     if(set_group_sr == 2){
         temp2 <- get.clark.data2(bpara,SR.coef=c(temp$Abs_low[1],temp$Abs[1],temp$Abs_high[1],1.386,2.079,2.773),is.ricker=c(F,F,F,T,T,T),steepness=c(temp$hs_low[1],temp$hs[1],temp$hs_high[1],0.5,0.7,0.8))
         pars_h <- rbind(pars_h,temp2) 
	 }
     if(set_group_sr == 1){
         temp2 <- get.clark.data2(bpara,SR.coef=c(temp$Abs_low[1],temp$Abs[1],temp$Abs_high[1]),is.ricker=c(F,F,F),steepness=c(temp$hs_low[1],temp$hs[1],temp$hs_high[1]))
         pars_h <- rbind(pars_h,temp2)
	 }
    }
    
	## out put Table
    result_group_h <- as.data.frame(pars_h, col.names = NULL)
	result_group_h <- round_df(result_group_h, 2)
    rownames(result_group_h) <- stocks
    result_group_h <- cbind(result_group_h,as.data.frame(stocks))
	
  if(set_group_sr == 1){
	 colnames(result_group_h) <- c("Fmax","SPRmmy","Fmmy","Fmsy_BH05","Fmsy_BH07","Fmsy_BH08","SPRmsy_BH05","SPRmsy_BH07","SPRmsy_BH08","Ymmy_Ymsy_BH05","Ymmy_Ymsy_BH07","Ymmy_Ymsy_BH08","Y35SPR_Ymsy_BH05","Y35SPR_Ymsy_BH07","Y35SPR_Ymsy_BH08","Y40SPR_Ymsy_BH05","Y40SPR_Ymsy_BH07","Y40SPR_Ymsy_BH08","h_low","h_med","h_hi","stock")
	 readr::write_csv(result_group_h,path=str_c(output_folder,"/result_group.csv"))
	}
  if(set_group_sr == 2){
	 colnames(result_group_h)<-c("Fmax","SPRmmy","Fmmy","Fmsy_BH05","Fmsy_BH07","Fmsy_BH08","Fmsy_RI05","Fmsy_RI07","Fmsy_RI08","SPRmsy_BH05","SPRmsy_BH07","SPRmsy_BH08","SPRmsy_RI05","SPRmsy_RI07","SPRmsy_RI08","Ymmy_Ymsy_BH05","Ymmy_Ymsy_BH07","Ymmy_Ymsy_BH08","Ymmy_Ymsy_RI05","Ymmy_Ymsy_RI07","Ymmy_Ymsy_RI08","Y35SPR_Ymsy_BH05","Y35SPR_Ymsy_BH07","Y35SPR_Ymsy_BH08","Y35SPR_Ymsy_RI05","Y35SPR_Ymsy_RI07","Y35SPR_Ymsy_RI08","Y40SPR_Ymsy_BH05","Y40SPR_Ymsy_BH07","Y40SPR_Ymsy_BH08","Y40SPR_Ymsy_RI05","Y40SPR_Ymsy_RI07","Y40SPR_Ymsy_RI08","h_low","h_med","h_hi","h_low2","h_med2","h_hi2","stock")
	 readr::write_csv(result_group_h,path=str_c(output_folder,"/result_group_BHRI.csv"))
	}	
}
## End of FishLife family scenario ==============================================================================

if(set_analysis == "species"){
 ##3. Use h of species from FishLife  ==============================================================================
  # Install and load package
  #devtools::install_github("james-thorson/FishLife")
  library( FishLife )
  vignette("tutorial","FishLife")
 
  # get h pars
  h_pars <- c()
  sd_pars_low <- c()
  sd_pars_high <- c ()
 
  for(i in 1 : dim(Japan_dat)[[1]]){
     Predict = Plot_taxa( Search_species(Genus=Japan_dat$Genus[[i]],Species=Japan_dat$Species[[i]])$match_taxonomy, mfrow=c(3,2) )
     temp <- Predict[[1]]$Mean_pred[[14]]
     temp_low <- qnorm(set_low_quantile,mean=Predict[[1]]$Mean_pred[[14]],sd=sqrt(Predict[[1]]$Cov_pred[14,14]))
     temp_high <- qnorm(set_high_quantile,mean=Predict[[1]]$Mean_pred[[14]],sd=sqrt(Predict[[1]]$Cov_pred[14,14]))
    
     h_pars <- rbind(h_pars,temp)
     sd_pars_low <- rbind(sd_pars_low,temp_low)
     sd_pars_high <- rbind(sd_pars_high,temp_high)
    }
  
  h_pars <- as.data.frame(h_pars) 
  sd_pars_low <- as.data.frame(sd_pars_low)
  sd_pars_high <- as.data.frame(sd_pars_high)
  
  h_pars <- 0.2 + 0.8*plogis(h_pars$V1)
  hs_low <- 0.2 + 0.8*plogis(sd_pars_low$V1)
  hs_high <- 0.2 + 0.8*plogis(sd_pars_high$V1)
  
  hs <- h_pars
  
  Abs <- (hs-0.2)/(0.8*hs)
  Abs_low <- (hs_low-0.2)/(0.8*hs_low)
  Abs_high <- (hs_high-0.2)/(0.8*hs_high)
  
  Ars <- log(-4*hs/(hs-1))
  Ars_low <- log(-4*hs_low/(hs_low-1))
  Ars_high <- log(-4*hs_high/(hs_high-1))
  
  Japan_dat_hs_Abs <- cbind(Japan_dat,hs_low,hs,hs_high,Abs_low,Abs,Abs_high,Ars_low,Ars,Ars_high)
  
  Stock_name <- Japan_dat_hs_Abs %>% distinct(Stock,.keep_all=TRUE)
  stocks <- Stock_name$Stock
  
   pars <- pars_h <- c()
  
  for(i in 1:length(stocks)){
    temp<-Japan_dat_hs_Abs %>% filter(Stock==stocks[[i]])
    
    bpara <- data.frame(waa=temp$waa,
                        maa=temp$maa,
                        M=temp$M,
                        F=temp$F,
                        Age=temp$Age,
                        Stock=temp$Stock)
    
    bpara$waa.mid <- bpara$waa
    bpara$maa.mid <- bpara$maa
    
	if(set_group_sr == 2){
         temp2<-get.clark.data2(bpara,SR.coef=c(temp$Abs_low[1],temp$Abs[1],temp$Abs_high[1],1.386,2.079,2.773),is.ricker=c(F,F,F,T,T,T),steepness=c(temp$hs_low[1],temp$hs[1],temp$hs_high[1],0.5,0.7,0.8))
         pars_h<- rbind(pars_h,temp2)
	}
	
	if(set_group_sr == 1){
        temp2<-get.clark.data2(bpara,SR.coef=c(temp$Abs_low[1],temp$Abs[1],temp$Abs_high[1]),is.ricker=c(F,F,F),steepness=c(temp$hs_low[1],temp$hs[1],temp$hs_high[1]))
        pars_h<- rbind(pars_h,temp2)	
	}
  }

  ## out put Table
    result_species_h<-as.data.frame(pars_h, col.names = NULL)
    result_species_h<-round_df(result_species_h, 2)
    rownames(result_species_h)<-stocks
    result_species_h<-cbind(result_species_h,as.data.frame(stocks))
	
	if(set_group_sr == 1){
	 colnames(result_species_h)<-c("Fmax","SPRmmy","Fmmy","Fmsy_BH05","Fmsy_BH07","Fmsy_BH08","SPRmsy_BH05","SPRmsy_BH07","SPRmsy_BH08","Ymmy_Ymsy_BH05","Ymmy_Ymsy_BH07","Ymmy_Ymsy_BH08","Y35SPR_Ymsy_BH05","Y35SPR_Ymsy_BH07","Y35SPR_Ymsy_BH08","Y40SPR_Ymsy_BH05","Y40SPR_Ymsy_BH07","Y40SPR_Ymsy_BH08","h_low","h_med","h_hi","stock")
	 readr::write_csv(result_species_h,path=str_c(output_folder,"/result_species.csv"))
	}
	if(set_group_sr == 2){
	 colnames(result_species_h)<-c("Fmax","SPRmmy","Fmmy","Fmsy_BH05","Fmsy_BH07","Fmsy_BH08","Fmsy_RI05","Fmsy_RI07","Fmsy_RI08","SPRmsy_BH05","SPRmsy_BH07","SPRmsy_BH08","SPRmsy_RI05","SPRmsy_RI07","SPRmsy_RI08","Ymmy_Ymsy_BH05","Ymmy_Ymsy_BH07","Ymmy_Ymsy_BH08","Ymmy_Ymsy_RI05","Ymmy_Ymsy_RI07","Ymmy_Ymsy_RI08","Y35SPR_Ymsy_BH05","Y35SPR_Ymsy_BH07","Y35SPR_Ymsy_BH08","Y35SPR_Ymsy_RI05","Y35SPR_Ymsy_RI07","Y35SPR_Ymsy_RI08","Y40SPR_Ymsy_BH05","Y40SPR_Ymsy_BH07","Y40SPR_Ymsy_BH08","Y40SPR_Ymsy_RI05","Y40SPR_Ymsy_RI07","Y40SPR_Ymsy_RI08","h_low","h_med","h_hi","h_low2","h_med2","h_hi2","stock")
     readr::write_csv(result_species_h,path=str_c(output_folder,"/result_species_BHRI.csv"))
	}
 }
##End of FishLife species scenario  ==============================================================================