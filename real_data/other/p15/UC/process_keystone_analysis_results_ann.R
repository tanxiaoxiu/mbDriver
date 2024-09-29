# Vanni Bucci, Ph.D.
# Assistant Professor
# Department of Biology
# Room: 335A
# University of Massachusetts Dartmouth
# 285 Old Westport Road
# N. Dartmouth, MA 02747-2300
# Phone: (508)999-9219
# Email: vbucci@umassd.edu
# Web: www.vannibucci.org
#-------------------------------------------------------------------------
rm(list=ls())
# function to check for installed packages
pkgTest <- function(x){
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE,repos="http://cran.rstudio.com/")
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

# function to check for installed packages using source
pkgTest_source <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE,type="source")
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

# test for packages. If present load them
pkgTest("ggplot2")
pkgTest("RColorBrewer")
pkgTest("ggthemes")
pkgTest("gridExtra")
pkgTest("scales")
pkgTest("doBy")
pkgTest("reshape2")
pkgTest_source("viridis")

args <- commandArgs(trailingOnly=TRUE)
print(args)
my_input_file <- args[1]
my_output_file <-args[2]


setwd("~/mbDriver/real_data/other/p15/UC/UC")
myd<-read.csv("output_diet/BVS.results.keystone_analysis.txt",sep="\t")

uPerturbationID <- unique(myd$PerturbationID)

for (ip in seq(1,length(uPerturbationID))){
  myd.sub.pert<-subset(myd,PerturbationID==uPerturbationID[ip])
  uParentID<-unique(myd.sub.pert$ParentID)
  my.list.2.g<-list()
  my.colours<-viridis
  my.colours.2<-colorRampPalette(brewer.pal(7, "Reds"), space="Lab")
  for (j in seq(1,length(uParentID))){
    my.list<-list()
    parentState<-myd.sub.pert[which(myd.sub.pert$StateID==uParentID[j]),]
    childStates<-myd.sub.pert[which(myd.sub.pert$ParentID==uParentID[j] & myd.sub.pert$StateID!=uParentID[j]),]
    parentState.densities<-parentState[,5:ncol(parentState)]
    childStates.densities<-childStates[,5:ncol(parentState)]
    cnt<-1
    my.list[[cnt]]<-NA
    for (i in seq(1,nrow(childStates.densities))){
      tr<-which(parentState.densities!=0 & childStates.densities[i,]==0)
      p<-parentState.densities 
      c<-childStates.densities[i,]
      p<-p[-tr]
      c<-c[-tr]
      error=p-c 
      rmse=sqrt(mean(as.matrix(error)^2))
      cnt<-cnt+1
      my.list[[cnt]]<-rmse
    }
    RMSE<-do.call('rbind',my.list)
    sorted.RMSE.childStates<-
      sort(RMSE[2:nrow(RMSE)],index.return=TRUE)
    ordered.states<-childStates$StateID[sorted.RMSE.childStates$ix]
    myd.out<-rbind(parentState,childStates)
    myd.out<-cbind(myd.out,RMSE)
    sorted.strains<-sort(colMeans(myd.out[,7:ncol(myd.out)-1]),index.return=TRUE) 
    myd.out.m<-melt(myd.out,id.vars="StateID")
    myd.out.m$StateID<-factor(myd.out.m$StateID,levels=rev(append(rev(ordered.states),parentState$StateID,after=0)))
    myd.out.m.strains<-myd.out.m[which(myd.out.m$variable!="ParentID" & 
                                         myd.out.m$variable!="N_species" &
                                         myd.out.m$variable!="RMSE" &
                                         myd.out.m$variable!="frequency" &
                                         myd.out.m$variable!="PerturbationID"),] 
    myd.out.m.strains$variable<-factor(myd.out.m.strains$variable,levels=names(myd.out)[sorted.strains$ix+5]) 
    myd.out.m.strains$value[myd.out.m.strains$value==0]<-NA 
    g1<-ggplot(myd.out.m.strains,aes(x=variable,y=StateID,
                                     fill=log10(value+1)))+
      geom_tile()+
      geom_tile(color="black")+
      scale_fill_gradientn(name ="log10", 
                           colours = my.colours(11), na.value = "transparent",space="Lab")+
      theme_few()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
      xlab("")+
      ylab("")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(panel.background=element_blank(),panel.border=element_blank())+
      theme(axis.line=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank())
    
    myd.out.m.rmse<-myd.out.m[which(myd.out.m$variable=="RMSE"),]
    g2<-ggplot(myd.out.m.rmse,aes(x=variable,y=StateID,
                                  fill=log10(value+1)))+
      geom_tile()+ 
      geom_tile(color="black")+
      scale_fill_gradientn(name ="log10", 
                           colours = my.colours.2(11), na.value = "transparent",space="Lab")+
      theme_few()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
      xlab("")+
      ylab("")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(panel.background=element_blank(),panel.border=element_blank())+
      theme(axis.line=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank())
    
    # put all together
    gA <- ggplotGrob(g1)
    gB <- ggplotGrob(g2)
    maxh = grid::unit.pmax(gA$heights[2:5], gB$heights[2:5])
    gA$heights[2:5] <- as.list(maxh)
    gB$heights[2:5] <- as.list(maxh)
    pdf(sprintf("%s.pert.%s.set.%d.pdf",my_output_file,ip,j),width=8, height=6)
    grid.arrange(gA, gB,  widths=c(2.0,1.0),ncol=2)
    dev.off()
    myd.out.m.rmse_name <- paste("myd.out.m.rmse",ip,j,".txt",sep = "_")
    write.table(myd.out.m.rmse,file = myd.out.m.rmse_name,row.names = F,col.names = T, sep = "\t",quote = F)
    
    ############Keystoneness
    library(dplyr)
    myd_out <- myd.out %>% select(-RMSE)
    StateID_microbe <- data.frame(StateID = integer(), microbe = character(), stringsAsFactors = FALSE)
    state_ids <- unique(myd_out$StateID)
    all_zero_columns <- names(myd_out[, -which(names(myd_out) == "StateID")])[colSums(myd_out[, -which(names(myd_out) == "StateID")]) == 0]
    first_state_id <- state_ids[1]  
    for (microbe in all_zero_columns) {
      StateID_microbe <- rbind(StateID_microbe, data.frame(StateID = first_state_id, microbe = microbe))
    }
    
    for (state_id in state_ids) {
      row_data <- myd_out[myd_out$StateID == state_id, -which(names(myd_out) %in% c("StateID", all_zero_columns))]
      zero_microbes <- names(row_data)[row_data == 0]
      for (microbe in zero_microbes) {
        StateID_microbe <- rbind(StateID_microbe, data.frame(StateID = state_id, microbe = microbe))
      }
    }
    
    #######
    myd.out.m.rmse$StateID <- as.integer(as.character(myd.out.m.rmse$StateID))
    microbe_Keystoneness <- left_join(StateID_microbe, myd.out.m.rmse, by = "StateID")
    microbe_Keystoneness <- microbe_Keystoneness[, -3]
    names(microbe_Keystoneness)[names(microbe_Keystoneness) == "value"] <- "keystoneness"
    microbe_Keystoneness <- arrange(microbe_Keystoneness, desc(keystoneness))
    microbe_Keystoneness_name <- paste("microbe_Keystoneness",ip,j,".txt",sep = "_")
    write.table(microbe_Keystoneness,file = microbe_Keystoneness_name,row.names = F,col.names = T, sep = "\t",quote = F)
  }
}









          
          
