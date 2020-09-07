library(deSolve)
library(FME)
library(ggplot2)
library(Cairo)
library(reshape2)
library(dplyr)
library(tidyverse)
library(latex2exp)

pars <- list(beta=0.19, sigma=0.5, gamma=0.2, Hu_Tot=300 )

solveCCHFV <- function(pars) {
  derivs <- function(t,state,pars) {    # returns rate of change
    with (as.list(c(state,pars)), {
      ## now code the model equations
      ######### Human ###############
      dHSdt <- -(beta*HS*HI)/(Hu_Tot)
      dHEdt <- (beta*HS*HI)/(Hu_Tot)-sigma*HE
      dHIdt <- sigma*HE-gamma*HI
      dHRdt <- gamma*HI
      
      return(list(c(dHSdt, dHEdt,dHIdt,dHRdt)))
    })
  }
  
  #state <- c(TS=800,TE=150,TI=50, LS=80, LE=15, LI=5, LR=0, HS=999, HE=1, HI=0, HR=0)
  #state <- c(TS=800,TE=150,TI=50, LS=70, LE=20, LI=10, LR=0, HS=975, HE=15, HI=5, HR=0)
  state <- c(HS=950, HE=40, HI=5, HR=5)
  tout  <- seq(0, 90, by = 1)
  ## ode solves the model by integration ...
  return(as.data.frame(ode(y = state, times = tout, func = derivs,
                           parms = pars)))
}

#c(TS=800,TE=150,TI=50, LS=70, LE=25, LI=5, LR=0,HS=950, HE=45, HI=5, HR=0)
out <- solveCCHFV(pars)
outframe<-data.frame(out)
plot(out$time, out$HI)

Inf_Data <- cbind.data.frame(out$time,out$HS,out$HE, out$HI, out$HR)

names(Inf_Data) <- c("Time","Susceptible", "Exposed","Infected", "Recovered")

library(ochRe)
# ochre_palettes$dead_reef
ggplot(Inf_Data)+
  geom_line(aes(x=Inf_Data$Time, y=Inf_Data$Susceptible),size=1,color="#758388") +
  geom_line(aes(x=Inf_Data$Time, y=Inf_Data$Exposed),size=1,color="#3C4347")+
  geom_line(aes(x=Inf_Data$Time, y=Inf_Data$Infected),size=1,color="#607848")+
  geom_line(aes(x=Inf_Data$Time, y=Inf_Data$Recovered),size=1,color="#548495")+
  theme_bw() + 
  theme(axis.title = element_text(color="black", face="bold", size=25))+
  labs(x = TeX('Time'), y= TeX('Population'), size = 80)+
  theme(axis.text = element_text(colour = "black", size = 25))+
  theme(axis.title = element_text( face="bold", size=40))+
  theme(axis.text = element_text(colour = "black", size = 25))+
  theme(axis.title = element_text( face="bold", size=25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(),axis.title.x=element_text(vjust=0,size=35),
        axis.title.y=element_text(vjust=0,size=35) )+
  theme(legend.position="top")
#theme(aspect.ratio=1)+
##########################################################################
library(dplyr)
library(reshape)
out1<-outframe[seq(from=0, to=91,by=1),]
tidy1<-melt(out1,id.vars = "time")
names(tidy1)<-c("Time","Population","Value")
ar1<-filter(tidy1,grepl("HS",Population))
ar2<-filter(tidy1,grepl("HE",Population))
ar3<-filter(tidy1,grepl("HI",Population))
ar4<-filter(tidy1,grepl("HR",Population))
Infect_Dat <- rbind.data.frame(ar1,ar2,ar3,ar4)

p1<-ggplot(Infect_Dat) +
  geom_line(aes(x=Time,y=Value,color=Population),size=2)+
  scale_colour_manual(values=c("#758388","#AB7E37","#d84860","#B9ACA3"))+
  theme_bw() + 
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=25))+
  scale_color_manual(labels = c(expression(paste("", Susceptible)), expression(paste("", Exposed)),
                                expression(paste("", Infected)),expression(paste("", Recovered))), 
                     values = c("#758388","#AB7E37","#d84860","#B9ACA3"))+
  theme(axis.title = element_text(color="black", face="bold", size=25))+
  labs(x = TeX('Time'), y= TeX('Population'), size = 80)+
  theme(axis.text = element_text(colour = "black", size = 25))+
  theme(axis.title = element_text( face="bold", size=40))+
  theme(axis.text = element_text(colour = "black", size = 25))+
  theme(axis.title = element_text( face="bold", size=25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_text(vjust=0,size=35),
        axis.title.y=element_text(vjust=0,size=35) )+
  theme(legend.position="bottom")
ggsave(filename = 'HIvsTIHu_TI_model.pdf',plot = last_plot(),width = 25, height = 25, units = "cm")
