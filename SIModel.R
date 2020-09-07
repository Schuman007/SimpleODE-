library(deSolve) # Solvers for Initial Value Problems of Differential Equations('ODE', 'DAE', 'DDE')
library(FME) # A Flexible Modelling Environment for Inverse Modelling, Sensitivity, Identifiability and Monte Carlo Analysis
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(Cairo) # R Graphics Device using Cairo Graphics Library for Creating High-Quality Bitmap (PNG, JPEG, TIFF), Vector (PDF, SVG,PostScript) and Display (X11 and Win32) Output
library(reshape2) # Flexibly Reshape Data: A Reboot of the Reshape Package
library(dplyr) # A Grammar of Data Manipulation
library(tidyverse) # Easily Install and Load the 'Tidyverse'
library(latex2exp) # Use LaTeX Expressions in Plots
library(ochRe) # Australia-Themed Color Palettes
library(reshape) # Flexibly Reshape Data
library("plot3D") # Three dimensional Plots 

pars <- list(beta=0.19, N=300 )

solveSIR <- function(pars) {
  derivs <- function(t,state,pars) {    # returns rate of change
    with (as.list(c(state,pars)), {
      ## now code the model equations
      dSdt <- -(beta*S*I)/N
      dIdt <- (beta*S*I)/N
      return(list(c(dSdt,dIdt)))
    })
  }
  
  state <- c(S=950, I=50)
  tout  <- seq(0, 90, by = 1)
  ## ode solves the model by integration ...
  return(as.data.frame(ode(y = state, times = tout, func = derivs,
                           parms = pars)))
}
out <- solveSIR(pars)
outframe<-data.frame(out)


out1<-outframe[seq(from=0, to=90,by=1),]
tidy1<-melt(out1,id.vars = "time")
names(tidy1)<-c("Time","Population","Value")
ar1<-filter(tidy1,grepl("S",Population))
ar2<-filter(tidy1,grepl("I",Population))
Infect_Dat <- rbind.data.frame(ar1,ar2)

p1<-ggplot(Infect_Dat) +
  geom_line(aes(x=Time,y=Value,color=Population),size=2)+
  scale_colour_manual(values=c("#758388","#AB7E37"))+
  theme_bw() + 
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=25))+
  scale_color_manual(labels = c(expression(paste("", Susceptible)),
                                expression(paste("", Infected))),
                     values = c("#758388","#AB7E37","#d84860"))+
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
ggsave(filename = 'SI.pdf',plot = last_plot(),width = 25, height = 25, units = "cm")
