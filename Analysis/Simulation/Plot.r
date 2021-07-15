##Plot for simulation of transmitted variants in each blastocysts

raw<-read.csv("Permutation_10000_13_in_14.csv",header=T)
library(ggplot2)
ggplot(raw,aes(x=MAX))+
	geom_density()+
	geom_vline(xintercept=5,col="red",linetype=2)+
	labs(title="P=0.0199")+
	xlim(0,8)+
	ylim(0,2)+
	theme_bw()+
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	
ggplot(raw,aes(x=SD))+
	geom_density()+
	geom_vline(xintercept=1.49173547429654,col="red",linetype=2)+
	labs(title="P=0.0034")+
	xlim(0,2)+
	ylim(0,4)+
	theme_bw()+
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(raw,aes(x=X.2.1))+
	geom_density()+
	geom_vline(xintercept=2)+
	labs(x=">2",title="P=0.1696")+
	theme_bw()

ggplot(raw,aes(x=X.1.1))+
	geom_density()+
	geom_vline(xintercept=3)+
	labs(x=">1",title="P=0.8275")+
	theme_bw()


raw<-read.csv("Permutation_10000_4_in_6.csv",header=T)
library(ggplot2)
ggplot(raw,aes(x=MAX))+
	geom_density()+
	geom_vline(xintercept=3,col="red",linetype=2)+
	labs(title="P=0.0954")+
	xlim(0,5)+
	ylim(0,3)+
	theme_bw()+
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	
ggplot(raw,aes(x=SD))+
	geom_density()+
	geom_vline(xintercept=1.211060141639,col="red",linetype=2)+
	labs(title="P=0.0954")+
	xlim(0,2)+
	ylim(0,8)+
	theme_bw()+
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


raw<-read.csv("Permutation_10000_2_in_10.csv",header=T)
library(ggplot2)
ggplot(raw,aes(x=MAX))+
	geom_density()+
	geom_vline(xintercept=1,col="red",linetype=2)+
	labs(title="P=1")+
	xlim(0,3)+
	ylim(0,10)+
	theme_bw()+
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	
ggplot(raw,aes(x=SD))+
	geom_density()+
	geom_vline(xintercept=0.421637021355784,col="red",linetype=2)+
	labs(title="P=1")+
	xlim(0,0.8)+
	ylim(0,50)+
	theme_bw()+
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
