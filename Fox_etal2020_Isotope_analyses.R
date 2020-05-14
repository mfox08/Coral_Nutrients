######## Fox et al. 2020 - Differential resistance and acclimation of two coral species to chronic nutrient
########                   enrichment reflect life history traits 

####### Code file 4 of 4 -- Stable isotope analyses for Figure 4
####### Required data files: 1) Isotopes_raw.csv
#                            2) Isotope nubbins.csv

####### Summary: This file completes all of the data summarization, analysis, and plotting for main text
#######          figure 4 (d15N, d13C, and C:N data)

####### Figures produced: Fig. 4

library(dplyr)
library(ggplot2)
library(cowplot)
require(lmerTest)
require(RColorBrewer)
require(emmeans) # for tukey tests of lmer model 
library(car)

#create plotting themes
newtheme <- theme_classic() + theme(text = element_text(size=11))+
  theme(axis.text.x = element_text(size=12,colour="black"), axis.text.y = element_text(size=12,colour="black"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(axis.line = element_line(size = 0))



nubs<-read.csv(file="Isotope nubbins.csv",header=T)
topes<-read.csv(file="Isotopes_raw.csv",header=T)

#join the data based on the sample ID 
chain<- left_join(nubs, topes)

# take means and SE
means<-chain %>% na.omit(chain) %>% #remove samples with NAs
  group_by(Species,Nutrient.Level,Tissue) %>% 
  summarize(
    N=sum(!is.na(d15N)),
    mean13C = mean(d13C,na.rm=T),
    mean15N = mean(d15N,na.rm=T),
    meanCN = mean(C.N,na.rm=T),
    se13C = sd(d13C)/sqrt(N),
    se15N = sd(d15N)/sqrt(N),
    seCN = sd(C.N)/sqrt(N))

#mean plots
cols<-c("navajowhite3","seagreen3") #set colors
labels<-c(POC="Pocillopora acuta",POR="Porites compressa") #set labels
means$Species<-ordered(means$Species,levels=c("POR","POC"))


# d15N
a<-ggplot(means,aes(x=Nutrient.Level,y=mean15N,ymin=mean15N-se15N,ymax=mean15N+se15N,fill=Tissue,group=Tissue,facets=Species))+
  geom_errorbar(position=position_dodge(0.9),width=0.0)+facet_wrap(~Species)+facet_grid(.~Species,labeller=labeller(Species=labels))+
  geom_point(aes(fill=Tissue),color="black",pch=21,size=2.5,position=position_dodge(0.9))+
  ylab(expression({delta}^15*N~'\u2030'))+
  scale_fill_manual(values=cols)+newtheme
a

#d13C
b<-ggplot(means,aes(x=Nutrient.Level,y=mean13C,ymin=mean13C-se13C,ymax=mean13C+se13C,group=Tissue,fill=Tissue,facets=Species))+
  geom_errorbar(position=position_dodge(0.9),width=0.0)+facet_wrap(~Species)+facet_grid(.~Species,labeller=labeller(Species=labels))+
  geom_point(aes(fill=Tissue),color="black",pch=21,size=2.5,position=position_dodge(0.9))+
  ylab(expression({delta}^13*C~'\u2030'))+
  scale_y_continuous(breaks=seq(-19,-16,by=1))+expand_limits(y = c(-19, -16))+
  scale_fill_manual(values=cols)+newtheme
b

# CN
c<-ggplot(means,aes(x=Nutrient.Level,y=meanCN,ymin=meanCN-seCN,ymax=meanCN+seCN,group=Tissue,fill=Tissue,facets=Species))+
  geom_errorbar(position=position_dodge(0.9),width=0.0)+facet_wrap(~Species)+facet_grid(.~Species,labeller=labeller(Species=labels))+
  geom_point(aes(fill=Tissue),color="black",pch=21,size=2.5,position=position_dodge(0.9))+
  ylab("C:N")+
  ylim(5,8.5)+
  scale_fill_manual(values=cols)+newtheme
 
c

#put the plots together to create main text Fig.4 
d<-plot_grid(a,b,c,ncol=1,align="v",axis='2')
quartz()
print(d) #export 500x575


### run seperate models for each species * tissue --- these are reported in the supplemental table S2
###################d15N

#POR
# 
d15N.porT<-lmer(d15N~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),data=subset(chain,Species=="POR"))
summary(d15N.porT)
anova(d15N.porT)
sum(resid(d15N.porT)^2) #22.89

emm = emmeans(d15N.porT, ~ Tissue)
pairs(emm)

#POC
d15N.pocT<-lmer(d15N~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
                data=subset(chain,Species=="POC"))
summary(d15N.pocT)
anova(d15N.pocT)
sum(resid(d15N.pocT)^2) #1O.92

emm = emmeans(d15N.pocT, ~ Nutrient.Level*Tissue)
pairs(emm)

###################d13C

#POR
d13C.porT<-lmer(d13C~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
                data=subset(chain,Species=="POR"))
summary(d13C.porT)
anova(d13C.porT)
sum(resid(d13C.porT)^2) #12.18
emm = emmeans(d13C.porT, ~ Nutrient.Level*Tissue)
pairs(emm)

#POC
d13C.pocT<-lmer(d13C~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
                data=subset(chain,Species=="POC"))
summary(d13C.pocT)
anova(d13C.pocT)
sum(resid(d13C.pocT)^2) #9.31
emm = emmeans(d13C.pocT, ~ Nutrient.Level*Tissue)
pairs(emm)

###################CN
#POR
CN.porT<-lmer(C.N~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
              data=subset(chain,Species=="POR"))
summary(CN.porT)
anova(CN.porT)
sum(resid(CN.porT)^2) #5.92
emm = emmeans(CN.porT, ~ Nutrient.Level*Tissue)
pairs(emm)

#POC
CN.pocT<-lmer(C.N~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
              data=subset(chain,Species=="POC"),REML=FALSE)
summary(CN.pocT)
anova(CN.pocT)
sum(resid(CN.pocT)^2) #18.25
emm = emmeans(CN.pocT, ~ Nutrient.Level*Tissue)
pairs(emm)














