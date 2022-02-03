library(readxl)
library(RColorBrewer)
library(Hmisc)
library(nlme)
library(afex)
require(lsmeans)
setwd("K:/CCBD/Nellan/CART ependymoma/Data")
source('S:/Shared Material/Shared Code/R/temp_table1.R')
#survival curves

#data prep:
vol_MAF811nsg<- data.frame(read_excel("MRI Tumor Volume and Survival Curves.xls",sheet="MAF811 NSG MRI Tumor Volume"))
vol_MAF928nsg<- data.frame(read_excel("MRI Tumor Volume and Survival Curves.xls",sheet="MAF928 NSG MRI Tumor Volume"))
vol_MAF811nbsgw<- data.frame(read_excel("MRI Tumor Volume and Survival Curves.xls",sheet="MAF811 NBSGW MRI Tumor Volume"))
vol_MAF811nbsgtak<- data.frame(read_excel("MRI Tumor Volume and Survival Curves.xls",sheet="MAF811 NSG TAK MRI Tumor Volume"))

#wide to long:
vol_MAF811nsg_l<-reshape(vol_MAF811nsg,
                    direction = "long",
                    varying = c(list(names(vol_MAF811nsg)[3:7])),
                    v.names = c("Volume"),
                    idvar = c("Mouse"),
                    timevar = "Days",
                    times = c(0,11,31,45,73))
vol_MAF811nsg_l<-vol_MAF811nsg_l[order(vol_MAF811nsg_l$Mouse,vol_MAF811nsg_l$Days),]
vol_MAF811nsg_l$lty<-ifelse(vol_MAF811nsg_l$Treatment=="HER2 CAR",1,2)
vol_MAF811nsg_l$col<-ifelse(vol_MAF811nsg_l$Treatment=="HER2 CAR",brewer.pal(3,"Set1")[1],brewer.pal(3,"Set1")[2])

vol_MAF928nsg_l<-reshape(vol_MAF928nsg,
                         direction = "long",
                         varying = c(list(names(vol_MAF928nsg)[3:6])),
                         v.names = c("Volume"),
                         idvar = c("Mouse"),
                         timevar = "Days",
                         times = c(0,30,51,81))
vol_MAF928nsg_l<-vol_MAF928nsg_l[order(vol_MAF928nsg_l$Mouse,vol_MAF928nsg_l$Days),]
vol_MAF928nsg_l$lty<-ifelse(vol_MAF928nsg_l$Treatment=="HER2 CAR",1,2)
vol_MAF928nsg_l$col<-ifelse(vol_MAF928nsg_l$Treatment=="HER2 CAR",brewer.pal(3,"Set1")[1],brewer.pal(3,"Set1")[2])

vol_MAF811nbsgw_l<-reshape(vol_MAF811nbsgw,
                         direction = "long",
                         varying = c(list(names(vol_MAF811nbsgw)[3:6])),
                         v.names = c("Volume"),
                         idvar = c("Mouse"),
                         timevar = "Days",
                         times = c(0,14,28,49))
vol_MAF811nbsgw_l<-vol_MAF811nbsgw_l[order(vol_MAF811nbsgw_l$Mouse,vol_MAF811nbsgw_l$Days),]
vol_MAF811nbsgw_l$lty<-ifelse(vol_MAF811nbsgw_l$Treatment=="HER2 CAR",1,2)
vol_MAF811nbsgw_l$col<-ifelse(vol_MAF811nbsgw_l$Treatment=="HER2 CAR",brewer.pal(3,"Set1")[1],brewer.pal(3,"Set1")[2])

vol_MAF811nbsgtak_l<-reshape(vol_MAF811nbsgtak,
                         direction = "long",
                         varying = c(list(names(vol_MAF811nbsgtak)[3:6])),
                         v.names = c("Volume"),
                         idvar = c("Mouse"),
                         timevar = "Days",
                         times = c(0,11,25,40))
vol_MAF811nbsgtak_l<-vol_MAF811nbsgtak_l[order(vol_MAF811nbsgtak_l$Mouse,vol_MAF811nbsgtak_l$Days),]
vol_MAF811nbsgtak_l$lty<-ifelse(vol_MAF811nbsgtak_l$Treatment=="HER2 CAR",1,
                                ifelse(vol_MAF811nbsgtak_l$Treatment=="TAK-779",2,3))
vol_MAF811nbsgtak_l$col<-ifelse(vol_MAF811nbsgtak_l$Treatment=="HER2 CAR",brewer.pal(3,"Set2")[1],
                                ifelse(vol_MAF811nbsgtak_l$Treatment=="TAK-779",brewer.pal(3,"Set2")[2],brewer.pal(3,"Set2")[3]))


####MOUSE TYPE 1
vol_MAF811nsg_mod<-subset(vol_MAF811nsg_l,!is.na(vol_MAF811nsg_l$Volume))
vol_MAF811nsg_mod<-subset(vol_MAF811nsg_mod,vol_MAF811nsg_mod$Days<73)
mod_MAF811nsg<-lme(Volume~factor(Days)*Treatment
                     ,random=~1|Mouse,data=vol_MAF811nsg_mod)

MAF811nsg_sum<-summary(mod_MAF811nsg)
MAF811nsg_sum<-MAF811nsg_sum$tTable[,c(1,2,5)]
MAF811nsg_sum<-as.data.frame(MAF811nsg_sum)
MAF811nsg_sum$Value<-round(MAF811nsg_sum$Value,3)
MAF811nsg_sum$Std.Error<-round(MAF811nsg_sum$Std.Error,3)
MAF811nsg_sum$`p-value`<-round(MAF811nsg_sum$`p-value`,3)

anova_MAF811nsg<-anova(mod_MAF811nsg)
#test for differences between groups at each time point:
ref_MAF811nsg <- lsmeans(mod_MAF811nsg, c("Treatment", "Days"))

c_list_MAF811nsg <- list(c_0 = c(-1, 1, 0, 0, 0, 0, 0, 0),
                 c_11 = c(0, 0, -1, 1, 0, 0, 0, 0),
                 c_31 = c(0, 0, 0, 0, -1, 1, 0, 0),
                 c_45 = c(0, 0, 0, 0, 0, 0, -1, 1),
                 c_firstlast = c(1, -1, 0, 0, 0, 0, -1, 1)
)
contrasts_MAF811nsg<-summary(contrast(ref_MAF811nsg, c_list_MAF811nsg))
contrasts_MAF811nsg<-contrasts_MAF811nsg[,c(1,2,3,6)]
contrasts_MAF811nsg$estimate<-round(contrasts_MAF811nsg$estimate,3)
contrasts_MAF811nsg$SE<-round(contrasts_MAF811nsg$SE,3)
contrasts_MAF811nsg$p.value<-round(contrasts_MAF811nsg$p.value,3)
contrasts_MAF811nsg$p.value_corr<-contrasts_MAF811nsg$p.value*nrow(contrasts_MAF811nsg)

####MOUSE TYPE 2
vol_MAF928nsg_mod<-subset(vol_MAF928nsg_l,!is.na(vol_MAF928nsg_l$Volume))
vol_MAF928nsg_mod<-subset(vol_MAF928nsg_mod,vol_MAF928nsg_mod$Days<81)
mod_MAF928nsg<-lme(Volume~factor(Days)*Treatment
                   ,random=~1|Mouse,data=vol_MAF928nsg_mod)

MAF928nsg_sum<-summary(mod_MAF928nsg)
MAF928nsg_sum<-MAF928nsg_sum$tTable[,c(1,2,5)]
MAF928nsg_sum<-as.data.frame(MAF928nsg_sum)
MAF928nsg_sum$Value<-round(MAF928nsg_sum$Value,3)
MAF928nsg_sum$Std.Error<-round(MAF928nsg_sum$Std.Error,3)
MAF928nsg_sum$`p-value`<-round(MAF928nsg_sum$`p-value`,3)

anova_MAF928nsg<-anova(mod_MAF928nsg)
#test for differences between groups at each time point:
ref_MAF928nsg <- lsmeans(mod_MAF928nsg, c("Treatment", "Days"))

c_list_MAF928nsg <- list(c_0 = c(-1, 1, 0, 0, 0, 0),
                         c_30 = c(0, 0, -1, 1, 0, 0),
                         c_51 = c(0, 0, 0, 0, -1, 1),
                         c_firstlast = c(1, -1, 0, 0, -1, 1)
)
contrasts_MAF928nsg<-summary(contrast(ref_MAF928nsg, c_list_MAF928nsg))
contrasts_MAF928nsg<-contrasts_MAF928nsg[,c(1,2,3,6)]
contrasts_MAF928nsg$estimate<-round(contrasts_MAF928nsg$estimate,3)
contrasts_MAF928nsg$SE<-round(contrasts_MAF928nsg$SE,3)
contrasts_MAF928nsg$p.value<-round(contrasts_MAF928nsg$p.value,3)
contrasts_MAF928nsg$p.value_corr<-contrasts_MAF928nsg$p.value*nrow(contrasts_MAF928nsg)

####MOUSE TYPE 3
vol_MAF811nbsgw_mod<-subset(vol_MAF811nbsgw_l,!is.na(vol_MAF811nbsgw_l$Volume))
vol_MAF811nbsgw_mod<-subset(vol_MAF811nbsgw_mod,vol_MAF811nbsgw_mod$Days<49)
mod_MAF811nbsgw<-lme(Volume~factor(Days)*Treatment
                   ,random=~1|Mouse,data=vol_MAF811nbsgw_mod)

MAF811nbsgw_sum<-summary(mod_MAF811nbsgw)
MAF811nbsgw_sum<-MAF811nbsgw_sum$tTable[,c(1,2,5)]
MAF811nbsgw_sum<-as.data.frame(MAF811nbsgw_sum)
MAF811nbsgw_sum$Value<-round(MAF811nbsgw_sum$Value,3)
MAF811nbsgw_sum$Std.Error<-round(MAF811nbsgw_sum$Std.Error,3)
MAF811nbsgw_sum$`p-value`<-round(MAF811nbsgw_sum$`p-value`,3)

anova_MAF811nbsgw<-anova(mod_MAF811nbsgw)
#test for differences between groups at each time point:
ref_MAF811nbsgw <- lsmeans(mod_MAF811nbsgw, c("Treatment", "Days"))

c_list_MAF811nbsgw <- list(c_0 = c(-1, 1, 0, 0, 0, 0),
                         c_14 = c(0, 0, -1, 1, 0, 0),
                         c_28 = c(0, 0, 0, 0, -1, 1),
                         c_firstlast = c(1, -1, 0, 0, -1, 1)
)
contrasts_MAF811nbsgw<-summary(contrast(ref_MAF811nbsgw, c_list_MAF811nbsgw))
contrasts_MAF811nbsgw<-contrasts_MAF811nbsgw[,c(1,2,3,6)]
contrasts_MAF811nbsgw$estimate<-round(contrasts_MAF811nbsgw$estimate,3)
contrasts_MAF811nbsgw$SE<-round(contrasts_MAF811nbsgw$SE,3)
contrasts_MAF811nbsgw$p.value<-round(contrasts_MAF811nbsgw$p.value,3)
contrasts_MAF811nbsgw$p.value_corr<-contrasts_MAF811nbsgw$p.value*nrow(contrasts_MAF811nbsgw)


####MOUSE TYPE 4
vol_MAF811nbsgtak_mod<-subset(vol_MAF811nbsgtak_l,!is.na(vol_MAF811nbsgtak_l$Volume))
# vol_MAF811nbsgtak_mod<-subset(vol_MAF811nbsgtak_mod,vol_MAF811nbsgtak_mod$Days<49)
mod_MAF811nbsgtak<-lme(Volume~factor(Days)*Treatment
                     ,random=~1|Mouse,data=vol_MAF811nbsgtak_mod)

MAF811nbsgtak_sum<-summary(mod_MAF811nbsgtak)
MAF811nbsgtak_sum<-MAF811nbsgtak_sum$tTable[,c(1,2,5)]
MAF811nbsgtak_sum<-as.data.frame(MAF811nbsgtak_sum)
MAF811nbsgtak_sum$Value<-round(MAF811nbsgtak_sum$Value,3)
MAF811nbsgtak_sum$Std.Error<-round(MAF811nbsgtak_sum$Std.Error,3)
MAF811nbsgtak_sum$`p-value`<-round(MAF811nbsgtak_sum$`p-value`,3)

anova_MAF811nbsgtak<-anova(mod_MAF811nbsgtak)
#test for differences between groups at each time point:
ref_MAF811nbsgtak <- lsmeans(mod_MAF811nbsgtak, c("Treatment", "Days"))

##########PICK UP HERE############
c_list_MAF811nbsgtak <- list(
                             # c_11.1 = c( 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0),
                             # c_11.2 = c( 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0),
                             # c_11.3 = c( 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0),
                             # 
                             # c_25.1 = c(0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0),
                             # c_25.2 = c(0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0),
                             # c_25.3 = c(0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0),
                             # 
                             # c_40.1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0),
                             # c_40.2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1),
                             # c_40.3 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1),
                             # #1, -1, 0, 0, -1, 1
                             c_fl.1 = c(1, -1, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0),
                             c_fl.2 = c(0, 1, -1, 0, 0, 0, 0, 0, 0, 0, -1, 1),
                             c_fl.3 = c(1, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 1))

contrasts_MAF811nbsgtak<-summary(contrast(ref_MAF811nbsgtak, c_list_MAF811nbsgtak))
contrasts_MAF811nbsgtak<-contrasts_MAF811nbsgtak[,c(1,2,3,6)]
contrasts_MAF811nbsgtak$estimate<-round(contrasts_MAF811nbsgtak$estimate,3)
contrasts_MAF811nbsgtak$SE<-round(contrasts_MAF811nbsgtak$SE,3)
contrasts_MAF811nbsgtak$p.value<-round(contrasts_MAF811nbsgtak$p.value,3)
contrasts_MAF811nbsgtak$p.value_corr<-contrasts_MAF811nbsgtak$p.value*nrow(contrasts_MAF811nbsgtak)


#plot for predicted lines:
ref_MAF811nsg<-data.frame(ref_MAF811nsg)
ref_MAF928nsg<-data.frame(ref_MAF928nsg)
ref_MAF811nbsgw<-data.frame(ref_MAF811nbsgw)
ref_MAF811nbsgtak<-data.frame(ref_MAF811nbsgtak)

jpeg("K:/CCBD/Nellan/CART ependymoma/Results/Plots/Tumor_volume.jpg",width=9,height=9,units="in",
     res=100)
par(mfrow=c(2,2),mar=c(3,3,2,2))
plot(c(0,45),c(0,65),type="n",main="MAF811nsg",xaxt="n",xlab="",ylab="")
axis(1,at=unique(vol_MAF811nsg_mod$Days))
mtext("Tumor volume",side=2,line=2)
mtext("Days",side=1,line=2)
for (i in unique(vol_MAF811nsg_mod$Mouse)){
  lines(vol_MAF811nsg_mod$Days[vol_MAF811nsg_mod$Mouse==i],vol_MAF811nsg_mod$Volume[vol_MAF811nsg_mod$Mouse==i],
        #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
        col=vol_MAF811nsg_mod$col[vol_MAF811nsg_mod$Mouse==i],lty=2)
}
lines(ref_MAF811nsg$Days[ref_MAF811nsg$Treatment=="CD19 CAR"],
      ref_MAF811nsg$lsmean[ref_MAF811nsg$Treatment=="CD19 CAR"],
        #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
        col="#377EB8",lwd=2)
points(ref_MAF811nsg$Days[ref_MAF811nsg$Treatment=="CD19 CAR"],
      ref_MAF811nsg$lsmean[ref_MAF811nsg$Treatment=="CD19 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#377EB8",pch=19)
lines(ref_MAF811nsg$Days[ref_MAF811nsg$Treatment=="HER2 CAR"],
      ref_MAF811nsg$lsmean[ref_MAF811nsg$Treatment=="HER2 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#E41A1C",lwd=2)
points(ref_MAF811nsg$Days[ref_MAF811nsg$Treatment=="HER2 CAR"],
       ref_MAF811nsg$lsmean[ref_MAF811nsg$Treatment=="HER2 CAR"],
       #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
       col="#E41A1C",pch=19)
contrasts_MAF811nsg
text(11,55,"*",font=2,cex=1.5)
text(31,55,"*",font=2,cex=1.5)
text(45,55,"*",font=2,cex=1.5)
points(c(0,45),c(60,60),pch=3)
lines(c(0,45),c(60,60),pch=3)
text(45/2,62,"**",cex=1.5)
# legend("topleft",levels(as.factor(vol_MAF811nbsgw_l$Treatment)), col=brewer.pal(3,"Set1")[c(1,2)],
#        lty=1,lwd=2)

plot(c(0,51),c(0,45),type="n",main="MAF928nsg",xaxt="n",xlab="Days",ylab="Tumor volume")
axis(1,at=unique(vol_MAF928nsg_mod$Days))
mtext("Tumor volume",side=2,line=2)
mtext("Days",side=1,line=2)
for (i in unique(vol_MAF928nsg_mod$Mouse)){
  lines(vol_MAF928nsg_mod$Days[vol_MAF928nsg_mod$Mouse==i],vol_MAF928nsg_mod$Volume[vol_MAF928nsg_mod$Mouse==i],
        #lty=vol_MAF928nsg_l$lty[vol_MAF928nsg_l$Mouse==i],
        col=vol_MAF928nsg_mod$col[vol_MAF928nsg_mod$Mouse==i],lty=2)
}
lines(ref_MAF928nsg$Days[ref_MAF928nsg$Treatment=="CD19 CAR"],
      ref_MAF928nsg$lsmean[ref_MAF928nsg$Treatment=="CD19 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#377EB8",lwd=2)
points(ref_MAF928nsg$Days[ref_MAF928nsg$Treatment=="CD19 CAR"],
      ref_MAF928nsg$lsmean[ref_MAF928nsg$Treatment=="CD19 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#377EB8",pch=19)
lines(ref_MAF928nsg$Days[ref_MAF928nsg$Treatment=="HER2 CAR"],
      ref_MAF928nsg$lsmean[ref_MAF928nsg$Treatment=="HER2 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#E41A1C",lwd=2)
points(ref_MAF928nsg$Days[ref_MAF928nsg$Treatment=="HER2 CAR"],
      ref_MAF928nsg$lsmean[ref_MAF928nsg$Treatment=="HER2 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#E41A1C",pch=19)
contrasts_MAF928nsg
text(51,38,"*",font=2,cex=1.5)
points(c(0,51),c(41,41),pch=3)
lines(c(0,51),c(41,41),pch=3)
text(51/2,43,"**",cex=1.5)

plot(c(0,28),c(0,6),type="n",main="MAF811nbsgw",xaxt="n",xlab="Days",ylab="Tumor volume")
axis(1,at=unique(vol_MAF811nbsgw_mod$Days))
mtext("Tumor volume",side=2,line=2)
mtext("Days",side=1,line=2)

for (i in unique(vol_MAF811nbsgw_mod$Mouse)){
  lines(vol_MAF811nbsgw_mod$Days[vol_MAF811nbsgw_mod$Mouse==i],vol_MAF811nbsgw_mod$Volume[vol_MAF811nbsgw_mod$Mouse==i],
        #lty=vol_MAF811nbsgw_l$lty[vol_MAF811nbsgw_l$Mouse==i],
        col=vol_MAF811nbsgw_mod$col[vol_MAF811nbsgw_mod$Mouse==i],lty=2)
}
lines(ref_MAF811nbsgw$Days[ref_MAF811nbsgw$Treatment=="CD19 CAR"],
      ref_MAF811nbsgw$lsmean[ref_MAF811nbsgw$Treatment=="CD19 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#377EB8",lwd=2)
points(ref_MAF811nbsgw$Days[ref_MAF811nbsgw$Treatment=="CD19 CAR"],
      ref_MAF811nbsgw$lsmean[ref_MAF811nbsgw$Treatment=="CD19 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#377EB8",pch=19)
lines(ref_MAF811nbsgw$Days[ref_MAF811nbsgw$Treatment=="HER2 CAR"],
      ref_MAF811nbsgw$lsmean[ref_MAF811nbsgw$Treatment=="HER2 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#E41A1C",lwd=2)
points(ref_MAF811nbsgw$Days[ref_MAF811nbsgw$Treatment=="HER2 CAR"],
      ref_MAF811nbsgw$lsmean[ref_MAF811nbsgw$Treatment=="HER2 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col="#E41A1C",pch=19)
contrasts_MAF811nbsgw
points(c(0,28),c(5.5,5.5),pch=3)
lines(c(0,28),c(5.5,5.5),pch=3)
text(28/2,5.7,"**",cex=1.5)

plot(c(0,40),c(0,20),type="n",main="MAF811nbsgtak",xaxt="n",xlab="Days",ylab="Tumor volume")
axis(1,at=unique(vol_MAF811nbsgtak_mod$Days))
mtext("Tumor volume",side=2,line=2)
mtext("Days",side=1,line=2)
for (i in unique(vol_MAF811nbsgtak_mod$Mouse)){
  lines(vol_MAF811nbsgtak_mod$Days[vol_MAF811nbsgtak_mod$Mouse==i],vol_MAF811nbsgtak_mod$Volume[vol_MAF811nbsgtak_mod$Mouse==i],
        #lty=vol_MAF811nbsgtak_l$lty[vol_MAF811nbsgtak_l$Mouse==i],
        col=vol_MAF811nbsgtak_mod$col[vol_MAF811nbsgtak_mod$Mouse==i],lty=2)
}
lines(ref_MAF811nbsgtak$Days[ref_MAF811nbsgtak$Treatment=="COMBO"],
      ref_MAF811nbsgtak$lsmean[ref_MAF811nbsgtak$Treatment=="COMBO"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col=vol_MAF811nbsgtak_l$col[vol_MAF811nbsgtak_l$Treatment=="COMBO"][1],lwd=2)
points(ref_MAF811nbsgtak$Days[ref_MAF811nbsgtak$Treatment=="CD19 CAR"],
       ref_MAF811nbsgtak$lsmean[ref_MAF811nbsgtak$Treatment=="CD19 CAR"],
       #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
       col=vol_MAF811nbsgtak_l$col[vol_MAF811nbsgtak_l$Treatment=="COMBO"][1],pch=19)
lines(ref_MAF811nbsgtak$Days[ref_MAF811nbsgtak$Treatment=="HER2 CAR"],
      ref_MAF811nbsgtak$lsmean[ref_MAF811nbsgtak$Treatment=="HER2 CAR"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col=vol_MAF811nbsgtak_l$col[vol_MAF811nbsgtak_l$Treatment=="HER2 CAR"][1],lwd=2)
points(ref_MAF811nbsgtak$Days[ref_MAF811nbsgtak$Treatment=="HER2 CAR"],
       ref_MAF811nbsgtak$lsmean[ref_MAF811nbsgtak$Treatment=="HER2 CAR"],
       #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
       col=vol_MAF811nbsgtak_l$col[vol_MAF811nbsgtak_l$Treatment=="HER2 CAR"][1],pch=19)
lines(ref_MAF811nbsgtak$Days[ref_MAF811nbsgtak$Treatment=="TAK-779"],
      ref_MAF811nbsgtak$lsmean[ref_MAF811nbsgtak$Treatment=="TAK-779"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col=vol_MAF811nbsgtak_l$col[vol_MAF811nbsgtak_l$Treatment=="TAK-779"][1],lwd=2)
points(ref_MAF811nbsgtak$Days[ref_MAF811nbsgtak$Treatment=="TAK-779"],
      ref_MAF811nbsgtak$lsmean[ref_MAF811nbsgtak$Treatment=="TAK-779"],
      #lty=vol_MAF811nsg_l$lty[vol_MAF811nsg_l$Mouse==i],
      col=vol_MAF811nbsgtak_l$col[vol_MAF811nbsgtak_l$Treatment=="TAK-779"][1],pch=19)
contrasts_MAF811nbsgtak

points(c(0,40),c(18.5,18.5),pch=3)
lines(c(0,40),c(18.5,18.5),pch=3)
text(40/2,19.3,"***",cex=1.5)
# legend("topleft",levels(as.factor(vol_MAF811nbsgtak_l$Treatment)), col=brewer.pal(3,"Set2"),lty=1,
#        lwd=2)
dev.off()
