# what can we do using the Drake & Griffen D. magna experiment
# http://doi.org/10.1038/nature09389

load("drake_griffen_data/drake_griffen.RData")

library(mgcv)
library(ggplot2)

# what's going on  here?
populations$Treatment <- "Control"
populations$Treatment[populations$deteriorating] <- "Deteriorating"
p <- ggplot(populations) +
  geom_point(aes(x=day, y=Nhat, colour=Treatment, group=ID)) +
  geom_vline(xintercept=154) +
  geom_text(aes(x=x,y=y, label=label), size=9,
            data=data.frame(x=230, y=300, label="Start deterioration")) +
  labs(y="Count")
print(p)


# fitting two separate models, looks like the start is similar in both
# populations (duh!) but then then end is less so (double duh!)
dmagna_happy <- gam(Nhat~s(day,  k=50), data=pop_happy, method="REML")
dmagna_unhappy <- gam(Nhat~s(day, k=50), data=pop_unhappy, method="REML")

par(mfrow=c(1,2))
plot(dmagna_happy, main="happy zooplankton", shift=mean(pop_happy$Nhat), scale=0, shade=TRUE)
points(pop_happy[,c("day", "Nhat")], pch=19, cex=0.3, alpha=0.3, col="darkgrey")
abline(h=0)
plot(dmagna_unhappy, main="unhappy zooplankton", shift=mean(pop_unhappy$Nhat), scale=0, shade=TRUE)
abline(h=0)
points(pop_unhappy[,c("day", "Nhat")], pch=19, cex=0.3, col="darkgrey")


# let's be FANCY
populations$Treatment <- as.factor(populations$Treatment)

# some cuter plotting...
cute_plotter <- function(model, main=""){
  dat <- data.frame(day = c(seq(0, 416, by=1),seq(0, 416, by=1)),
                    Treatment = c(rep("Deteriorating", 417), rep("Control", 417)))
  dat$Treatment <- as.factor(dat$Treatment)
  pp <- predict(model, newdata=dat, type="response", se.fit=TRUE)

  # deteriorating
  det_ind <- 1:417
  plot(dat$day[det_ind], pp$fit[det_ind], type="l", ylim=c(0,160),
       xlab="Day", ylab="Abundance", main=main, col="red")
  polygon(c(dat$day[det_ind], rev(dat$day[det_ind])),
          c(pp$fit[det_ind]+2*pp$se.fit[det_ind],
            rev(pp$fit[det_ind]-2*pp$se.fit[det_ind])), col=rgb(1,0,0,0.4),
          border=NA)
  # control
  con_ind <- 418:834
  lines(dat$day[con_ind], pp$fit[con_ind], type="l")
  polygon(c(dat$day[con_ind], rev(dat$day[con_ind])),
          c(pp$fit[con_ind]+2*pp$se.fit[con_ind],
            rev(pp$fit[con_ind]-2*pp$se.fit[con_ind])), col=rgb(0,0,0,0.4),
          border=NA)

}

# okay let's get serious about these models and use Eric's grid...

k <- 40
# models that don't make sense
# 1. nope -- need the groups

# Model 2 Global smoother w/ group-specific deviations, with equal smoothness penalties for all groups
dmagna_m2 <- gam(Nhat~s(day, k=k) + s(day, Treatment, bs="fs"),
                 data=populations, select=TRUE, method="REML", family=tw())

# model 3 Global smooth term w/ group-specific deviations, with different smoothness penalties for each group
dmagna_m3 <- gam(Nhat~s(day, k=k) + s(day,by=Treatment)+s(Treatment,bs="re"),
                 data=populations, select=TRUE, method="REML", family=tw())



# model 4 No global smoother, group-level smooths share a penalty
dmagna_m4 <- gam(Nhat~s(day, Treatment, bs="fs", k=k),
                 data=populations, method="REML", family=tw())

# Model 5: No shared information between groups
dmagna_m5 <- gam(Nhat~ s(day, by=Treatment, k=k)+s(Treatment,bs="re"),
                 data=populations, select=TRUE, method="REML", family=tw())




par(mfrow=c(4,5))
cute_plotter(dmagna_m2, main="global, group dev, eq sp")
plot(dmagna_m2)
plot(1:10,1:10, type="n", axes=FALSE, xlab="", ylab="")
plot(1:10,1:10, type="n", axes=FALSE, xlab="", ylab="")
cute_plotter(dmagna_m3, main="global, group dev, diff sp")
plot(dmagna_m3)
cute_plotter(dmagna_m4, main="group dev, eq sp")
plot(dmagna_m4)
plot(1:10,1:10, type="n", axes=FALSE, xlab="", ylab="")
plot(1:10,1:10, type="n", axes=FALSE, xlab="", ylab="")
plot(1:10,1:10, type="n", axes=FALSE, xlab="", ylab="")
cute_plotter(dmagna_m5, main="group dev, diff sp")
plot(dmagna_m5)



# observation: m2/m3 and m4/m5 have v. similar check plots
