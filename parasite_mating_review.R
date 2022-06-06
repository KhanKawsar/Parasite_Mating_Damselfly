
library(bootES) 
library(betareg)
library(yarrr)

## calculate p valuee for beta 

mat <- read.csv("data/parasite_mating_dataset.csv", stringsAsFactors = FALSE)
bmat <- read.csv("data/parasite_mating_dataset_betareg.csv", stringsAsFactors = FALSE)
head(mat,10)
head(bmat,10)

## one sample t-test to test if beta of parasitised male is lower than 0.5   
t.test(mat$parasitised_beta, mu=0.5, alternative="less", conf.level=0.95)

# calculate cohen's d
d <- (0.5-(mat$parasitised_beta))/sd(mat$parasitised_beta)
bootES(d, R = 10000)


## beta reg with intercept only model 

bmat2 <- bmat%>%
  filter(male == "prasitised")
head(bmat2,10)

b2 <- betareg(beta ~ 1, link = "log", data = bmat2)
summary(b2)
plot(b2) ##to check moddel

## to find the mean value of the intercept 
exp(-1.1505) 

## calculate confidence interval 

confint(b2) 
#log transformed confidence interva to fit with the mean
exp(-1.4280495)
exp( -0.8729633)

###plot proportion of parasitism and proportion of mating 
plot(bmat$parasite_proportion, bmat$mate_proportion, pch = 19, cex = 1.5, 
     col = "red", xlab = "Percentage of parasitism", ylab  = "Percentage of mating", 
     ylim = c(0,1), xlim = c(0,1), cex.lab = 1.5, cex.axis = 1.5, las = 1)

### paarsite intensity 

intens <- read.csv("data/intensity_0605.csv")
head(intens)
tail(intens)

mean(intens$parasite_no)
length(intens$parasite_no)

sd(intens$parasite_no)

se <- (sd(intens$parasite_no)/sqrt(length(intens$parasite_no)))
bootES(intens$parasite_no)

## calculate Mann whitney test 

wilcox.test(parasite_no ~ mate_status, data = intens)

##calculate CI of cohens'd 

set.seed(12345)
bootES(intens, R =10000, data.col = "parasite_no", 
            contrast = c("non_mated", "mated"),
            group.col = "mate_status",
            effect.type = "cohens.d", 
       ci.type = "bca", 
       ci.conf = 0.95,
       plot = FALSE)

### plot Figure 


pageWidthLarge<- 7.08661
pageHeightLarge <- pageWidthLarge * 1.5
pagePaper <- 'special'
fontFamily <- 'Times'

tiff("output/Figure1.tiff", width=pageWidthLarge, height= pageHeightLarge, 
     family=fontFamily, units = "in", res = 300)

a <- layout(matrix(c(1,1,1,1,2,2,3,3,0,4,4,0), nrow=3,ncol=4,byrow =  TRUE))
layout.show(a)

##figure 1a and 1b


image<- readJPEG("data/picture1.jpg") 
pixelWidth <- dim(image)[2]
pixelHeight <- dim(image)[1]
par(mar=c(0,0,0,0))
#Setup a plot region with domain and range from 0 to 2, and no axes or labels
plot(NULL, xlim=c(0,pixelWidth), ylim=c(0,pixelHeight), axes=FALSE,xlab= "" , ylab="", asp= 1)
#Draw the raster image across (half y space and whole x space); 
#trial and error based, the photo has long x axis but small y axis
#....rasterImage(image, xleft, ybottom, xright, ytop,)
rasterImage(image, 0,0,pixelWidth,pixelHeight) 
mtext("(a)", side = 3, line = -3, adj= 0.09, cex=1.0)
mtext("(b)", side = 3, line = -3, adj= 0.55, cex=1.0)

##fig 1c
par(mar=c(3.5,5,0,2))
boxplot(mat$parasitised_proportion, mat$parasitised_mated_proportion, 
        ylab = "Proportion of parasitised males", 
        cex.lab= 1.5,cex.axis = 1.5,
        outline= FALSE, las = 1.0,
        border = yarrr::transparent(c("blue", "red"), trans.val = .4), mgp= c(2.5,1,0),col = "white",
        boxwex=0.5,staplewex= 0.2, ylim = c(0,0.4))

stripchart(data.frame(mat$parasitised_proportion, mat$parasitised_mated_proportion), 
           method = "jitter", pch = 19, vertical = TRUE, add = TRUE, cex = 1.2, 
           col = yarrr::transparent(c("blue", "red"), trans.val = .3))


mtext(c("Free", "Copula"), at= c(1, 2), side =1, line = 1.2, cex = 1.5*par()$cex)
mtext("(c)", side = 3, line = 0.5, adj=0.07, cex = 1.5*par()$cex)

## fig 1d 
par(mar=c(5,5,0,2))

plot(mat$parasitised_proportion, mat$parasitised_beta, pch = 19, cex = 1.5, 
     col = yarrr::transparent(c("red"), trans.val = 0.3), xlab = "Proportion of parasitism", ylab  = "Manly beta", 
     ylim = c(0,1), xlim = c(0,0.5), cex.lab = 1.5, cex.axis = 1.5, las = 1)
abline(h= 0.5, col= "black", lwd= 2.0, lty= 2)
mtext("(d)", side = 3, line = 0.5, adj=0.07, cex = 1.5*par()$cex)

##1e

par(mar=c(4,5,2,2))
boxplot(parasite_no ~ mate_status, data = intens,  
        ylab= "Number of parasites", 
        ylim = c(0,10),
        cex.lab= 1.5,cex.axis = 1.5,
        outline= FALSE, las = 1.0,
        border = yarrr::transparent(c( "red", "blue"), trans.val = .4), mgp= c(2.5,1,0),col = "white",
        boxwex=0.5,staplewex= 0.2, xlab = "",
        names = FALSE)

stripchart(parasite_no ~ mate_status, data = intens, 
           method = "jitter", pch = 19, vertical = TRUE, add = TRUE, 
           cex = 1.2, col = yarrr::transparent(c( "red", "blue"), trans.val = .4))

mtext(c( "Copula", "Free"), at= c(1, 2), side =1, line = 1.2, cex = 1.5*par()$cex)
mtext("(e)", side = 3, line = 0.5, adj=0.07, cex = 1.5*par()$cex)

dev.off()


