

## ANALYSES FOR:
# Reynolds AZ, Wander K, Sum CY, Su M, Emery Thompson M, Hooper PL, Li H, Shenk MK, Starkweather KE, Blumefield T, Mattison SM. (in review) Matriliny reverses gender disparities in inflammation and hypertension among the Mosuo of China. 


#-------------------------------------------------
# READ DATA FROM FILE

data <- read.csv("Dataset_S1.csv")



#-------------------------------------------------
# BAYESIAN LOGISTIC REGRESSION MODELS 

require(rethinking)



#------------------------------------
# SUBSET ONLY COMPLETE DATA
d2 <- data[!is.na(data$elevatedCRP),]

# CRP MAIN MODEL, REF: PAT FEMALES
start <- list(
	a = mean(d2$elevatedCRP),
	b.Age = 0,
	b.Male = 0,
	b.Matriliny = 0,
	b.Male.Matriliny = 0
)
flist <- alist(
	elevatedCRP ~ dbinom(1,p),
	logit(p) <- a + b.Male*Male + b.Matriliny*Matriliny + b.Male.Matriliny*Male*Matriliny + b.Age*Age,
	a ~ dnorm(-10,10),
	b.Age ~ dnorm(0,10),
	b.Male ~ dnorm(0,10),
	b.Matriliny ~ dnorm(0,10),
	b.Male.Matriliny ~ dnorm(0,10)
)
mCRP.patF <- map(flist, data=d2, start=start)
precis(mCRP.patF, prob=0.95)
post.patF <- extract.samples(mCRP.patF, 1e4)
precis(post.patF, prob=.95)
dens(post.patF)

# what is the posterior prob that the age param is < 0?
sum(post.patF$b.Age < 0) / 1e4
# what is the posterior probability that the Mat parameter is > zero?
sum(post.patF$b.Matriliny > 0) / 1e4 
# what is the posterior prob that male param is > 0
sum(post.patF$b.Male > 0) / 1e4
# what is the posterior prob that the male-Matriliny param is < 0?
sum(post.patF$b.Male.Matriliny < 0) / 1e4




#---------------------------------
# CRP HOH MODEL, REF. PAT FEMALES

## SUBSET DATA
d3 <- d2
d3 <- d3[!is.na(d3$Head),]

start <- list(
	a = mean(d2$elevatedCRP),
	b.Age = 0,
	b.Male = 0,
	b.Matriliny = 0,
	b.Male.Matriliny = 0,
	b.Head = 0)
flist <- alist(
	elevatedCRP ~ dbinom(1,p),
	logit(p) <- a + b.Age*Age + b.Head*Head + b.Male*Male + b.Matriliny*Matriliny + b.Male.Matriliny*Male*Matriliny,
	a ~ dnorm(-10,10),
	b.Age ~ dnorm(0,10),
	b.Male ~ dnorm(0,10),
	b.Matriliny ~ dnorm(0,10),
	b.Male.Matriliny ~ dnorm(0,10),	
	b.Head ~ dnorm(0,10))
mCRP.hoh <- map(flist, data=d3, start=start)
precis(mCRP.hoh, prob=0.90)
post.hoh <- extract.samples(mCRP.hoh, 1e4)
precis(post.hoh, prob=.95)
# dens(post.hoh)

sum(post.hoh$b.Age < 0) / 1e4
sum(post.hoh$b.Male > 0) / 1e4
sum(post.hoh$b.Matriliny > 0) / 1e4
sum(post.hoh$b.Male.Matriliny < 0) / 1e4
sum(post.hoh$b.Head > 0) / 1e4







#----------------------------------
# HTN MAIN MODEL, REF: PAT FEMALES
dataHTN <- data[!is.na(data$htn) & !is.na(data$Male) & !is.na(data$Matriliny) & !is.na(data$Age),]

with(dataHTN, table(htn, Male, Matriliny))

start <- list(
	a = mean(dataHTN$htn),	
	b.Age = 0,
	b.Male = 0,
	b.Matriliny = 0,
	b.Male.Matriliny = 0
)
flist <- alist(
	htn ~ dbinom(1,p),
	logit(p) <- a + b.Age*Age + b.Male*Male + b.Matriliny*Matriliny + b.Male.Matriliny*Male*Matriliny,
	a ~ dnorm(-10,10),
	b.Age ~ dnorm(0,10),
	b.Male ~ dnorm(0,10),
	b.Matriliny ~ dnorm(0,10),
	b.Male.Matriliny ~ dnorm(0,10)
)
mHTN <- map(flist, data=dataHTN, start=start)

precis(mHTN, prob=0.95)
post.htn <- extract.samples(mHTN, 1e4)
precis(post.htn, prob=.95)
dens(post.htn)
# HPDI(post.htn$b.Pat, prob=.90)  # At prob=0.90 is the highest we can go while excluding the possibility of zero for the pat term

# what is the posterior prob that the Age param is > 0?
sum(post.htn$b.Age > 0) / 1e4
sum(post.htn$b.Age < 0) / 1e4

# what is the posterior prob that male param is > 0
sum(post.htn$b.Male > 0) / 1e4
sum(post.htn$b.Male < 0) / 1e4	

# what is the posterior probability that the Matriliny parameter is > zero?
sum(post.htn$b.Matriliny >  0) / 1e4 
sum(post.htn$b.Matriliny < 0) / 1e4 

# what is the posterior prob that the male-pat param is < 0?
sum(post.htn$b.Male.Matriliny < 0) / 1e4
sum(post.htn$b.Male.Matriliny > 0) / 1e4


#---------------------------------------
# HTN HOH MODEL, REF: PAT FEMALES

dataHTN3 <- dataHTN[!is.na(dataHTN$Head),]

start <- list(
	a = mean(dataHTN$htn),
	b.Age = 0,
	b.Male = 0,
	b.Matriliny = 0,
	b.Male.Matriliny = 0,
	b.Head = 0)
flist <- alist(
	htn ~ dbinom(1,p),
	logit(p) <- a + b.Age*Age + b.Head*Head + b.Male*Male + b.Matriliny*Matriliny + b.Male.Matriliny*Male*Matriliny,
	a ~ dnorm(-10,10),
	b.Age ~ dnorm(0,10),
	b.Male ~ dnorm(0,10),
	b.Matriliny ~ dnorm(0,10),
	b.Male.Matriliny ~ dnorm(0,10),	
	b.Head ~ dnorm(0,10))
mHTN.head4 <- map(flist, data=dataHTN3, start=start)
precis(mHTN.head4, prob=0.90)
post4 <- extract.samples(mHTN.head4, 1e4)
precis(post4, prob=.95)

sum(post4$b.Age < 0) / 1e4
sum(post4$b.Male > 0) / 1e4
sum(post4$b.Matriliny > 0) / 1e4
sum(post4$b.Male.Matriliny < 0) / 1e4
sum(post4$b.Head < 0) / 1e4



### TABLE S1
#---------------------------------
# CRP BMI MODEL, REF: PAT FEMALE

d6 <- d2[!is.na(d2$OverObese) & !is.na(d2$Underweight),]

start <- list(
	a = mean(d6$elevatedCRP),
	b.Age = 0,
	b.Male = 0,
	b.Matriliny = 0,
	b.Male.Matriliny = 0,
	b.OverObese = 0,
	b.Underweight = 0
)
flist <- alist(
	elevatedCRP ~ dbinom(1,p),
	logit(p) <- a + b.Male*Male + b.Matriliny*Matriliny + b.Male.Matriliny*Male*Matriliny + b.Age*Age + b.OverObese*OverObese + b.Underweight*Underweight,
	a ~ dnorm(-10,10),
	b.Age ~ dnorm(0,10),
	b.Male ~ dnorm(0,10),
	b.Matriliny ~ dnorm(0,10),
	b.Male.Matriliny ~ dnorm(0,10),
	b.OverObese ~ dnorm(0,10),
	b.Underweight ~ dnorm(0,10)
)
mCRP.met <- map(flist, data=d6, start=start)
precis(mCRP.met, prob=0.95)
post.met <- extract.samples(mCRP.met, 1e4)
precis(post.met)

# what is the posterior prob that the Age param is > 0?
sum(post.met$b.Age > 0) / 1e4
sum(post.met$b.Age < 0) / 1e4

# what is the posterior probability that the pat parameter is > zero?
sum(post.met$b.Matriliny > 0) / 1e4 
sum(post.met$b.Matriliny < 0) / 1e4 

# what is the posterior prob that male param is > 0
sum(post.met$b.Male > 0) / 1e4
sum(post.met$b.Male < 0) / 1e4

# what is the posterior prob that the male-pat param is < 0?
sum(post.met$b.Male.Matriliny < 0) / 1e4
sum(post.met$b.Male.Matriliny > 0) / 1e4

# what is the posterior prob that the OverObese param is > 0?
sum(post.met$b.OverObese > 0) / 1e4
sum(post.met$b.OverObese < 0) / 1e4

# what is the posterior prob that the Underweight param is > 0?
sum(post.met$b.Underweight > 0) / 1e4
sum(post.met$b.Underweight < 0) / 1e4



#-------------------------------
# HTN BMI MODEL, REF: PAT FEMALE
 
dataHTN6 <- dataHTN[!is.na(dataHTN$OverObese) & !is.na(dataHTN$Underweight),]

start <- list(
	a = mean(dataHTN6$htn),
	b.Age = 0,
	b.Male = 0,
	b.Matriliny = 0,
	b.Male.Matriliny = 0,
	b.OverObese = 0,
	b.Underweight = 0
)
flist <- alist(
	htn ~ dbinom(1,p),
	logit(p) <- a + b.Male*Male + b.Matriliny*Matriliny + b.Male.Matriliny*Male*Matriliny + b.Age*Age + b.OverObese*OverObese + b.Underweight*Underweight,
	a ~ dnorm(-10,10),
	b.Age ~ dnorm(0,10),
	b.Male ~ dnorm(0,10),
	b.Matriliny ~ dnorm(0,10),
	b.Male.Matriliny ~ dnorm(0,10),
	b.OverObese ~ dnorm(0,10),
	b.Underweight ~ dnorm(0,10)
)
mHTN.bmi <- map(flist, data=dataHTN6, start=start)

precis(mHTN.bmi, prob=0.95)
post.htn.bmi <- extract.samples(mHTN.bmi, 1e4)
precis(post.htn.bmi, prob=.95)
# dens(post.htn)

# what is the posterior prob that the Age param is > 0?
sum(post.htn.bmi$b.Age > 0) / 1e4
sum(post.htn.bmi$b.Age < 0) / 1e4

# what is the posterior probability that the pat parameter is > zero?
sum(post.htn.bmi$b.Matriliny >  0) / 1e4 
sum(post.htn.bmi$b.Matriliny < 0) / 1e4 

# what is the posterior prob that male param is > 0
sum(post.htn.bmi$b.Male > 0) / 1e4
sum(post.htn.bmi$b.Male < 0) / 1e4

# what is the posterior prob that the male-pat param is < 0?
sum(post.htn.bmi$b.Male.Matriliny < 0) / 1e4
sum(post.htn.bmi$b.Male.Matriliny > 0) / 1e4

# what is the posterior prob that the OverObese param is > 0?
sum(post.htn.bmi$b.OverObese > 0) / 1e4
sum(post.htn.bmi$b.OverObese < 0) / 1e4

# what is the posterior prob that the Underweight param is > 0?
sum(post.htn.bmi$b.Underweight > 0) / 1e4
sum(post.htn.bmi$b.Underweight < 0) / 1e4







#------------------------------------------
# PREDICT AND DRAW FIGURE 1

require(ggplot2)
require(ggsci)
require(ggpubr)


### CRP, MAIN MODEL

new_data <- data.frame(
	Age = 45,
	Male = c(0, 1, 0, 1),
	Matriliny = c(0, 0, 1, 1)
)

p <- link(mCRP.patF, data=new_data)	
p_mean <- apply(p,2,mean) 
p_HPDI <- apply(p,2,HPDI,prob=0.68) 
p2 <- rbind(p_mean,p_HPDI)

pdata <- new_data
p2 <- t(p2)
colnames(p2) <- c("mean", "lower", "upper")
preds <- cbind(pdata, p2)

preds$Gender <- c("Women", "Men", "Women", "Men")
preds$Matriliny <- as.factor(preds$Matriliny)


### HTN, MAIN MODEL

new_data <- data.frame(
	Age = 45,
	Male = c(0, 1, 0, 1),
	Matriliny = c(0, 0, 1, 1)
)

p <- link(mHTN, data=new_data)	
p_mean <- apply(p,2,mean) 
p_HPDI <- apply(p,2,HPDI,prob=0.68) 
p2 <- rbind(p_mean,p_HPDI)

pdata <- new_data
p2 <- t(p2)
colnames(p2) <- c("mean", "lower", "upper")
preds.htn <- cbind(pdata, p2)

preds.htn$Gender <- c("Women", "Men", "Women", "Men")
preds.htn$Matriliny <- as.factor(preds.htn$Matriliny)


### AND PLOT

pd <- position_dodge(.2)

dev.new(height=5, width=11)

CRP <- ggplot(preds, aes(x=Matriliny, y=mean, colour=Gender, group=Gender)) +
	geom_errorbar(aes(ymin=lower, ymax=upper), size=.9, width=.2, position=pd) +
	geom_line(position=pd, size=.9) +
	geom_point(position=pd, size=3, aes(shape=Gender)) + 
	ggtitle(element_blank()) + 
	scale_x_discrete(labels=c("Patrilineal", "Matrilineal")) + 
	ylab("Predicted probability of CRP 3-5 mg/L") + 
	# scale_y_continuous(limits=c(0, .10)) + 
	scale_color_aaas() +
	scale_shape_manual(values=c(17, 16, 17, 16)) + 
	theme_classic(base_size=16) +  
	theme(legend.position = "none",
		plot.margin=unit(c(.5,.5,.5,.5), "cm"),
		axis.title.y=element_text(margin=margin(0,15,0,0)),
		axis.title.x=element_blank())
		
HTN <- ggplot(preds.htn, aes(x=Matriliny, y=mean, colour=Gender, group=Gender)) +
	geom_errorbar(aes(ymin=lower, ymax=upper), size=.9, width=.2, position=pd) +
	geom_line(position=pd, size=.9) +
	geom_point(position=pd, size=3, aes(shape=Gender)) + 
	ggtitle(element_blank()) +   
	scale_x_discrete(labels=c("Patrilineal", "Matrilineal")) + 
	ylab("Predicted probability of hypertension") + 
	# scale_y_continuous(limits=c(0, .10)) + 
	scale_color_aaas() + 
	scale_shape_manual(values=c(17, 16, 17, 16)) + 
	theme_classic(base_size=16) + 
	theme(legend.position = "none",
		plot.margin=unit(c(.5,.5,.5,.5), "cm"),
		axis.title.y=element_text(margin=margin(0,15,0,0)),
		axis.title.x=element_blank())

ggarrange(CRP, HTN, 
		  labels = c("A", "B"),
		  ncol = 2, nrow = 1,
		  common.legend = TRUE, legend = "right")

