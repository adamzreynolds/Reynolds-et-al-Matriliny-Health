

## ANALYSES FOR:
# Reynolds AZ, Wander K, Sum CY, Su M, Emery Thompson M, Hooper PL, Li H, Shenk MK, Starkweather KE, Blumenfield T, Mattison SM. (2020) Matriliny reverses gender disparities in inflammation and hypertension among the Mosuo of China.

#-------------------------------------------------
# READ DATA FROM FILE

data <- read.csv("Dataset_S1.csv")



#-------------------------------------------------
# REPORT SAMPLE SIZES AND PREVALENCES

# N, Base Models
nrow(data[!is.na(data$htn),])
nrow(data[!is.na(data$elevatedCRP),])

# N, HoH Models
nrow(data[!is.na(data$htn) & !is.na(data$Head),])
nrow(data[!is.na(data$elevatedCRP) & !is.na(data$Head),])

# N, BMI Models
nrow(data[!is.na(data$htn) & !is.na(data$Head) & !is.na(data$OverObese) & !is.na(data$Underweight),])
nrow(data[!is.na(data$elevatedCRP) & !is.na(data$Head) & !is.na(data$OverObese) & !is.na(data$Underweight),])

# TABLES OF CRP PREVALENCES
data_mat <- data[data$Matriliny==1,]
data_pat <- data[data$Matriliny==0,]
mytable <- ftable(xtabs(~elevatedCRP+Male, data=data_mat))
prop.table(mytable, 2) * 100
mytable <- ftable(xtabs(~elevatedCRP+Male, data=data_pat))
prop.table(mytable, 2) * 100

# TABLES OF HTN PREVALENCES
mytable <- ftable(xtabs(~htn+Male, data=data_mat))
prop.table(mytable, 2) * 100
mytable <- ftable(xtabs(~htn+Male, data=data_pat))
prop.table(mytable, 2) * 100


#-------------------------------------------------
# BAYESIAN LOGISTIC REGRESSION MODELS 

require(rethinking)
# set.seed(39485)

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
# precis(mCRP.patF, prob=0.95)
post.patF <- extract.samples(mCRP.patF, 1e4)
precis(post.patF, prob=.95)
# dens(post.patF)

# what is the posterior prob that the age param is < 0?
sum(post.patF$b.Age < 0) / 1e4

# what is the posterior prob that Male param is > 0
sum(post.patF$b.Male > 0) / 1e4

# what is the posterior probability that the Mat parameter is > 0?
sum(post.patF$b.Matriliny > 0) / 1e4 

# what is the posterior prob that the Male-Matriliny param is < 0?
sum(post.patF$b.Male.Matriliny < 0) / 1e4




#---------------------------------
# CRP HOH MODEL, REF. PAT FEMALES

## SUBSET DATA
d3 <- d2[!is.na(d2$Head),]

start <- list(
	a = mean(d2$elevatedCRP),	# begin with same value as base model
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
# precis(mCRP.hoh, prob=0.95)
post.hoh <- extract.samples(mCRP.hoh, 1e4)
precis(post.hoh, prob=.95)
# dens(post.hoh)

# What is the posterior probability that the Age param is < 0?
sum(post.hoh$b.Age < 0) / 1e4

# What is the posterior probability that the Male param is > 0
sum(post.hoh$b.Male > 0) / 1e4

# What is the posterior probability that the Matriliiny param is > 0?
sum(post.hoh$b.Matriliny > 0) / 1e4

# What is the posterior probability that the Male-Matriliny param is < 0?
sum(post.hoh$b.Male.Matriliny < 0) / 1e4

# What is the posterior probability that the Head param is > 0?
sum(post.hoh$b.Head > 0) / 1e4







#----------------------------------
# HTN MAIN MODEL, REF: PAT FEMALES
dataHTN <- data[!is.na(data$htn) & !is.na(data$Male) & !is.na(data$Matriliny) & !is.na(data$Age),]

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
# precis(mHTN, prob=0.95)
post.htn <- extract.samples(mHTN, 1e4)
precis(post.htn, prob=.95)

# what is the posterior prob that the Age param is < 0?
sum(post.htn$b.Age < 0) / 1e4

# what is the posterior prob that male param is > 0
sum(post.htn$b.Male > 0) / 1e4

# what is the posterior probability that the Matriliny param is > zero?
sum(post.htn$b.Matriliny >  0) / 1e4 

# what is the posterior prob that the Male-Mat param is < 0?
sum(post.htn$b.Male.Matriliny < 0) / 1e4



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
# precis(mHTN.head4, prob=0.95)
post4 <- extract.samples(mHTN.head4, 1e4)
precis(post4, prob=.95)

# What is the posterior prob that the age param is < 0?
sum(post4$b.Age < 0) / 1e4

# What is the posterior prob that the Male param is > 0?
sum(post4$b.Male > 0) / 1e4

# What is the posterior prob that the matriliny param is > 0?
sum(post4$b.Matriliny > 0) / 1e4

# What is the posterior prob that the male-matriliny param is < 0?
sum(post4$b.Male.Matriliny < 0) / 1e4

# What is the posterior prob that the head param is < 0?
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
	b.Underweight = 0,
	b.OverObese = 0
)
flist <- alist(
	elevatedCRP ~ dbinom(1,p),
	logit(p) <- a + b.Male*Male + b.Matriliny*Matriliny + b.Male.Matriliny*Male*Matriliny + b.Age*Age + b.Underweight*Underweight + b.OverObese*OverObese,
	a ~ dnorm(-10,10),
	b.Age ~ dnorm(0,10),
	b.Male ~ dnorm(0,10),
	b.Matriliny ~ dnorm(0,10),
	b.Male.Matriliny ~ dnorm(0,10),
	b.Underweight ~ dnorm(0,10),
	b.OverObese ~ dnorm(0,10)
)
mCRP.met <- map(flist, data=d6, start=start)
# precis(mCRP.met, prob=0.95)
post.met <- extract.samples(mCRP.met, 1e4)
precis(post.met, prob=0.95)

# what is the posterior prob that the Age param is < 0?
sum(post.met$b.Age < 0) / 1e4

# what is the posterior prob that male param is > 0
sum(post.met$b.Male > 0) / 1e4

# what is the posterior probability that the mat parameter is > 0?
sum(post.met$b.Matriliny > 0) / 1e4 

# what is the posterior prob that the male-mat param is < 0?
sum(post.met$b.Male.Matriliny < 0) / 1e4

# what is the posterior prob that the Underweight param is < 0?
sum(post.met$b.Underweight < 0) / 1e4

# what is the posterior prob that the OverObese param is < 0?
sum(post.met$b.OverObese < 0) / 1e4



#-------------------------------
# HTN BMI MODEL, REF: PAT FEMALE
 
dataHTN6 <- dataHTN[!is.na(dataHTN$OverObese) & !is.na(dataHTN$Underweight),]

start <- list(
	a = mean(dataHTN6$htn),
	b.Age = 0,
	b.Male = 0,
	b.Matriliny = 0,
	b.Male.Matriliny = 0,
	b.Underweight = 0,
	b.OverObese = 0
)
flist <- alist(
	htn ~ dbinom(1,p),
	logit(p) <- a + b.Male*Male + b.Matriliny*Matriliny + b.Male.Matriliny*Male*Matriliny + b.Age*Age + b.Underweight*Underweight + b.OverObese*OverObese,
	a ~ dnorm(-10,10),
	b.Age ~ dnorm(0,10),
	b.Male ~ dnorm(0,10),
	b.Matriliny ~ dnorm(0,10),
	b.Male.Matriliny ~ dnorm(0,10),
	b.Underweight ~ dnorm(0,10),
	b.OverObese ~ dnorm(0,10)
)
mHTN.bmi <- map(flist, data=dataHTN6, start=start)
# precis(mHTN.bmi, prob=0.95)
post.htn.bmi <- extract.samples(mHTN.bmi, 1e4)
precis(post.htn.bmi, prob=.95)


# what is the posterior prob that the Age param is < 0?
sum(post.htn.bmi$b.Age < 0) / 1e4

# what is the posterior prob that male param is > 0
sum(post.htn.bmi$b.Male > 0) / 1e4

# what is the posterior probability that the mat parameter is > 0?
sum(post.htn.bmi$b.Matriliny >  0) / 1e4 

# what is the posterior prob that the male-mat param is < 0?
sum(post.htn.bmi$b.Male.Matriliny < 0) / 1e4

# what is the posterior prob that the Underweight param is > 0?
sum(post.htn.bmi$b.Underweight > 0) / 1e4

# what is the posterior prob that the OverObese param is < 0?
sum(post.htn.bmi$b.OverObese < 0) / 1e4








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

# dev.new(height=5, width=11)

CRP <- ggplot(preds, aes(x=Matriliny, y=mean, colour=Gender, group=Gender)) +
	geom_errorbar(aes(ymin=lower, ymax=upper), size=.9, width=.2, position=pd) +
	geom_line(position=pd, size=.9) +
	geom_point(position=pd, size=4, aes(shape=Gender)) + 
	ggtitle(element_blank()) + 
	scale_x_discrete(labels=c("Patrilineal", "Matrilineal")) + 
	ylab("Predicted probability of CRP 3-5 mg/L") + 
	# scale_y_continuous(limits=c(0, .10)) + 
	scale_color_aaas() +
	scale_shape_manual(values=c(17, 16, 17, 16)) + 
	theme_classic(base_size=16) +  
	theme(legend.position = "none",
		legend.title=element_text(color="black", size=20),
		legend.text=element_text(size=20),
		plot.margin=unit(c(.5,.5,.5,.5), "cm"),
		axis.title.y=element_text(size=rel(1.3), margin=margin(0,15,0,0)),
		axis.title.x=element_blank(),
		axis.text.x=element_text(color="black", size=20),
		axis.text.y=element_text(color="black", size=20))
		
HTN <- ggplot(preds.htn, aes(x=Matriliny, y=mean, colour=Gender, group=Gender)) +
	geom_errorbar(aes(ymin=lower, ymax=upper), size=.9, width=.2, position=pd) +
	geom_line(position=pd, size=.9) +
	geom_point(position=pd, size=4, aes(shape=Gender)) + 
	ggtitle(element_blank()) +   
	scale_x_discrete(labels=c("Patrilineal", "Matrilineal")) + 
	ylab("Predicted probability of hypertension") + 
	# scale_y_continuous(limits=c(0.10, .35)) + 
	scale_color_aaas() + 
	scale_shape_manual(values=c(17, 16, 17, 16)) + 
	theme_classic(base_size=16) + 
	theme(legend.position = "none",
		plot.margin=unit(c(.5,.5,.5,.5), "cm"),
		axis.title.y=element_text(size=rel(1.3), margin=margin(0,15,0,0)),
		axis.title.x=element_blank(),
		axis.text.x=element_text(color="black", size=20),
		axis.text.y=element_text(color="black", size=20))

pdf("Figure 1 Revised.pdf", height=6, width=13, onefile=FALSE)
# dev.new(width=11.5, height=5.5)
ggarrange(CRP, HTN, 
		  labels = c("A", "B"),
		  font.label = list(size=20),
		  ncol = 2, nrow = 1,
		  common.legend = TRUE, legend = "right")
dev.off()


# Print Predicted Probabilities from CRP Model
print(preds)

# Print Predicted Probabilities from HRN Model
print(preds.htn)




#-------------------------------------------------
## MODEL W PAT MALES AS REFERENCE CATEGORY:

# # CRP MAIN MODEL, REF: PAT MALES

# d2 <- data[!is.na(data$elevatedCRP) & !is.na(data$Male) & !is.na(data$Matriliny) & !is.na(data$Age),]

# d2$Female <- (d2$Male - 1) * -1

# # CRP MAIN MODEL, REF: PAT FEMALES
# start <- list(
	# a = mean(d2$elevatedCRP),
	# b.Age = 0,
	# b.Female = 0,
	# b.Matriliny = 0,
	# b.Female.Matriliny = 0
# )
# flist <- alist(
	# elevatedCRP ~ dbinom(1,p),
	# logit(p) <- a + b.Female*Female + b.Matriliny*Matriliny + b.Female.Matriliny*Female*Matriliny + b.Age*Age,
	# a ~ dnorm(-10,10),
	# b.Age ~ dnorm(0,10),
	# b.Female ~ dnorm(0,10),
	# b.Matriliny ~ dnorm(0,10),
	# b.Female.Matriliny ~ dnorm(0,10)
# )
# mCRP.matF <- map(flist, data=d2, start=start)
# # precis(mCRP.matF, prob=0.95)
# post.matF <- extract.samples(mCRP.matF, 1e4)
# precis(post.matF, prob=.95)
# # dens(post.matF)

# # what is the posterior prob that the age param is < 0?
# sum(post.matF$b.Age < 0) / 1e4

# # what is the posterior prob that Female param is > 0
# sum(post.matF$b.Female > 0) / 1e4

# # what is the posterior probability that the Mat parameter is > 0?
# sum(post.matF$b.Matriliny > 0) / 1e4 

# # what is the posterior prob that the Female-Matriliny param is < 0?
# sum(post.matF$b.Female.Matriliny < 0) / 1e4


# # HTN MAIN MODEL, REF: PAT MALES

# dataHTN <- data[!is.na(data$htn) & !is.na(data$Male) & !is.na(data$Matriliny) & !is.na(data$Age),]

# dataHTN$Female <- (dataHTN$Male - 1) * -1

# start <- list(
	# a = mean(dataHTN$htn),	
	# b.Age = 0,
	# b.Female = 0,
	# b.Matriliny = 0,
	# b.Female.Matriliny = 0
# )
# flist <- alist(
	# htn ~ dbinom(1,p),
	# logit(p) <- a + b.Age*Age + b.Female*Female + b.Matriliny*Matriliny + b.Female.Matriliny*Female*Matriliny,
	# a ~ dnorm(-10,10),
	# b.Age ~ dnorm(0,10),
	# b.Female ~ dnorm(0,10),
	# b.Matriliny ~ dnorm(0,10),
	# b.Female.Matriliny ~ dnorm(0,10)
# )
# mHTN <- map(flist, data=dataHTN, start=start)
# # precis(mHTN, prob=0.95)
# post.htn <- extract.samples(mHTN, 1e4)
# precis(post.htn, prob=.95)


# # what is the posterior prob that the Age param is < 0?
# sum(post.htn$b.Age < 0) / 1e4

# # what is the posterior prob that the female param is < 0
# sum(post.htn$b.Female >< 0) / 1e4

# # what is the posterior probability that the matriliny parameter is < 0?
# sum(post.htn$b.Matriliny < 0) / 1e4 

# # what is the posterior prob that the male-pat param is > 0?
# sum(post.htn$b.Female.Matriliny > 0) / 1e4





#-------------------------------------------------------
# PREVALENCES FOR TABLE S2

# MAKE AGE CATEGORIES
data$Cohort <- NA
data$Cohort[data$Age < 20] <- "16-19"
data$Cohort[data$Age >= 20 & data$Age < 25] <- "20-24"
data$Cohort[data$Age >= 25 & data$Age < 30] <- "25-29"
data$Cohort[data$Age >= 30 & data$Age < 35] <- "30-34"
data$Cohort[data$Age >= 35 & data$Age < 40] <- "35-39"
data$Cohort[data$Age >= 40 & data$Age < 45] <- "40-44"
data$Cohort[data$Age >= 45 & data$Age < 50] <- "45-49"
data$Cohort[data$Age >= 50 & data$Age < 55] <- "50-54"
data$Cohort[data$Age >= 55 & data$Age < 60] <- "55-59"
data$Cohort[data$Age >= 60 & data$Age < 65] <- "60-64"
data$Cohort[data$Age >= 65 & data$Age < 70] <- "65-69"
data$Cohort[data$Age >= 70] <- "70+"


## CRP COLUMN

# Entire Sample
mean(data$elevatedCRP, na.rm=TRUE) * 100

# Kinship System
mytable <- ftable(xtabs(~Matriliny+elevatedCRP, data=data))
prop.table(mytable, 1) * 100

# Gender
mytable <- ftable(xtabs(~Male+elevatedCRP, data=data))
prop.table(mytable, 1) * 100

# Head
mytable <- ftable(xtabs(~Head+elevatedCRP, data=data))
prop.table(mytable, 1) * 100

# Age Cohort
mytable <- ftable(xtabs(~Cohort+elevatedCRP, data=data))
prop.table(mytable, 1) * 100


## HTN COLUMN

# Entire Sample
mean(data$htn, na.rm=TRUE) * 100

# Kinship System
mytable <- ftable(xtabs(~Matriliny+htn, data=data))
prop.table(mytable, 1) * 100

# Gender
mytable <- ftable(xtabs(~Male+htn, data=data))
prop.table(mytable, 1) * 100

# Head
mytable <- ftable(xtabs(~Head+htn, data=data))
prop.table(mytable, 1) * 100

# Age Cohort
mytable <- ftable(xtabs(~Cohort+htn, data=data))
prop.table(mytable, 1) * 100




