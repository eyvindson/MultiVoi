rm(list=ls())
library(compiler)
library(mgcv)
library(ggplot2)
library(plyr)
library(ggpubr)
enableJIT(3)

set.seed(20180507)

#PATH WHERE FIGURES WILL BE PLACED
#Identify your personal path here.
path = "U:/My Documents/Inventory_Uncertainty_JYU/"


#DATA
#number of classes:
J <-4
#Conservation costs
b_cost <- c(48000,66000,40200,35600)
#Inventory costs
c_cost <- c(5100,5800,6200,5600)

#data_v : number of species found in class
#data_r : stepwise cumulative distribution function
data_v_1 <-c(6,7,8,9,10,11,12,13,14,15)
data_r_1 <-c(0.142857143,0.285714286,0.5,0.642857143,0.714285714,0.785714286,0.857142857,0.857142857,0.857142857)
data_v_2 <-c(2,5,6,7,8,9,10,11,13,14,15,16)
data_r_2 <-c(0.058823529,0.117647059,0.235294118,0.352941176,0.411764706,0.470588235,0.588235294,0.705882353,0.823529412,0.882352941,0.941176471)
data_v_3 <-c(5,6,7,8,9,10,11,12,13,21)
data_r_3 <-c(0.055555556,0.222222222,0.333333333,0.388888889,0.611111111,0.722222222,0.833333333,0.888888889,0.944444444)
data_v_4 <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16)
data_r_4 <-c(0.047619048,0.047619048,0.19047619,0.238095238,0.333333333,0.428571429,0.523809524,0.571428571,0.619047619,0.714285714,0.80952381,0.857142857,0.904761905,0.952380952)

#SINGLE-STAND PROBLEM (Stand 2): poisson distribution, lambda = 9.6, inventory cost = 5800, conservation cost = 66000

lambda <- 9.6
y0 <- 12
b <- b_cost[2] #Conservation cost
c <- c_cost[2] #Inventory cost

Fy0 <-  ppois(y0,lambda)
k <- 0:(y0-1)
EYcond <- (lambda-sum((k*dpois(k,lambda))))/(1-ppois((y0-1),lambda))

VOI1 <- -(1-Fy0)*(b+c) - Fy0*c
VOI2 <- (1-Fy0)*EYcond 

#Multicriteria value of information: y0=12
cat(VOI1,VOI2,"\n")

y0 <- 7
Fy0 <-  ppois(y0,lambda)
k <- 0:(y0-1)
EYcond <- (sum((k*dpois(k,lambda))))/(ppois((y0-1),lambda))

VOI1 <- -(1-Fy0)*(c) + Fy0*(b-c)
VOI2 <- -Fy0*EYcond 

#Multicriteria value of information: y0=7
cat(VOI1,VOI2,"\n")

#MULTI-STAND PROBLEM

#Fitting the data to a monotonic spline
prob_1 <- diff(c(0,data_r_1,1))
prob_2 <- diff(c(0,data_r_2,1))
prob_3 <- diff(c(0,data_r_3,1))
prob_4 <-diff(c(0,data_r_4,1))
dataobj <- list( list(x=data_v_1, prob=prob_1), list(x=data_v_2, prob=prob_2), 
			 list(x=data_v_3, prob=prob_3), list(x=data_v_4, prob=prob_4))

lambda <- c(sum(data_v_1*prob_1), sum(data_v_2*prob_2), 
		sum(data_v_3*prob_3), sum(data_v_4*prob_4))
datavariance <- c(sum((data_v_1 - lambda[1])^2 * prob_1), sum((data_v_2 - lambda[2])^2 * prob_2), 
			  sum((data_v_3 - lambda[3])^2 * prob_3), sum((data_v_4 - lambda[4])^2 *prob_4))
nbsize <- lambda^2 / (pmax(0.0001,datavariance - lambda))

#The parameters of the negative binomial distribution are estimated by the method of moments because it is simple to implement
#For stand 2 the variance is smaller than the mean and consequently the size parameter becomes negative and the negative binomial 
#distribution is not defined. I have solved this by an adhoc correction (pmax). In practice, this means that the distribution for 
#stand 2 is Poisson(lambda[2]).
#From help:
#An alternative parametrization (often used in ecology) is by the mean mu (see above), and size, the dispersion parameter, 
#where prob = size/(size+mu). The variance is mu + mu^2/size in this parametrization.

#Checking the model fit
offset <- 0.1
for(j in 1:J)
{  
domain <- seq(0,1.5*max(dataobj[[j]]$x))
fittedprob <- dnbinom(domain,mu=lambda[j], size=nbsize[j])
}  

#A list describing the options for conservation decisions or inventory decisions
combilist <- list()
for(j in 1:J)
{  
combilist[[j]] <- c(0,1)
}    
combidf <- expand.grid(combilist)
colnames(combidf) <- paste0("Stand",1:J)
JJ <- nrow(combidf)

#Conservation Decision names
DEC_NAMES <- c("A1","B1","B2","C1","B3","C2","C3","D1","B4","C4","C5","D2","C6","D3","D4","E1")
#Inventory scheme names
INV_NAMES <- c("V1","W1","W2","X1","W3","X2","X3","Y1","W4","X4","X5","Y2","X6","Y3","Y4","Z1")

#Function inner() found from the web
require(dplyr)
inner = function(a, b, FUN=prod) {
apply(a, 1, function(a_col) {
apply(b, 2, function(b_row) {
  sapply(seq_along(a_col), function(i) {
	FUN(a_col[i], b_row[i])
  }) %>%
	sum()
})
}) %>%
t()
}

#An implementation for finding non-dominant solutions.
nondominant <- function(x,returnlist=T)
{
n <- nrow(x)
m <- ncol(x)
larger <- function(z,y)
{
return(as.numeric(z > y))
}
A <- inner(x,t(x),larger) + (inner(x,t(x),identical)==m)
domind <- as.logical(rowSums(A==0))
if(returnlist)
{
return(list(index=!domind,values=x[!domind,]))  
} else {  
return(x[!domind,])
}  
}  

#Creating a combination of all options
planexpand <- function(x,combidf)
{
if(is.data.frame(x)) return( as.matrix(combidf) %*% t(as.matrix(x)) ) 
else return( as.matrix(combidf) %*% cbind(x) )
}  

#Costs as a matrix 
cons_costs = as.matrix(combidf) %*% cbind(-b_cost)
inv_costs = as.matrix(-1*planexpand(c_cost,combidf))

#Calculating the outcomes for each scenario
newsolution <- function(nspecies_known,nspecies_true,nspecies_blue,origsolution_ind,jjr)
{
#input: vectors of length of JJ
nspecies_blue_update <- cbind(bluef(cons_costs+inv_costs[jjr]))
delta <- nspecies_known - nspecies_blue_update 
deltagain <- delta - delta[origsolution_ind]
new_ind <- which.max(deltagain) #Identifying the updated decision
truegain <- nspecies_true[new_ind] - nspecies_blue_update[new_ind]  #Real difference to the indifference curve
costs <- cons_costs[new_ind]+inv_costs[jjr]
return(list(ind = new_ind, deltagain = deltagain[new_ind], truegain = truegain, nspecies = nspecies_true[new_ind],costs = costs))
}  

nspeciesknown <- function(nspecies_true,nspecies_exp,combidf)
{
#input: vectors of length of J
nspecies_true_mat <- cbind(rep(1,nrow(combidf))) %*% rbind(nspecies_true)
nspecies_exp_mat <- cbind(rep(1,nrow(combidf))) %*% rbind(nspecies_exp)
nspecies_known <- combidf * nspecies_true_mat + (1 - combidf) * nspecies_exp_mat 
return(nspecies_known)
}  

#A function to fit a monotonic spline that goes through the indifference points given by the DM
bluecurve <- function(indiffind,indiffvalue,altm,solutionind=NULL,nondomind=NULL)
{
indiff <- cbind(altm[indiffind,1],indiffvalue)
if(!is.null(solutionind)) indiff <- rbind(altm[solutionind,], indiff)
return(splinefun(indiff,method="monoH.FC"))
}  

bsum <- sum(b_cost)
lambdasum <- sum(lambda)

inialt <- cbind(as.matrix(combidf) %*% cbind(-b_cost), as.matrix(combidf) %*% cbind(lambda))
#PREFERENCE INFORMATION
origsolution_ind <- 13 #identifying the most preferred solution
bluef <- bluecurve(c(16,14,6,9,1), c(46,28,19,9,3), altm=inialt, solutionind=origsolution_ind)
indiff_solution <- c(16,14,6,9,1)
indiff_values <-c(46,28,19,9,3)
ininondom <- nondominant(inialt)

xpoints <- sort( c( seq(min(inialt[,1]),0, length.out = 50), inialt[,1]))
blueypoints <- bluef(xpoints)

#Simulation for finding optimal inventory decision
yblue <- cbind(bluef(inialt[,1]))
M <- 10000 #Number of simulation runs (datasets generated)
ygenall <- matrix(NA,nrow=M,ncol=J,byrow=T)
for(j in 1:J)
{
ygenall[,j] <- rnbinom(M,size=nbsize[j],mu=lambda[j]) 
}

kk = 0
results_temp <- data.frame( simu = rep((1+kk*M/20):((kk+1)*M/20), each=JJ), invplan = rep(1:JJ, M/20), newind_1 = NA, deltagain_1 = NA, truegain_1 = NA, newind_2 = NA, deltagain_2 = NA, truegain_2 = NA)
results = NULL
row <- 1
ll = 0

for(m in 1:M)
{  
ll =ll+1
yknown <- nspeciesknown(ygenall[m,],lambda,combidf)
for(jjr in 1:JJ)
{
newsol <- newsolution(planexpand(yknown[jjr,],combidf), planexpand(ygenall[m,],combidf), yblue, origsolution_ind,jjr)

results_temp$newind[row] <- newsol$ind
results_temp$deltagain[row] <- newsol$deltagain
results_temp$truegain[row] <- newsol$truegain
results_temp$nspecies[row] <- newsol$nspecies
results_temp$costs[row] <- newsol$costs
row <- row + 1
}  
#Loop to identify how long the process will take.
if(ll == M/20){
ll = 0
kk = kk+1
results = rbind(results, results_temp)
print(m)
row = 1
results_temp <- data.frame( simu = rep((1+kk*M/20):((kk+1)*M/20), each=JJ), invplan = rep(1:JJ, M/20), newind_1 = NA, deltagain_1 = NA, truegain_1 = NA, newind_2 = NA, deltagain_2 = NA, truegain_2 = NA)
}
}

#Summarizing simulation results
results$change <- results$newind != origsolution_ind   

#ADDING THE EXPECTED NUMBER OF SPECIES IN THE SUMMARY
ressummary <- ddply(results,.(invplan),summarize,truegain = mean(truegain), E_Cost = mean(costs), E_nspecies = mean(nspecies), probchange = mean(change),deltagain = mean(deltagain))
ressummary$invcost <- as.vector(planexpand(c_cost,combidf))
ressummary$invdetails <- apply(combidf,1,paste,collapse="")

#ADDED TO PROVIDE MORE DETAIL
colnames(inialt) <-c("Cost","lambda")

#Added to evaluate the number of conservation solutions.
combns = expand.grid(seq(1,16,1), seq(1,16,1))
combns = cbind (combns, apply(combns, 1, function(x)sum(results$newind==x[2] & results$invplan==x[1])))
colnames(combns) <- c("X","Y","Z")
indepth <- matrix(combns$Z,nrow=16,ncol=16)

colnames(indepth) <- c("P1","P2","P3,","P4,","P5,","P6,","P7","P8,","P9","P10","P11","P12","P13","P14","P15","P16")
rownames(indepth) <- rownames(ressummary)
ressummary <- ressummary[,c("invplan","invdetails","invcost","truegain","probchange","deltagain","E_Cost","E_nspecies")]
ressummary<- cbind(ressummary, indepth) 

print(ressummary)

#EVALUATING THE EXPECTED COSTS, RATHER THAN ONLY THE INVENTORY COSTS
#Plotting the alternatives for the inventory decision
invnondom <- nondominant(cbind(ressummary$E_Cost, ressummary$truegain))
xpointsinv <- sort( c( seq(min(ressummary$E_Cost),max(ressummary$E_Cost), length.out = 50), ressummary$E_Cost))
blueypointsinv <- bluef(xpointsinv) - inialt[origsolution_ind,2]

#Non dominated solutions, comparing expected cost vs true gain
E_costnondom <- nondominant(cbind(ressummary$E_Cost, ressummary$truegain))

#Plotting expected costs vs improvement from indifference curve
E_Cost_vs_Improvement <- ggplot(ressummary, aes(x =(E_Cost-E_Cost[1])/1000, y = truegain)) + 
annotate("text",x=(ressummary$E_Cost[!E_costnondom$index]-ressummary$E_Cost[1])/1000, y=ressummary$truegain[!E_costnondom$index], label=INV_NAMES[1:JJ][!E_costnondom$index], col="grey")+
annotate("text",x=(ressummary$E_Cost[E_costnondom$index]-ressummary$E_Cost[1])/1000, y=ressummary$truegain[E_costnondom$index], label=INV_NAMES[1:JJ][E_costnondom$index], col="black")+theme_bw()+
xlab("Change in expected costs (thousands €)")+ylab("Change from indifference curve")

#PRINTING FIGURES
blue_line <-data.frame(xpoints, blueypoints)

#Stand  data figure
data1<-data.frame(x =dataobj[[1]]$x,prob=dataobj[[1]]$prob,z="Actual")
fittedprob <- dnbinom(domain,mu=lambda[1], size=nbsize[1])
data1<-rbind(data1,data.frame(x=domain,prob=fittedprob,z="Fitted"))
for(j in 0:length(subset(data1,z =="Fitted")$x))
{if(!j %in% subset(data1,z == "Actual")$x)
{data1<- rbind(data1,c(j,"0","Actual"))}
}
data1[1] <- sapply(data1[1], as.numeric)
data1[2] <- sapply(data1[2], as.numeric)
ST1 <-ggplot(data1, aes(x= x, y =prob,fill = z))+geom_bar(stat="identity",position=position_dodge())+theme_bw()+ scale_fill_manual(values=c('Black','Grey'))+ggtitle(c("Age group 1: <=80 years"))+xlab(c("Number of Species"))+ylab(c("Occurance probability"))+labs(fill="")+scale_y_continuous(expand = c(0, 0), limits = c(0,0.23))
ST1<-ST1+annotate("label", x = Inf, y = Inf, label = "Inventoried\n stands: 14",hjust = "inward",vjust = "inward",size = 4)

data2<-data.frame(x =dataobj[[2]]$x,prob=dataobj[[2]]$prob,z="Actual")
fittedprob <- dnbinom(domain,mu=lambda[2], size=nbsize[2])
data2<-rbind(data2,data.frame(x=domain,prob=fittedprob,z="Fitted"))
for(j in 0:length(subset(data2,z =="Fitted")$x))
{if(!j %in% subset(data2,z == "Actual")$x)
{data2<- rbind(data2,c(j,"0","Actual"))}
}
data2[1] <- sapply(data2[1], as.numeric)
data2[2] <- sapply(data2[2], as.numeric)
ST2 <-ggplot(data2, aes(x= x, y =prob,fill = z))+geom_bar(stat="identity",position=position_dodge())+theme_bw()+ scale_fill_manual(values=c('Black','Grey'))+ggtitle(c("Age group 2: 81-95 years"))+xlab(c("Number of Species"))+ylab(c("Occurance probability"))+labs(fill="")+scale_y_continuous(expand = c(0, 0), limits = c(0,0.23))
ST2<-ST2+annotate("label", x = Inf, y = Inf, label = "Inventoried\n stands: 17",hjust = "inward",vjust = "inward",size = 4)

data3<-data.frame(x =dataobj[[3]]$x,prob=dataobj[[3]]$prob,z="Actual")
fittedprob <- dnbinom(domain,mu=lambda[3], size=nbsize[3])
data3<-rbind(data3,data.frame(x=domain,prob=fittedprob,z="Fitted"))
for(j in 0:length(subset(data3,z =="Fitted")$x))
{if(!j %in% subset(data3,z == "Actual")$x)
{data3<- rbind(data3,c(j,"0","Actual"))}
}
data3[1] <- sapply(data3[1], as.numeric)
data3[2] <- sapply(data3[2], as.numeric)
ST3 <-ggplot(data3, aes(x= x, y =prob,fill = z))+geom_bar(stat="identity",position=position_dodge())+theme_bw()+ scale_fill_manual(values=c('Black','Grey'))+ggtitle(c("Age group 3: 96-110 years"))+xlab(c("Number of Species"))+ylab(c("Occurance probability"))+labs(fill="")+scale_y_continuous(expand = c(0, 0), limits = c(0,0.23))
ST3<-ST3+annotate("label", x = Inf, y = Inf, label = "Inventoried\n stands: 18",hjust = "inward",vjust = "inward",size = 4)

data4<-data.frame(x =dataobj[[4]]$x,prob=dataobj[[4]]$prob,z="Actual")
fittedprob <- dnbinom(domain,mu=lambda[4], size=nbsize[4])
data4<-rbind(data4,data.frame(x=domain,prob=fittedprob,z="Fitted"))
for(j in 0:length(subset(data4,z =="Fitted")$x))
{if(!j %in% subset(data4,z == "Actual")$x)
{data4<- rbind(data4,c(j,"0","Actual"))}
}
data4[1] <- sapply(data4[1], as.numeric)
data4[2] <- sapply(data4[2], as.numeric)
ST4 <-ggplot(data4, aes(x= x, y =prob,fill = z))+geom_bar(stat="identity",position=position_dodge())+theme_bw()+ scale_fill_manual(values=c('Black','Grey'))+ggtitle(c("Age group 4: > 110 years"))+xlab(c("Number of Species"))+ylab(c("Occurance probability"))+labs(fill="")+scale_y_continuous(expand = c(0, 0), limits = c(0,0.23))
ST4<-ST4+annotate("label", x = Inf, y = Inf, label = "Inventoried\n stands: 21",hjust = "inward",vjust = "inward",size = 4)

#Creating Figure 1: Stand distributions
postscript(file=paste(path,"Figure_1.eps",sep=""),horizontal=FALSE, onefile=FALSE,width=9,height=9,pointsize=12)
ggarrange(ST1, ST2,ST3,ST4, ncol = 2,nrow = 2, common.legend = TRUE, legend= "right")
dev.off()

#Creating Figure 2
postscript(file=paste(path,"Figure_2.eps",sep=""),horizontal=FALSE, onefile=FALSE,width=6,height=6,pointsize=12)
xpoints <- sort( c( seq(min(inialt[,1]),0, length.out = 50), inialt[,1]))
blueypoints <- bluef(xpoints)
plot(inialt[,1]/1000, inialt[,2], type="n", ylim=c(0, max(blueypoints)), 
 xlab = "Conservation cost (thousands €)", ylab = "Expected number of species")
text(inialt[!ininondom$index,1]/1000, inialt[!ininondom$index,2], DEC_NAMES[1:JJ][!ininondom$index], col="grey")
text(inialt[ininondom$index,1]/1000, inialt[ininondom$index,2], DEC_NAMES[1:JJ][ininondom$index])
ord <- order(inialt[ininondom$index,1])
lines(inialt[ininondom$index,1][ord]/1000, inialt[ininondom$index,2][ord], type="S")
lines(xpoints/1000, blueypoints, col="black", lty = 2)
points(inialt[indiff_solution]/1000, indiff_values, col="black")
points(inialt[origsolution_ind]/1000,  inialt[origsolution_ind,][2], col="black")
text(inialt[origsolution_ind,1]/1000, inialt[origsolution_ind,2], DEC_NAMES[1:JJ][origsolution_ind],  font=2)
legend("topright", legend=c("Pareto front","Indifference curve"), lty=c(1,2),col=c("black","black"))
dev.off()

#Figure 3 : Expected conservation benefits and inventory costs
postscript(file=paste(path,"Figure_3.eps",sep=""),horizontal=FALSE, onefile=FALSE,width=6,height=6,pointsize=12)
E_Cost_vs_Improvement
dev.off()

#Creating Figure 4
postscript(file=paste(path,"Figure_4.eps",sep=""),horizontal=FALSE, onefile=FALSE,width=8,height=6,pointsize=10)
par(mfrow=c(1,3))  
pldata = subset(results,invplan == 1)
pldata$costs <- pldata$costs/1000
proportion = table(pldata$costs)
proportion  <- proportion/1000
boxplot(pldata$nspecies ~ pldata$costs,boxwex = 10, width = proportion ,at = sort(unique(pldata$costs)), xlim=c(-200, 0),ylim=c(0,50),xaxt="n",lty="solid",
main = "Scheme V1",  xlab = "Cost (thousands)", ylab = "Number of species")
axis(1, at = seq(0,-200,-20) , labels = seq(0,-200,-20) )
lines(blue_line$xpoints/1000, blue_line$blueypoints, col="black", lty=2)

pldata = subset(results,invplan == 10)
pldata$costs <- pldata$costs/1000
proportion = table(pldata$costs)
proportion  <- proportion/1000
boxplot(pldata$nspecies ~ pldata$costs,boxwex = 5, width = proportion ,at = sort(unique(pldata$costs)), xlim=c(-200, 0),ylim=c(0,50),xaxt="n",lty="solid",
main = "Scheme X4",  xlab = "Cost (thousands)", ylab = "Number of species")
axis(1, at = seq(0,-200,-20) , labels = seq(0,-200,-20) )
lines(blue_line$xpoints/1000, blue_line$blueypoints, col="black", lty=2)

pldata = subset(results,invplan == 9)
pldata$costs <- pldata$costs/1000
proportion = table(pldata$costs)
 proportion  <- proportion/1000
boxplot(pldata$nspecies ~ pldata$costs,boxwex = 5, width = proportion ,at = sort(unique(pldata$costs)), xlim=c(-200, 0),ylim=c(0,50),xaxt="n",lty="solid",
main = "Scheme W4",  xlab = "Cost (thousands)", ylab = "Number of species")
axis(1, at = seq(0,-200,-20) , labels = seq(0,-200,-20) )
lines(blue_line$xpoints/1000, blue_line$blueypoints, col="black", lty=2)
dev.off()
