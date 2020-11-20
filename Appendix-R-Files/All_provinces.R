# This script fits the SEIR model to each province separately and then a hierarchical model to canada as a whole. 

library(rjags)
library(runjags)
library(modeest) 

# Specify the assumed values for infectious period and latent period
w <- 1/5; 
g <- 1/2;

# Give the filename for the output
file.name=c("filename.RData")


####################################
# Read in and format the data for analysis
z1 <- read.csv("cases.csv")
Province.Data <- read.csv("Province_data.csv")
Prov.Name <- c("AB.RData","MB.RData","NB.RData","NL.RData","NS.RData","ON.RData","PE.RData","QC.RData","SK.RData","BC.RData")

#end.day <- as.Date("2020-03-25") # March 25 end date
end.day <- as.Date("2020-04-1"); N.days<-72; # April 1 end date and the total number of days for longest province (ON)

# read in the data and filter by province

N.Prov <- length(Province.Data$Province);
Start.Days <- c()
Province.Data$L_TS<-0;

C<-matrix(-1,nrow=N.days,ncol=10);

for(j in 1:(N.Prov-1)){
	z<-subset(z1,province==as.character(Province.Data$Province[j]))
	z$date <- as.Date(z$date_report,format="%d-%m-%Y")
	start.day<-as.Date(Province.Data$Start_Date[j]) 
	z$day <- julian(z$date,start.day)
	Province.Data$L_TS[j]<-julian(end.day,start.day) 
	for(i in 1:Province.Data$L_TS[j]){C[i,j] <- length(subset(z,day==i)$day);}
	Province.Data$L_T0[j] <- Province.Data$L_TS[j]-2; #Length of the vector of possible spill-over days to search over
}
	j<-10 #to put BC last
	z<-subset(z1,province==as.character(c("BC")))
	z$date <- as.Date(z$date_report,format="%d-%m-%Y")
	start.day<-as.Date(Province.Data$Start_Date[j]) 
	z$day <- julian(z$date,start.day)
	Province.Data$L_TS[j]<-julian(end.day,start.day) 
	for(i in 1:Province.Data$L_TS[j]){C[i,j] <- length(subset(z,day==i)$day);}
	Province.Data$L_T0[j] <- Province.Data$L_TS[j]-2; #Length of the vector of possible spill-over days to search over

p.Tc <- matrix(0,nrow=max(Province.Data$L_TS),ncol=N.Prov)
for(i in 1:N.Prov){p.Tc[2:(Province.Data$L_T0[i]+1),i]<-1/Province.Data$L_T0[i]}

Sundays <- seq(3,66,by=7); # for the BC data

N <- Province.Data$N;
L_TS <- Province.Data$L_TS;

#####################################################

#####################################################
# Controls for Bayesian fitting of the SEIR model


# Arrays to record parameter value results
R0.L <- array(0,dim=c(11,3))
T0.L <- array(0,dim=c(11,3))
psi.L <- array(0,dim=c(11,3))
B.L <- array(0,dim=c(11,3))
theta.L <- array(0,dim=c(11,3))

# Fit the model to the first 9 provinces
for(i in 1:9){

	print(Prov.Name[i])

	# Fit the model
	L_T0<-L_TS[i]-2;
	fit<-run.jags(model="Canada.bug", 	data=list("N"=N[i],"L_TS"=L_TS[i],"L_T0"=L_T0,"C"=C[,i],"w"=w,"g"=g),monitor=c("psi","B","T0"),n.chains=16,method="rjags", burnin=60000, adapt=60000,sample=1000);

	m<-summary(fit)
	psi.L[i,1:3] <- m[1,1:3];
	B.L[i,1:3] <- m[2,1:3];
	T0.L[i,1:3] <- m[3,1:3];
	R0.L[i,1:3] <- B.L[i,1:3]/w;
	theta.L[i,1:3] <- 1 - 1/R0.L[i,1:3];

	save(m, fit, file=Prov.Name[i])		
}

# Fit the model to BC data (must be done separately due to no reporting of cases on Sundays)
#BC
for(i in 10:10){

	# Fit the model
	L_T0<-L_TS[i]-2;
	fit<-run.jags(model="BC.bug", 	data=list("N"=N[i],"L_TS"=L_TS[i],"L_T0"=L_T0,"C"=C[,i],"w"=w,"g"=g,"Sundays"=Sundays),monitor=c("psi","B","T0"),n.chains=16,method="rjags", burnin=60000, adapt=60000,sample=1000);

	m<-summary(fit)
	psi.L[i,1:3] <- m[1,1:3];
	B.L[i,1:3] <- m[2,1:3];
	T0.L[i,1:3] <- m[3,1:3];
	R0.L[i,1:3] <- B.L[i,1:3]/w;
	theta.L[i,1:3] <- 1 - 1/R0.L[i,1:3];

	save(m, fit, file=Prov.Name[i])		
}



#Canada Hierarchical 

fit<-run.jags(model="CA_heirarchical.bug", data=list("N"=N,"L_TS"=L_TS,"C"=C,"N.Prov"=N.Prov,"p.Tc"=p.Tc,"N.days"=N.days,"Sundays"=Sundays,"w"=w,"g"=g),monitor=c("psi.mu","psi.sigma","B.mu","B.sigma","Tc"),n.chains=16,method="rjags", burnin=60000, adapt=60000,sample=1000);

	i<-11;
	m<-summary(fit)
	psi.L[i,1:3] <- exp(m[1,1:3]);
	B.L[i,1:3] <- exp(m[3,1:3]);
	R0.L[i,1:3] <- B.L[i,1:3]/w;
	theta.L[i,1:3] <- 1 - 1/R0.L[i,1:3];
	
	m.CA<-m; fit.CA<-fit;

save(m.CA,fit.CA,R0.L,T0.L,psi.L,B.L,theta.L,file=file.name)




