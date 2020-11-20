# SEIR Model BUGS

model {
	
#	w <- 1/7; # adjust this parameter here for each scenario
#	g <- 1/5;	

	psi ~ dunif(0,3) # immigration rate of infection
	B ~ dunif(0,2) #transmission coefficient
	
	# Day of spillover from immigration to community transmission
	for(i in 1:L_T0){p[i]<-1/L_T0}
	T0 ~ dcat(p[])
	Tc <- T0+1

	# State variables, initial conditions
	T[1] <- psi/w
	S[1] <- N
	E[1] <- 0
	I[1] <- 0

	y_hat[1]<-0; #Expectation

	for(t in 2:L_TS){
			
		T[t] <- T[t-1] + psi - w*T[t-1]
		S[t] <- S[t-1] - B*S[t-1]*(I[t-1])/N
		E[t] <- E[t-1] + B*S[t-1]*(I[t-1])/N - g*E[t-1]
		I[t] <- ifelse(t==Tc, 1, I[t-1] + g*E[t-1] - w*I[t-1])
	
		y_hat[t] <- (T[t]+I[t]); #Expectation
		C[t] ~ dpois(y_hat[t]); #Likelihood
					
	}		
}