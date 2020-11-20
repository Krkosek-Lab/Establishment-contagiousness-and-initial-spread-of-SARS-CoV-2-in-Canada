# SEIR Model BUGS

model {
	
	# Parameters and priors

	psi.mu ~ dunif(-5,5) # immigration rate of infection
	psi.sigma ~ dunif(0,1000000) 
	psi.tau <- pow(psi.sigma, -2)
	
	B.mu ~ dunif(-5,5) #transmission coefficient
	B.sigma ~ dunif(0,1000000)
	B.tau <- pow(B.sigma, -2)
	

	# Spillover days
	for(i in 1:N.Prov){Tc[i] ~ dcat(p.Tc[1:N.days,i])}
	
	#Heirarchical parts
	for(i in 1:N.Prov){B.p[i] ~ dnorm(B.mu, B.tau)}
	for(i in 1:N.Prov){psi.p[i] ~ dnorm(psi.mu, psi.tau)}

	# State variables, initial conditions
	for(i in 1:N.Prov){
		T[1,i] <- exp(psi.p[i])/w
		S[1,i] <- N[i]
		E[1,i] <- 0
		I[1,i] <- 0
	}
	for(i in 1:N.Prov){
		T[2,i] <- exp(psi.p[i])/w
		S[2,i] <- N[i]
		E[2,i] <- 0
		I[2,i] <- 0
	}


	# Expectation 
	for(i in 1:N.Prov){y_hat[1,i]<-0} 


# Model and likelihood for all provinces except BC
for(j in 1:9){ 
	for(t in 3:L_TS[j]){
			
		T[t,j] <- T[t-1,j] + exp(psi.p[j]) - w*T[t-1,j]
#		T[t,j] <- T[t-1,j] + psi[j] - w*T[t-1,j]
		S[t,j] <- S[t-1,j] - exp(B.p[j])*S[t-1,j]*I[t-1,j]/N[j];
		E[t,j] <- E[t-1,j] + exp(B.p[j])*S[t-1,j]*I[t-1,j]/N[j] - g*E[t-1,j]
		I[t,j] <- ifelse(t==Tc[j], 1, I[t-1,j] + g*E[t-1,j] - w*I[t-1,j])
	
		y_hat[t,j] <- (T[t,j]+I[t,j]) #Expectation
		C[t,j] ~ dpois(y_hat[t,j]) #Likelihood					
		}		
	}
	
# Model and likelihood for all provinces except BC
for(j in 10:10){
	for(t in 3:L_TS[j]){
			
		T[t,j] <- T[t-1,j] + exp(psi.p[j]) - w*T[t-1,j]
#		T[t,j] <- T[t-1,j] + psi[j] - w*T[t-1,j]
		S[t,j] <- S[t-1,j] - exp(B.p[j])*S[t-1,j]*I[t-1,j]/N[j];
		E[t,j] <- E[t-1,j] + exp(B.p[j])*S[t-1,j]*I[t-1,j]/N[j] - g*E[t-1,j]
		I[t,j] <- ifelse(t==Tc[j], 1, I[t-1,j] + g*E[t-1,j] - w*I[t-1,j])
	
		y_hat[t,j] <- ifelse(min(abs(t-Sundays))==0, 0, ifelse(t==67, (T[t,j]+T[t-1,j]+I[t,j]+I[t-1,j]) ,ifelse(t==61, (T[t,j]+T[t-2,j]+I[t,j]+I[t-2,j]), ifelse(t==54, (T[t,j]+T[t-2,j]+I[t,j]+I[t-2,j]) ,(T[t,j]+I[t,j])))));

		C[t,j] ~ dpois(y_hat[t,j])
	}

}
}