********************************************************************************
************ Summer School - Structural Micro (Practical Sessions) *************
********************************************************************************

********************************************************************************
******************* SESSION 1 & 2 - Solving Value Functions ********************

* The following script contains the solution of bike maintanance problem 
* discussed during practical sessions. 
* It includes: 

*     1. Value Function Iteration, 
*     2. Construction of likelihood funcition
*     3. Estimation of A=5 structural model with infinite horizon
*     4. Estimation of T=4 structural model with finite horizon


* Prepared by Jacek Barszczewski (BSE and UAB)
* Contact: jacek.barszczewski[at]bse.eu

********************************************************************************


**************************************
*** Administrative Commands

cap log close                       
set more off                        
clear all                                 

*** set path to the folder with data
cd "C:\Users\Sergi\Dropbox\PhD\Teaching\Third Year\Structural Micro BSS\Structural Micro BSS\3. Materials"

log using Session1.log, replace text

********************************************************************************
*** PART I: LOAD THE DATA - INFINITE HORIZON

*** Import from excel file

import excel "PracticalSession1_2_data_inf.xls", sheet("data") clear

*** Rename the variables

*bike identifier, age of bike battery in each time period
rename (A B C D E F G) (id t0 t1 t2 t3 t4 t5)    

********************************************************************************
*** PART II: PREPARE THE DATA - INFINITE HORIZON

*** Define matrix of state variables (age of bike battery)

mata st_view(a=., ., "t0 t1 t2 t3 t4 t5")

*** Define number of states (T = inf)

mata A = max(a)

*** Number of observation periods 

mata P = cols(a)

*** Number of observations (bikes) in our data

mata nobs=rows(a)

*** Define decision matrix d with dimensions (nobs x (P - 1))

mata d = J(nobs, P-1, 0)

*** Discrete choice (1 if battery replace, 0 if do not replace - at each time period)
*** Last period is skipped, since we only know state variable, but not decision

mata:
for (i=1; i<=P-1; i++) {
	for (j=1; j<=nobs; j++){	
		if (a[j,i+1]==1){
			d[j,i]=1
		}	
	} 
}
end

********************************************************************************
*** PART III: SOLVING VALUE FUNCTIONS - FIXED POINT ALGORITHM


*** Discount factor
mata beta=0.95

*** Theta initial values (parameters to estimate)
mata theta=J(1,3,0)  
  
* Elements:
* [1,1] - replacement cost (theta_R)
* [1,2] - maintenance cost (theta_M1)
* [1,3] - maintenance cost (theta_M2)

***Initial guess
mata EV0=J(A,1,0)

mata:
function fixed_point(theta, beta, EV0){
	
	/* Input
	* theta - vector of parameters 
	* beta - discount factor 
	* EV0 - initial guess 
	*/
	
	/* Output
	EV1 - Emax matrix 
	*/
	
	EV1=EV0 /* Outcome matrix */
	toler=1e-10 /* Tolerance level in value function iteration */
	rule=1 /* Convergence criteria */
	
	while (rule>toler){
		
		exp11 = exp(-theta[1,1] + beta*EV0[1,1])                    /*exp{v1}*/
		exp01 = exp(-theta[1,2] - theta[1,3]+ beta*(EV0[2,1]))          /*exp{v0(1)}*/
		exp02 = exp(-theta[1,2]*2 - theta[1,3]*4 + beta*EV0[3,1])     /*exp{v0(2)}*/
		exp03 = exp(-theta[1,2]*3 - theta[1,3]*9 + beta*EV0[4,1])     /*exp{v0(3)}*/
		exp04 = exp(-theta[1,2]*4 - theta[1,3]*16 + beta*EV0[5,1])    /*exp{v0(4)} -- EV0(5)==v1 */
		exp05 = exp(-theta[1,2]*5 - theta[1,3]*25 + beta*EV0[5,1]) 

		EV1[1,1]  = log(exp11 + exp01)
		EV1[2,1]  = log(exp11 + exp02)
		EV1[3,1]  = log(exp11 + exp03)
		EV1[4,1]  = log(exp11 + exp04)
		EV1[5,1]  = log(exp11 + exp05)
	   
		rule=norm(EV1-EV0) /* convergence criteria */

		EV0=EV1
		}
	return(EV1)	
	}	
	
end	

*** Value function for initial guess of theta
mata EV1 = fixed_point(theta, beta, EV0)
mata EV1

********************************************************************************
*** PART IV: ESTIMATE THE INFINITE HORIZON MODEL

*** Construct log-likelihood function	
mata: 
function loglike(todo, theta, beta, d, a, EV0, y, g, H){
	
	/* Input 
	* todo - a scalar message variable from optimize() that indicates where to compute 1st and 2nd order derivatives;
	* theta - parameters which are going to be estimated
	* beta - discount factor
	* d - matrix of decisions
	* a - matrix of state variable values
	* EV0 - initial matrix to find a fixed point of Emax 
	*/
	
	/* Output
	* y - the value of the objective function that is being optimized
	* g - the gradient vector
	* H - the Hessian matrix
	*/
	
	/* find fixed point for Emax */
	EV1 = fixed_point(theta, beta, EV0)
	
	/* construct conditional value functions */
	v1 = -theta[1,1] + beta*EV1[1,1] 
	v01 = -theta[1,2] - theta[1,3]+ beta*EV1[2,1]
	v02 = -theta[1,2]*2 - theta[1,3]*4 + beta*EV1[3,1]
	v03 = -theta[1,2]*3 - theta[1,3]*9 + beta*EV1[4,1]
	v04 = -theta[1,2]*4 - theta[1,3]*16 + beta*EV1[5,1]
	v05 = -theta[1,2]*5 - theta[1,3]*25 + beta*EV1[5,1]
	
	/* paysoff diffrences: base category no replacement */ 
	payoffdiff1 = v1 - v01 
	payoffdiff2 = v1 - v02 
	payoffdiff3 = v1 - v03 
	payoffdiff4 = v1 - v04 
	payoffdiff5 = v1 - v05
	payoffdiff = (payoffdiff1 \ payoffdiff2 \ payoffdiff3 \ payoffdiff4\payoffdiff5)
	
	/* using conditional choice probabilities construct log-likelihood function */
	f = J(rows(d),cols(d),0)

	for (i=1; i<=rows(d); i++) {
		for (j=1; j<=cols(d); j++){
				f[i,j] = d[i,j]*log(1 - 1/(1+exp(payoffdiff[a[i,j],1]))) + (1-d[i,j])*log(1/(1+exp(payoffdiff[a[i,j],1])))
		}
	}
	
	y = sum(f) 
	}

end

*** Call your optimization function
mata s=optimize_init()

*** Tell it where to find the evaluator function
mata optimize_init_evaluator(s,&loglike())

*** Define type of evaluator - default is d0 (scalar)
mata optimize_init_evaluatortype(s,"d0")

/*	Remaining evaluator types are e.g.:
	- "d1" - same as "d0" and returns gradient rowvector
	- "d2" - same as "d1" and returns Hessian matrix
	- "gf0" - function() returns colvector value
	- "gf1" - same as "gf0" and returns score matrix
	- "gf2" - same as "gf1" and returns Hessian matrix   */

	
*** Define starting parameters 
mata optimize_init_params(s,J(1,3,0))

*** Indicate to optimizer which variables include as remaining function arguments
mata optimize_init_argument(s,1,beta)
mata optimize_init_argument(s,2,d)
mata optimize_init_argument(s,3,a)
mata optimize_init_argument(s,4,EV1)

*** Choose optimization method - default is Newton-Rhapson ("nr")
mata optimize_init_technique(s,"nr")

*** Returns the value of theta which maximizes loglike(theta).
mata thetahat = optimize(s)
mata thetahat
*** Standard Errors

mata hessian = optimize_result_Hessian(s)

mata hessian=-hessian

*** Take the inverse of a matrix

mata H = luinv(hessian)

*** Compute Standard Errors

mata SE = sqrt(diagonal(H))

*** Display estimated coefficents and their standard errors

mata thetahat \ SE'

********************************************************************************
*** PART V: LOAD THE DATA - FORCED CHOICE AFTER 3 PERIODS

*** Import from excel file
import excel "PracticalSession1_2_data_fin.xls", sheet("data") clear

*** Rename the variables

*bike identifier, age the bike is replaced
rename (A B) (id dur)    

********************************************************************************
*** PART VI: PREPARE THE DATA - FORCED CHOICE AFTER 4 PERIODS

*** Define matrix of state variables (age of a bike)

mata st_view(dur=., ., "dur")

*** Define time horizon

mata T = max(dur)

*** Number of observations (bikes) in our data

mata nobs=rows(a)

*** Create state variable matrix
mata

a = J(nobs,T,0)

for (i=1; i<=T; i++) {
	for (j=1; j<=nobs; j++){	
		if (i <= dur[j,1]) {
			a[j,i] = i
		}
	}
}	

end

mata a

*** Discount factor
mata beta=0.95

*** Theta initial values (parameters to estimate)
mata theta=J(1,3,0) 

*** Define decision matrix d with dimensions (nobs x T)

mata d1 = J(nobs, T, 0)

*** Discrete choice (1 if bike replace, 0 if do not replace - at each time period)

mata:

for (i=1; i<=T; i++) {
	for (j=1; j<=nobs; j++){	
		
		if (i==dur[j,1]){
				d1[j,i] = 1
		}
	} 
}

end

*** Define auxiliary matrix selecting periods contributing to the loglikelihood
mata
D1 = a :> 0
end	

********************************************************************************
*** PART VII: ESTIMATE THE FINITE HORIZON MODEL
cap mata:mata drop s1
cap mata:mata drop loglike_fin()
mata: 

function loglike_fin(todo,theta, beta, d1, D1, y, g, H){

	/* Input 
	* todo - a scalar message variable from optimize() that indicates where to compute 1st and 2nd order derivatives;
	* theta - parameters which are going to be estimated
	* beta - discount factor
	* d1 - matrix of decisions
	* D1 - selection to likelihood
	*/
	
	/* Output
	* y - the value of the objective function that is being optimized
	* g - the gradient vector
	* H - the Hessian matrix
	*/
	
	/* construct conditional value functions */
	v1 = -theta[1,1] +  beta*10
	v03 = -theta[1,2]*3 - theta[1,3]*9 + beta*log(exp(v1))
	v02 = -theta[1,2]*2 - theta[1,3]*4 + beta*log(exp(v03) + exp(v1))
	v01 = -theta[1,2] - theta[1,3] + beta*log(exp(v02) + exp(v1))
	
	/* paysoff diffrences: base category no replacement */ 
	payoffdiff1= v1 - v01
	payoffdiff2= v1 - v02
	payoffdiff3= v1 - v03
	
	/* using conditional choice probabilities construct log-likelihood function */
	f1_fin=d1[.,1]:*log((1-1:/(1+exp(payoffdiff1))))+(1:-d1[.,1]):*log(1:/(1+exp(payoffdiff1)))
	f2_fin=d1[.,2]:*log((1-1:/(1+exp(payoffdiff2))))+(1:-d1[.,2]):*log(1:/(1+exp(payoffdiff2)))
	f3_fin=d1[.,3]:*log((1-1:/(1+exp(payoffdiff3))))+(1:-d1[.,3]):*log(1:/(1+exp(payoffdiff3)))
	
	f_fin = (f1_fin, f2_fin, f3_fin)
	
	y = sum(f_fin:*D1[.,1::3])
}
end


*** Call your optimization function
mata s1=optimize_init()

*** Tell it where to find the evaluator function
mata optimize_init_evaluator(s1,&loglike_fin())

*** Define starting parameters 
mata optimize_init_params(s1,J(1,3,0))

*** Indicate to optimizer which variables include as remaining function arguments
mata optimize_init_argument(s1,1,beta)
mata optimize_init_argument(s1,2,d1)
mata optimize_init_argument(s1,3,D1)

*** Choose optimization method - default is Newton-Rhapson ("nr")
mata optimize_init_technique(s1,"nr")

*** Returns the value of theta which maximizes loglike(theta).
mata thetahat_fin = optimize(s1)

*** Standard Errors

mata hessian_fin = optimize_result_Hessian(s1)

mata hessian_fin=-hessian_fin

*** Take the inverse of a matrix

mata H_fin = luinv(hessian_fin)

*** Compute Standard Errors
mata SE_fin = sqrt(diagonal(H_fin))

*** Display estimated coefficents and their standard errors
mata thetahat_fin \ SE_fin'

********************************************************************************
*** PART VIII: LOAD THE DATA - FINITE HORIZON
mata mata clear
*** Import from excel file

import excel "PracticalSession1_2_data_inf.xls", sheet("data") clear

*** Rename the variables

*bike identifier, age of bike battery in each time period
rename (A B C D E F G) (id t0 t1 t2 t3 t4 t5)    

********************************************************************************
*** PART IX: PREPARE THE DATA - FINITE HORIZON

*** Define matrix of state variables (age of bike battery)

mata st_view(a=., ., "t0 t1 t2 t3 t4 t5")

*** Define number of states (T = inf)

mata A = max(a)

*** Number of observation periods 

mata P = cols(a)

*** Number of observations (bikes) in our data

mata nobs=rows(a)

*** Define decision matrix d with dimensions (nobs x (P - 1))

mata d = J(nobs, P-1, 0)

*** Discrete choice (1 if battery replace, 0 if do not replace - at each time period)
*** Last period is skipped, since we only know state variable, but not decision

mata:
for (i=1; i<=P-1; i++) {
	for (j=1; j<=nobs; j++){	
		if (a[j,i+1]==1){
			d[j,i]=1
		}	
	} 
}
end


*** Discount factor
mata beta=0.95

*** Theta initial values (parameters to estimate)
mata theta=J(1,3,0)  
  
* Elements:
* [1,1] - replacement cost (theta_R)
* [1,2] - maintenance cost (theta_M1)
* [1,3] - maintenance cost (theta_M2)


********************************************************************************
*** PART X: GET THE CONTINUATION VALUE FOR EACH PERIOD AND STATE

* Define terminal condition
mata
T = J(5,1,.)
T[1,1] = 7
T[2,1] = 4
T[3,1] = 3
T[4,1] = 2
T[5,1] = 1
end


cap mata: mata drop get_continuation()
mata
function get_continuation(theta, beta, T){
	
	/* Input
	* thetas
	* Terminal continuation Value
	*/
	
	/* Output
	
	* C - Matrix with continuation values for each state and period
	
	*/
	
	C = J(5,5,.)

	// At period 5:
	C[1,5] = T[1,1]
	C[2,5] = T[2,1]
	C[3,5] = T[3,1]
	C[4,5] = T[4,1]
	C[5,5] = T[5,1]
	
	// For the other periods:
	for (t=4;t>=1;t--){
		v1 = -theta[1,1] + beta:*C[1,t+1] 
		v01 = -theta[1,2] - theta[1,3]+ beta:*C[2,t+1] 
		v02 = -theta[1,2]*2 - theta[1,3]*4 + beta:*C[3,t+1] 
		v03 = -theta[1,2]*3 - theta[1,3]*9 + beta:*C[4,t+1] 
		v04 = -theta[1,2]*4 - theta[1,3]*16 + beta:*C[5,t+1]
		
		C[1,t] = log(exp(v1) + log(exp(v01)))
		C[2,t] = log(exp(v1) + log(exp(v02)))
		C[3,t] = log(exp(v1) + log(exp(v03)))
		C[4,t] = log(exp(v1) + log(exp(v04)))
		}
	return(C)
	}
end

mata C = get_continuation(theta,beta,T)
mata C
********************************************************************************
*** PART XI: ESTIMATE THE FINITE HORIZON MODEL

mata a
*** Construct log-likelihood function	
mata: 
function loglike_finite(todo, theta, beta, d, a, T, y, g, H){
	
	/* Input 
	* todo - a scalar message variable from optimize() that indicates where to compute 1st and 2nd order derivatives;
	* theta - parameters which are going to be estimated
	* beta - discount factor
	* d - matrix of decisions
	* a - matrix of state variable values
	* T - Terminal Continuation Values 
	*/
	
	/* Output
	* y - the value of the objective function that is being optimized
	* g - the gradient vector
	* H - the Hessian matrix
	*/
	
	/* find a continuation value for each period */
	C = get_continuation(theta,beta,T)

	
	/* construct conditional value functions for each period*/

	v1 = J(1,5,.)
	v01 = J(1,5,.)
	v02 = J(1,5,.)
	v03 = J(1,5,.)
	v04 = J(1,5,.)	
	
	for (t=1;t<=4;t++){
		v1[1,t] = -theta[1,1] + beta*C[1,1] 
		v01[1,t] = -theta[1,2] - theta[1,3]+ beta*C[2,t+1]
		v02[1,t] = -theta[1,2]*2 - theta[1,3]*4 + beta*C[3,t+1]
		v03[1,t] = -theta[1,2]*3 - theta[1,3]*9 + beta*C[4,t+1]
		v04[1,t] = -theta[1,2]*4 - theta[1,3]*16 + beta*C[5,t+1]
	}
	
	
	/* paysoff diffrences: base category no replacement */ 

	payoffdiff1 = J(1,5,.)
	payoffdiff2 = J(1,5,.)
	payoffdiff3 = J(1,5,.)
	payoffdiff4 = J(1,5,.)
	
	for (t=1;t<=4;t++){
		
		payoffdiff1[1,t] = v1[1,t] - v01[1,t]
		payoffdiff2[1,t] = v1[1,t] - v02[1,t]
		payoffdiff3[1,t] = v1[1,t] - v03[1,t]
		payoffdiff4[1,t] = v1[1,t] - v04[1,t]
		
	}
	
	payoffdiff = (payoffdiff1 \ payoffdiff2 \ payoffdiff3 \ payoffdiff4)
	
	/* using conditional choice probabilities construct log-likelihood function */
	
	f = J(rows(d),cols(d)-1,0)

	for (i=1; i<=rows(d); i++) {
		for (t=1; t<=(cols(d)-1); t++){
			if (a[i,t] < 5) {
				f[i,t] = d[i,t]*log(1 - 1/(1+exp(payoffdiff[a[i,t],t]))) + (1-d[i,t])*log(1/(1+exp(payoffdiff[a[i,t],t])))
			}
		}
	}
	
	y = sum(f) 
	
	}

end



*** Call your optimization function
mata s=optimize_init()

*** Tell it where to find the evaluator function
mata optimize_init_evaluator(s,&loglike_finite())

*** Define type of evaluator - default is d0 (scalar)
mata optimize_init_evaluatortype(s,"d0")

/*	Remaining evaluator types are e.g.:
	- "d1" - same as "d0" and returns gradient rowvector
	- "d2" - same as "d1" and returns Hessian matrix
	- "gf0" - function() returns colvector value
	- "gf1" - same as "gf0" and returns score matrix
	- "gf2" - same as "gf1" and returns Hessian matrix   */

	
*** Define starting parameters 
mata optimize_init_params(s,J(1,3,0))

*** Indicate to optimizer which variables include as remaining function arguments
mata optimize_init_argument(s,1,beta)
mata optimize_init_argument(s,2,d)
mata optimize_init_argument(s,3,a)
mata optimize_init_argument(s,4,T)

*** Choose optimization method - default is Newton-Rhapson ("nr")
mata optimize_init_technique(s,"nr")

*** Returns the value of theta which maximizes loglike(theta).
mata thetahat = optimize(s)
mata thetahat
*** Standard Errors

mata hessian = optimize_result_Hessian(s)

mata hessian=-hessian

*** Take the inverse of a matrix

mata H = luinv(hessian)

*** Compute Standard Errors

mata SE = sqrt(diagonal(H))

*** Display estimated coefficents and their standard errors

mata thetahat \ SE'


log close


********************************************************************************
********************************************************************************
