
	**************************************************************************
	
	// Class Code
	
	

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

log close


********************************************************************************
********************************************************************************
