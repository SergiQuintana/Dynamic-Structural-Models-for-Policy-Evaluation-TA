********************************************************************************
************ Summer School - Structural Micro (Practical Sessions) *************
********************************************************************************

********************************************************************************
********* SESSION 4: - CCP Estimation: Finite Dependence

* The following script contains the solution of labor force participation
* problem for CCP framework, discussed during practical sessions. 
* It includes:  

*     1. Estimation of probabilities
*     2. Construction of likelihood function
*     3. Estimation of structural model

* Prepared by Jacek Barszczewski (BSE and UAB)
* Contact: jacek.barszczewski[at]bse.eu

**************************************
*** Administrative Commands

cap log close                       
set more off                        
clear all                 

*** set path to the folder with data
cd "C:\Users\Sergi\Dropbox\PhD\Teaching\Third Year\Structural Micro BSS\Structural Micro BSS\3. Materials"

log using Session4.log, replace text

********************************************************************************
*** PART I: LOAD THE DATA


*** Import decision matrix from excel file

import excel "PracticalSession4_data.xlsx", sheet("decision") firstrow clear

*** NOTE: each row corresponds to labor force participation decision of given 
*** individual for years 2008-2013

*** Define decision matrix in mata

mata d = st_data(.,("B","C","D","E","F"))

*** Import state matrix from excel file

import excel "PracticalSession4_data.xlsx", sheet("state") firstrow clear

*** Define state matrix

mata s = st_data(.,("B","C","D","E","F"))

*** Import gender information matrix from excel file

import excel "PracticalSession4_data.xlsx", sheet("gender") firstrow clear

*** Define state matrix

mata w = st_data(.,("sex"))

*** Define size variables

mata p = cols(d)

mata nobs = rows(d)

********************************************************************************
*** PART II: PREDICTED PROBABILITIES - FREQUENCIES

*** Set number of states
mata S = max(s) + 1

*** Create vector of states
mata state = (0::S-1)

*** Calculate conditional choice probability of working
mata pW = J(S,1,0) /* probability of working */ 

mata:
for(j=0;j<=length(state)-1;j++){          
	
	pW[j+1,1]=sum(d:*(s:==j))/sum((s:==j))

}	
end

*** Calculate conditional choice probability of staying home*/
mata pH = 1 :- pW
mata pH

********************************************************************************
*** PART III: SOLVE THE MODEL

*** VERSION 1 - NO HETEROGENEITY BY GENDER 

*** Discount factor
mata beta = 0.95

mata:
function loglikeH_FD(todo,phi,beta,d,s,pW,y,g,H){
 
	/* Input 
	* todo - a scalar message variable from optimize() that indicates where to compute 1st and 2nd order derivatives;
	* phi - parameters which are going to be estimated
	* beta - discount factor
	* d - matrix of decisions
	* s - matrix of state variables (working experience)
	* pW - estimated conditional choice probability of working
	*/
	
	/* Output
	* y - the value of the objective function that is being optimized
	* g - the gradient vector
	* H - the Hessian matrix
	*/
	
	P = cols(d) /* number of periods */
	
	/* select probability for each observation and period based on current state */
	
	pW_obs = pW[s[.,1]:+1,1] /* assign probability of working for all observation in the first period */
	
	/* use the loop to assign probabilities of working for remaining periods */
	
	for(i=2;i<=P;i++){
		pW_obs = pW_obs, pW[s[.,i]:+1,1]
	}
	
	/* calculate probability of staying home for all observations and all periods */
	pH = (1:-pW \ 1) 
	
	/* select probability for each observation and period based on current state */
	
	pH_obs = pH[s[.,1]:+1,1] /* assign probability of working for all observation in the first period */
	
	/* use the loop to assign probabilities of staying home for remaining periods */
	
	for(i=2;i<=P;i++){
		pH_obs = pH_obs, pH[s[.,i]:+1,1]
	}
	
	/* calculate log probabilities for the payoffs difference */
	lnpW_obs = log(pW_obs) /* log probability of working conditional on working experience*/
	lnpH_obs = log(pH_obs) /* log probability of staying home conditional on working experience*/
	
	/* utility of working */
	/* quadratic */
	/*u_w = phi[1,1] :+ phi[1,2]:*s + phi[1,3]:*s:^2 */
	/*linear */
	/*u_w = phi[1,1] :+ phi[1,2]:*s */
	/* cubic */
	u_w = phi[1,1] :+ phi[1,2]:*s + phi[1,3]:*s:^2 + phi[1,4]:*s:^3
	
	/* calculate payoffs difference */
	payoffsdiff = (1-beta):*u_w + beta:*(lnpW_obs - lnpH_obs)
	
	/* calculate likelihood function */
	f = ((1:-d):*log(1:/(1:+exp(payoffsdiff))) + d:*log(1 :- 1:/(1:+exp(payoffsdiff))))
	
	/* select only observations for exp < 40 */
	f = f:*(s:<40)
	
	/* sum likelihood across observations and periods */
	y = sum(f)
	
	}
end

*** Call your optimization function
mata s1=optimize_init()

*** Tell it where to find the evaluator function
mata optimize_init_evaluator(s1,&loglikeH_FD())

*** Define type of evaluator - default is d0 (scalar)
mata optimize_init_evaluatortype(s1,"d0")
	
*** Define starting parameters 
mata optimize_init_params(s1,J(1,4,0))

*** Indicate to optimizer which variables include as remaining function arguments
mata optimize_init_argument(s1,1,beta)
mata optimize_init_argument(s1,2,d)
mata optimize_init_argument(s1,3,s)
mata optimize_init_argument(s1,4,pW)

*** Choose optimization method - default is Newton-Rhapson ("nr")
mata optimize_init_technique(s1,"nr")

*** Returns the value of theta which maximizes loglike(theta).
mata phihat = optimize(s1)

*** Standard Errors

mata hessian = optimize_result_Hessian(s1)

mata hessian=-hessian

*** Take the inverse of a matrix

mata H = luinv(hessian)

*** Compute Standard Errors
mata SE = sqrt(diagonal(H))

mata (phihat \ SE')


*** VERSION 2 - HETEROGENEITY BY GENDER 

*** Discount factor
mata beta = 0.95

mata:
function loglikeH_FD_g(todo,phi,beta,d,s,w,pW,y,g,H){
 
	/* Input 
	* todo - a scalar message variable from optimize() that indicates where to compute 1st and 2nd order derivatives;
	* phi - parameters which are going to be estimated
	* beta - discount factor
	* d - matrix of decisions
	* s - matrix of state variables (working experience)
	* w - matrix of woman dummies 
	* pW - estimated conditional choice probability of working
	*/
	
	/* Output
	* y - the value of the objective function that is being optimized
	* g - the gradient vector
	* H - the Hessian matrix
	*/
	
	P = cols(d) /* number of periods */
	
	/* select probability for each observation and period based on current state */
	
	pW_obs = pW[s[.,1]:+1,1] /* assign probability of working for all observation in the first period */
	
	/* use the loop to assign probabilities of working for remaining periods */
	
	for(i=2;i<=P;i++){
		pW_obs = pW_obs, pW[s[.,i]:+1,1]
	}
	
	/* calculate probability of staying home for all observations and all periods */
	pH = (1:-pW \ 1) 
	
	/* select probability for each observation and period based on current state */
	
	pH_obs = pH[s[.,1]:+1,1] /* assign probability of working for all observation in the first period */
	
	/* use the loop to assign probabilities of staying home for remaining periods */
	
	for(i=2;i<=P;i++){
		pH_obs = pH_obs, pH[s[.,i]:+1,1]
	}
	
	/* calculate log probabilities for the payoffs difference */
	lnpW_obs = log(pW_obs) /* log probability of working conditional on working experience*/
	lnpH_obs = log(pH_obs) /* log probability of staying home conditional on working experience*/
	
	/* utility of working */
	/* quadratic */
	/* men */
	u_w_m = phi[1,1] :+ phi[1,2]:*s + phi[1,3]:*s:^2
	/* women */
	u_w_w = phi[1,4] :+ phi[1,5]:*s + phi[1,6]:*s:^2
	
	/* calculate payoffs difference */
	payoffsdiff = (1-beta):*(u_w_m:*(1:-w) + u_w_w:*w) + beta:*(lnpW_obs - lnpH_obs)
	
	/* calculate likelihood function */
	f = ((1:-d):*log(1:/(1:+exp(payoffsdiff))) + d:*log(1 :- 1:/(1:+exp(payoffsdiff))))
	
	/* select only observations for exp < 40 */
	f = f:*(s:<40)
	
	/* sum likelihood across observations and periods */
	y = sum(f)
	
	}
end

*** Call your optimization function
mata s1=optimize_init()

*** Tell it where to find the evaluator function
mata optimize_init_evaluator(s1,&loglikeH_FD_g())

*** Define type of evaluator - default is d0 (scalar)
mata optimize_init_evaluatortype(s1,"d0")
	
*** Define starting parameters 
mata optimize_init_params(s1,J(1,6,0))

*** Indicate to optimizer which variables include as remaining function arguments
mata optimize_init_argument(s1,1,beta)
mata optimize_init_argument(s1,2,d)
mata optimize_init_argument(s1,3,s)
mata optimize_init_argument(s1,4,w)
mata optimize_init_argument(s1,5,pW)

*** Choose optimization method - default is Newton-Rhapson ("nr")
mata optimize_init_technique(s1,"nr")

*** Returns the value of theta which maximizes loglike(theta).
mata phihat = optimize(s1)

*** Standard Errors

mata hessian = optimize_result_Hessian(s1)

mata hessian=-hessian

*** Take the inverse of a matrix

mata H = luinv(hessian)

*** Compute Standard Errors
mata SE = sqrt(diagonal(H))

mata (phihat \ SE')

log close


********************************************************************************
********************************************************************************
