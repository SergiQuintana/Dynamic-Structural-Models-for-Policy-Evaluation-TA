********************************************************************************
************ Summer School - Structural Micro (Practical Sessions) *************
********************************************************************************

********************************************************************************
********* SESSION 3 - Conditional Choice Probability (CCP) Estimation **********

* The following script contains the solution of bike maintanance problem 
* for CCP framework, discussed during practical sessions. 
* It includes: 

*     1. Standard NFXP algorithm, 
*     2. Estimation of probabilities
*     3. Construction of likelihood function
*     4. Estimation of A=5 structural model

* Prepared by Sergi Quitnana(BSE and UAB)



**************************************
*** Administrative Commands

cap log close                       
set more off                        
clear all                                  

*** set path to the folder with data
cd "C:\Users\Sergi\Dropbox\PhD\Teaching\Third Year\Structural Micro BSS\Structural Micro BSS\3. Materials"

log using Session3.log, replace text

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



********************************************************************************
*** PART IV: PREDICTED PROBABILITIES

*** Relative frequencies

* Declare vectors of frequencies p1 and p0
mata p1 = J(1,A,0)
mata p0 = J(1,A,0)

* Probability of replacement conditional on value of state variables  
* using frequency method
mata p1[1,1]=sum(d:*(a[.,1::5]:==1))/sum((a[.,1::5]:==1))
mata p0[1,1]=1-p1[1,1]

mata p1[1,2]=sum(d:*(a[.,1::5]:==2))/sum((a[.,1::5]:==2))
mata p0[1,2]=1-p1[1,2]

mata p1[1,3]=sum(d:*(a[.,1::5]:==3))/sum((a[.,1::5]:==3))
mata p0[1,3]=1-p1[1,3]

mata p1[1,4]=sum(d:*(a[.,1::5]:==4))/sum((a[.,1::5]:==4))
mata p0[1,4]=1-p1[1,4]

mata p1[1,5]=sum(d:*(a[.,1::5]:==5))/sum((a[.,1::5]:==5))
mata p0[1,5]=1-p1[1,5]

mata p1


********************************************************************************
*** PART V: SOLVE HOTZ & MILLER (1993)

*** Construct log-likelihood function
mata: 
function loglikeH(todo, theta, beta, d, a, p1, y, g, H){
	
	/* Input 
	* todo - a scalar message variable from optimize() that indicates where to compute 1st and 2nd order derivatives;
	* theta - parameters which are going to be estimated
	* beta - discount factor
	* d - matrix of decisions
	* a - matrix of state variable values
	* p1 - estimated conditional choice probabilities
	*/
	
	/* Output
	* y - the value of the objective function that is being optimized
	* g - the gradient vector
	* H - the Hessian matrix
	*/
	
	/* difference between conditional value functions: base category replacement */	

	payoffdiff1_H=theta[1,1]-theta[1,2]-theta[1,3]+beta*(log(p1[1,1])-log(p1[1,2])) /* v_01 - v1 */    
	payoffdiff2_H=theta[1,1]-theta[1,2]*2-theta[1,3]*4+beta*(log(p1[1,1])-log(p1[1,3])) /* v_02 - v1 */    
	payoffdiff3_H=theta[1,1]-theta[1,2]*3-theta[1,3]*9+beta*(log(p1[1,1])-log(p1[1,4])) /* v_03 - v1 */    
	payoffdiff4_H=theta[1,1]-theta[1,2]*4-theta[1,3]*16+beta*(log(p1[1,1]-log(p1[1,5]))) /* v_04 - v1 */   
	payoffdiff5_H=theta[1,1]-theta[1,2]*5-theta[1,3]*25+beta*(log(p1[1,1]-log(p1[1,5]))) /* v_05 - v1 */ 
	payoffdiff_H = (payoffdiff1_H \ payoffdiff2_H \ payoffdiff3_H \ payoffdiff4_H\payoffdiff5_H)
	
	/* using conditional choice probabilities construct log-likelihood function */
	h = J(rows(d),cols(d),0)

	for (i=1; i<=rows(d); i++) {
		for (j=1; j<=cols(d); j++){
				h[i,j] = d[i,j]*log(1/(1+exp(payoffdiff_H[a[i,j],1]))) + (1-d[i,j])*log(1-1/(1+exp(payoffdiff_H[a[i,j],1])))
		}
	}
	
	y = sum(h)
	
}
end

*** Call your optimization function
mata s1=optimize_init()

*** Tell it where to find the evaluator function
mata optimize_init_evaluator(s1,&loglikeH())

*** Define type of evaluator - default is d0 (scalar)
mata optimize_init_evaluatortype(s1,"d0")
	
*** Define starting parameters 
mata optimize_init_params(s1,J(1,3,0))

*** Indicate to optimizer which variables include as remaining function arguments
mata optimize_init_argument(s1,1,beta)
mata optimize_init_argument(s1,2,d)
mata optimize_init_argument(s1,3,a)
mata optimize_init_argument(s1,4,p1)

*** Choose optimization method - default is Newton-Rhapson ("nr")
mata optimize_init_technique(s1,"nr")

*** Using relative frequencies
mata thetahat_freq = optimize(s1)
mata thetahat_freq

*** Standard Errors

mata hessian_freq = optimize_result_Hessian(s1)

mata hessian_freq=-hessian_freq

*** Take the inverse of a matrix

mata H_freq = luinv(hessian_freq)

*** Compute Standard Errors
mata SE_freq = sqrt(diagonal(H_freq))


********************************************************************************
*** PART VI: AGUIRREGABIRIA AND MIRA: Update the CCPs

mata
function update_ccps(theta,oldccp,beta){
	
	/* Input 
	* theta --> current estimate of theta
	* oldccp --> previous ccps used to estimate theta
	* beta --> time discout factor
	*/
	
	/* Output
	* newccp --> ccps updated with our model
	*/
	
	payoffdiff1_H=theta[1,1]-theta[1,2]-theta[1,3]+beta*(log(oldccp[1,1])-log(oldccp[1,2])) /* v_01 - v1 */    
	payoffdiff2_H=theta[1,1]-theta[1,2]*2-theta[1,3]*4+beta*(log(oldccp[1,1])-log(oldccp[1,3])) /* v_02 - v1 */    
	payoffdiff3_H=theta[1,1]-theta[1,2]*3-theta[1,3]*9+beta*(log(oldccp[1,1])-log(oldccp[1,4])) /* v_03 - v1 */    
	payoffdiff4_H=theta[1,1]-theta[1,2]*4-theta[1,3]*16+beta*(log(oldccp[1,1]-log(oldccp[1,5]))) /* v_04 - v1 */   
	payoffdiff5_H=theta[1,1]-theta[1,2]*5-theta[1,3]*25+beta*(log(oldccp[1,1]-log(oldccp[1,5]))) /* v_05 - v1 */ 
	payoffdiff_H = (payoffdiff1_H \ payoffdiff2_H \ payoffdiff3_H \ payoffdiff4_H\payoffdiff5_H)
	
	newccps = 1:/(1:+exp(payoffdiff_H))
	
	return(newccps)
	}

end
mata newccps = update_ccps(thetahat_freq,p1,beta)
mata newccps

********************************************************************************
*** PART VI: AGUIRREGABIRIA AND MIRA: ALGORITHM


mata 

function optimizer(d,a,p1,theta0,beta){
	
			s1=optimize_init()
			// Tell it where to find the evaluator function
			optimize_init_evaluator(s1,&loglikeH())

			// Define type of evaluator - default is d0 (scalar)
			optimize_init_evaluatortype(s1,"d0")
				
			// Define starting parameters 
			optimize_init_params(s1,theta0)

			// Indicate to optimizer which variables include as remaining function arguments
			optimize_init_argument(s1,1,beta)
			optimize_init_argument(s1,2,d)
			optimize_init_argument(s1,3,a)
			optimize_init_argument(s1,4,p1)

			// Choose optimization method - default is Newton-Rhapson ("nr")
			optimize_init_technique(s1,"nr")

			// Using relative frequencies
			thetahat = optimize(s1)
			return(thetahat)
	
}

end



mata
function aguirremira(d,a,p_initial,theta,beta,iterations){

		/* Input 
		* p_initial --> initial estimate of ccps (data frequencies)
		* theta --> initial guess of theta
		* beta --> time discout factor
		* iterations --> number of iterations in the aguirregabira mira algorithm
		*/
	
		/* Output
		* theta --> estimated parameters improved
		*/
		
		// repeat the proces the desired amount of periods
		for (it =1; it<= iterations;it++){
			
			// Get ccps
			if (it==1){
				newccps = p_initial
			}
			if (it > 1){
				newccps = update_ccps(theta,p_initial,beta)
				newccps = newccps'
			}
			// Update the estimates of theta:			
			thetanew = optimizer(d,a,newccps,theta,beta)
			// Update theta
			theta = thetanew
			theta
		}
		return(thetanew)
}
end

mata theta0 = J(1,3,0)
mata thetahat = aguirremira(d,a,p1,theta0,beta,2)


log close


********************************************************************************
********************************************************************************
