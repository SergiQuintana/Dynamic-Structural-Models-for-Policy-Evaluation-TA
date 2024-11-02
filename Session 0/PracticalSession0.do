********************************************************************************
************ Summer School - Structural Micro (Practical Sessions) *************
********************************************************************************

********************************************************************************
******************* SESSION 0: Introduction to Mata

* The following script contains the brief introduction to Mata software. Due to
* time constraints we focus only on the features and functions useful in later
* sessions. 
* It includes: 

*     1. STATA <-> MATA
*     2. Matrix Operations
*	  3. Loops
*     4. Functions
*     5. Simple Example - Linear Regression

* Prepared by Sergi Quintana (BSE and UAB)
* Based on the work of Jacek Barszczewski (BSE and UAB)

********************************************************************************


**************************************
*** Administrative Commands

cap log close                       
set more off                        
clear all                 

*** set path to the folder with data
cd "C:\Users\Sergi\Dropbox\PhD\Teaching\Third Year\Structural Micro BSS\Structural Micro BSS\3. Materials"

log using Session0.log, replace text

********************************************************************************
*** PART I: Introduction

/*
	Two main ways of coding:
	1. Using .do files and preceding each command with 'mata'.
	2. Go to command window and type 'mata'

*/

*** Where to get help

help mata cholesky()


*** How to open mata

mata

end

*** Let's define matrix A 

mata A = (1, 2, 3 \ 4, 5, 6 \ 7, 8, 9)

*** Now do the same directly in mata command window

*** For the sequence of commands:

mata:
	A = (1, 2, 3 \ 4, 5, 6 \ 7, 8, 9)
	B = A*A
end

**************************************
***  Importing variables from Stata

*** Open the dataset in Stata

use "PracticalSession0_data.dta", clear
compress

describe

summarize

*** Export to Mata wage variable

mata wage = st_data(.,"ln_wage")

mata X = st_data(.,("age", "age2", "exper", "exper2", "educ", "female"))

*** Export only selected rows
* syntax matrixname = st_data(rows (if all, then .), varnames in Stata dataset)

mata X1 = st_data((1::10\300::359),("age", "age2", "exper", "exper2", "educ", "female"))

mata X2 = st_data((1::10\300::359),(3, 4, 5, 6, 7, 8))

*** To omit rows with missing value, add argument equal 0

mata X3 = st_data((1::10\300::359),(3, 4, 5, 6, 7, 8), 0)

*** Drop selected matrices

mata mata drop X3 X2

*** Delete all the variables

mata mata clear

*** Alternative command for transferring the data to Mata
* syntax: st_view(matrixname_., rows (if all, then .), varnames in Stata dataset)

mata st_view(wage=., ., "ln_wage")
mata wage

mata st_view(X=., ., ("age", "age2", "exper", "exper2", "educ", "female"))

mata st_view(X1=., (1::10\300::359), ("age", "age2", "exper", "exper2", "educ", "female"))
mata X1
*** Differences between st_data and st_view
* st_view does not copy the data from the dataset, so it uses less memory 
* st_view points to the stata columns 
* If you introduce changes to stata dataset, then values in the matrix also change
* stat_view is quicker and uses less memeory, but one needs to be careful
* when you use it

**************************************
***  Mata -> Stata

*** Create an interaction of female and education in Mata and save it in Stata

mata st_view(f=., ., "female")
mata st_view(e=., ., "educ")

mata educ_fem = f:*e

mata st_addvar("byte","educ_fem") // add variable to stata dataset (format, name)
mata st_store(.,"educ_fem",educ_fem) // asign values to existing stata var

* Stata numeric formats
* byte (max digit 2, only integers)
* int (max digit 4, only integers)
* long (max digit 9, only integers)
* float (max digit 7, real numbers)
* double (max digit 16, real numbers)
* by default stata creates all variables as flaot
* More prcise storage formats take up more memory -> the file is larger and Stata is slower

**************************************
***  Calling Mata inside Stata loop

*** Create new variable - experience relative to the mean

gen long exper3 = exper^3
rename exper exper1

forvalues i=1/3{
	
	mata st_view(exper=., ., "exper`i'")
	mata mean_exper = sum(exper)/rows(exper)
	mata exper`i'_avg = mean_exper:*J(rows(exper),1,1)
	mata st_addvar("float","exper`i'_avg") 
	mata st_store(.,"exper`i'_avg",exper`i'_avg) 
	gen float exper`i'_r = exper`i'/exper`i'_avg
}

********************************************************************************
*** PART II: Matrix Operations

**************************************
*** Defining matrices

*** Define row vector

mata a = (1, 2, 3, 4)

*** Define column vector

mata b = (1 \ 2 \ 3 \ 4)

*** Or simply transpose the row vector

mata b = a'

*** Column vector with series of numbers

mata c = (1::10)
mata c
*** Row vector with series of numbers

mata d = (1..10)
mata d

*** Back to matrices

mata A = (1, 2 \ 3, 4 )

mata B = J(2,2,0)

mata C = J(2,2,1)

mata D = diag(J(2,1,1))
mata D

*** or alternatively

mata D = I(2)

*** Construct E = | A B | 
*				  | C D |

mata E = (A, B \ C, D)

*** Or:

mata E = A, B \ C, D

*** Construct vector from matrix

mata e = vec(E)

*** Indexing

mata E[1,4] = 1

mata E[1,.] = J(1,4,1)


**************************************
*** Matrix Operations

*** Extract diagonal of given matrix

mata e = diagonal(E)

*** Extract lower- upper-triangle of given matrix

mata LowE = lowertriangle(E)

mata UpE = uppertriangle(E)

*** Sorting a matrix

mata sort(E,1) /*sort on a first column*/

mata sort(E,(1,2)) /*sort on a first and second column*/

*** Matrix sum
mata F = A+C

*** Sum of matrix and constant 

mata F = A:+2

*** Matrix multiplication

mata F = A*C

*** Element by element

mata F1 = A:*C
mata F2 = A:/F

*** Inverse of full rank, square matrix

mata A = (4,-3 \ -3,4)
mata A_inv = luinv(A)

*** Generalized inverse of positive-definite, symmetric matrix
 
mata A_inv = invsym(A)

*** Number of rows and columns

mata r = rows(E)

mata c = cols(E)

*** Computing the sum of rows and columns

mata r_sum = rowsum(E)

mata c_sum = colsum(E)

*** Sum of all elements of the matrix

mata e_sum = sum(E)

*** To see more usefull mata functions, one can use moremata package

* ssc install moremata
*help moremata 

mata mata describe

mata mata drop e_sum

mata mata clear



***************
* EXERCICE 1. 
webuse auto,clear

mata mata clear
mata st_view(price=., ., "price")
mata price2 = price:^2
mata price2[4] = 1000
mata st_addvar("float","price2")
mata st_store(.,"price2",price2)


****************
* EXERCICE 2
mata

// 1. 
a=2

A=J(4,4,2)
B=J(4,4,0) // first create matrix B only containing zeros

B[1,1]=2 // replace diagonal elements step-wise by 2
B[2,2]=2
B[3,3]=2
B[4,4]=2

// alternatively, checking out the mata standard functions (help m4_standard)
B=I(4)*2

// to display the elements, just type their name:
A
B
a

// 2. 
C=a*A
D=C+invsym(B)
E=(A,B)\(C,D)
b=E[1,1], E[1,3], E[1,5]

C
D
E
b

// 3. 
// Update matrix A:
A=A\((1::20),(2::21),(3::22),(4::23))
A

end 



********************************************************************************
*** PART III: Loops 

*** Fibonacci Sequence - define vector of missing values

mata N = 20

mata F = J(1,N,.)

*** Loop for

mata F[1,1] = 1

mata F[1,2] = 1

mata:
for(n=3;n<=N;n++){
	F[1,n] = F[1,n-1] + F[1,n-2]
}
end

*** Loop while

mata n=1

mata F = F, J(1,10,.)

mata:
while (n<=10){
	F[1,20+n] = F[1,20+n-1] + F[1,20+n-2]
	n++
}
end

*** Loop do while - sum of geometric series

mata:
geosum = 0 
r = 1/2 
n = 10 
i=1
end

mata:
do {
	geosum= geosum + r^(i-1)
	i++
	geosum
} while (i<=n)
end

*** Logical Conditions!

mata
n = 5
i = 1
while (1) {
printf("i=%g\n", i)
i++
if (i>n) {
break
}
}
printf("i=%g\n", i)
printf("done\n")
end

*** Exercices loops
mata:
A = J(10,20,1)
rows(A)
cols(A)
for(j=1;j<=rows(A);j++){
	for (p=1;p<=cols(A);p++){
		if ((j+p) <10){
			A[j,p] = (j+p)^2
		}
		else{
			A[j,p] = 0
		}
		
	}
}
A
end

mata 
A = J(10,15,2)
while (rows(A)<cols(A)){
	B = A[.,1::cols(A)-1]
	B[.,cols(A)-1] = B[.,cols(A)-1] + A[.,cols(A)]
	A = B
	A
}
end

********************************************************************************
*** PART IV: Functions

* Function without arguments or return
cap mata: mata drop hi()
mata
function hi(){
	printf("Hi!\n")
}
hi()
end

* Function with arguments

cap mata: mata drop func1()
mata
function func1(a){
	a+4
}
func1(3)
end

set matastrict off
* Function with return
cap mata: mata drop func2()
mata
function func2(a)
{
	b = a+2
	return(b)
}
b = func2(10)
b
end

* Function with declarations 
cap mata: mata drop zeros()
mata
real colvector zeros(real scalar c)
{
real colvector a
a = J(c,1,0)
return(a)
}
c = zeros(2)
c
end

*** EXERCICE!
*** Goal - define a function that returns a vector row vector of n independent 
*** random variables drawn from the Bernoulli distribution given success 
*** probability theta

*** Define vector size and theta

mata n = 10000

mata theta = 0.75

*** Erase function from the memory if it was specified before 
cap mata: mata drop bernoulli()

mata: 
function bernoulli(n,theta){
	A = uniform(1, n) 
	bern = J(1,n,0)
	for(i=1;i<=n;i++){
		if (A[1,i]<=theta){
			bern[1,i] = 1
		}
	}
	return(bern)
}
end

mata bern_vec = bernoulli(n,theta)

*** Compute fraction of ones
mata frac_one = sum(bern_vec)/n

*** Alternatively, one can create separate .do file with the function

*** Optinally, one can set seed to obtain the same draw of pseudo-random variables, i.e.:
set seed 31415926

********************************************************************************
*** PART V: Simple example - linear regression

/*
	Let's follow our example and try to estimate simple mincerian regression.
	As a regressors let's use all available variables that describe age, 
	education and its powers, as well as education and gender. 

*/

use "PracticalSession0_data.dta", clear

*** Check if all of the variables are non-missing

gen byte sample=(age!=. & age2!=. & exper!=. & exper2!=. & educ!=. & female!=.)

*** Let's 'export' data to Mata
*   Note: only rows where sample = 1 will be transferred to Mata

mata st_view(x=0,.,("age", "age2", "exper", "exper2", "educ", "female"), "sample")

mata st_view(y=0,.,("ln_wage"), "sample")

*** Include constant term

mata x=x,J(rows(x),1,1)

*** Compute OLS coefficients

mata b=invsym(x'*x)*x'*y

*** Compute standard errors
*   First: find residuals

mata e=y-x*b

* Variance-covariance matrix

mata v=(e'*e)/(rows(x)-cols(x))*invsym(x'*x)

* Compute standard errors

mata se=sqrt(diagonal(v))

* Obtain t-statistics

mata t=b:/se

* Compute p-value using ttail function, which takes as arguments the degrees 
* of freedom and the t-statistic and returns the probability

mata p=2*ttail(rows(x)-cols(x),abs(t))


*** Present your results in Mata

mata b,se,t,p


*** Compare your results with Stata output

reg ln_wage age age2 exper exper2 educ female


*** Exercice Function OLS

cap mata: mata drop myols()
mata
function ols(y,X){
	
	
}
end


********************************************************************************
*** PART VI: Optimization

help mf_optimize


cap mata: mata drop myfunc()
mata:mata clear
mata
   void myfunc(todo, x, y, g, H) {
      y = -4*x^3 + 3*x^2 + 25*x + 6
      }
maxval = optimize_init()
optimize_init_which(maxval, "max")
optimize_init_evaluator(maxval, &myfunc())
optimize_init_params(maxval, 1)
xmax = optimize(maxval)
xmax
end

** Plot

mata
  void myfunc(todo, x, y, g, H) {
   y = -4*x^3 + 3*x^2 + 25*x + 6
  }
    // maximum point
    maxval = optimize_init()
    optimize_init_which(maxval, "max")
    optimize_init_evaluator(maxval, &myfunc())
    optimize_init_params(maxval, 1)
       xmax = optimize(maxval)
    
    // minimum point
    minval = optimize_init()
    optimize_init_which(minval, "min")
    optimize_init_evaluator(minval, &myfunc())
    optimize_init_params(minval, -1)
       xmin = optimize(minval)
       xmin, xmax  // display the values
       // pass the values back to Stata
       st_local("minx", strofreal(xmin))
       st_local("maxx", strofreal(xmax))
end
twoway (function -4*x^3 + 3*x^2 + 25*x + 6, range(-4 4)), ///
 yline(0) xline(0) ///
 xline(`minx', lc(red)) xline(`maxx', lc(blue)) aspect(1) xsize(1) ysize(1) xlabel(-4(1)4)

log close

********************************************************************************
********************************************************************************
