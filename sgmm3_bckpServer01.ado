*==========================*
*  GMM PROGRAM             *
*==========================*
***************** DECKER BUILD AUGUST 1, 2012 ***************************************
* this will do system GMM for any number of equations (including 1 equation)
* 
* TO USE : must define global lists with all variables in equations 
*          dep.variables must be in lists y1 y2...
*          regressors must be in x1 x2 ... and instruments in z1 z2 ...
*          must define global macro G which has number of equations
* for example :
* global y1="ik"                             /* EQ 1: dep.var */
* global x1="const l.ik sk"                  /* EQ 1: rhs variables */
* global z1="const l.ik l.sk l2.ik l2.sk"    /* EQ 1: instruments (include exog vars here)  */  
* global G=1
*
* after exit this program will leave behind matrices: b2sls bgmm var std
* results will be posted to Stata so tests can be done using Stata's command
* for ex: test [eq1_]sk=[eq2_]h_sk   
*

capture program drop sgmm3_bckpServer01
program define sgmm3_bckpServer01
version 6.0
*set trace on

di "GMM started : $S_TIME  "
/* getting number of variables in each equation */


*di "numcountries..."
*di $nc







local l = 0       /* total number of instruments for the system */
local k = 0       /*                 regressors                 */
local g=1
while `g'<=$G {
  local l`g' : word count ${z`g'} 
  local k`g' : word count ${x`g'} 
  local l=`l'+`l`g''
  local k=`k'+`k`g''
  local rownms = "`rownms' ${x`g'}"
  *di "g = `g'"
  *di " ${x`g'}"
  local g=`g'+1 
}
*di "rownms = `rownms'"
*di "locals l1=`l1' l2=`l2' k1=`k1' k2=`k2' total l= `l' k=`k'" 

/* accumultaing matrices */

tempname zy temp
*di in b "accumulating matrices equation "
*di in b "accumulating matrices equation " _c
local g=1
while `g'<=$G {
*di "`g',"
*di "`g'," _c
tempname zz`g' zx`g' zy`g'
qui mat accum `temp'= ${z`g'} ${x`g'}, nocon       /* equation g */
mat `zz`g''=`temp'[1..`l`g'', 1..`l`g'']
mat `zx`g''=`temp'[1..`l`g'', `l`g''+1...]
qui mat accum `temp'= ${y`g'} ${z`g'}, nocon
mat `zy`g''=`temp'[2...,1]
local listzx="`listzx' `zx`g''"   /* accumulate lists of all zx and zz matrices */
local listzz="`listzz' `zz`g''"
mat `zy'= nullmat(`zy') \ `zy`g''    /* accumulating ZY - it is not diagonal */
local g=`g'+1 
}
*di "xg =  ${x1}"
*mata
*zz = st_matrix(listzx)
*end

*di in b "calculating b2sls..."

/* System 2SLS - should be equal to equation by equation */

tempname zx zz invz
diag `zx' = `listzx'
diag `zz' = `listzz'
mat `invz'=invsym(`zz')

*mat b2sls=invsym(`zx'' * `invz' * `zx' ) * `zx'' * `invz' * `zy'
* b2sls se obtendrá con el MATA
*di "getting varlist..."
local varlist = "$vrlst"
*di "`varlist'"
di in b "Mata beginnt..."
* The b is calculated with the pertinent Z :)
*di "-------------------------------------------------"
*di "-------------------------------------------------"
**************************************************MATA CODE call function
mata: getb2sls("`varlist'")
**************************************************
*di "-------------------------------------------------"
*di "-------------------------------------------------"
di in b "Mata done..."
*numObs proviene de Mata (se usó local). ;)
*di "numobs = $numobs"
*di "numobs = `numobs'"
*matrix list b2M
*di "hnames = ${hnames}"
	*matrix rownames b2Mata = ${x1}
	*matrix colnames b2Mata = ${hnames}
matrix rownames b2M = `rownms'
matrix b2sls = b2M

* y generamos el bgmm (sin los errores robustos) ...
mat bgmm = b2M
* construimos var igual sin "robustez", a la manera de R-D aunq es ineficiente...
local g = 1
	*	while `g'<=$G {
		*	local namet = "matvar"
		*	local listvs "`listvs' `namet'" 
		*	local g=`g'+1
	*	}
	*di "listvs = `listvs'"
	*	tempname tempv
	*	diag `tempv' = `listvs'
	*	mat var = `tempv'
mat var = matvar
*** if u want to see matrix var...
*di "matrix list var = "
*mat list var

local nused = `numobs'
global T=`nused'


mat `temp'=bgmm

local listu=""
local g = 1

while `g'<=$G {
  tempname b`g'           /* extract coefficients for equation g into b`g' */
  tempvar t`g' u`g'
  mat `b`g''=`temp'[1..`k`g'',.]'           /* note bg's are row vectors! */
  capture mat `temp'=`temp'[`k`g''+1...,.] 
  mat score `t`g''=`b`g''
  qui gen `u`g''=${y`g'}-`t`g''    
  capture drop u`g'
  qui gen u`g'=${y`g'}-`t`g''    /* NOTE - added permenent variable with residual */
	*di "mat list ug"
	*list u`g' if id == 1

*di " score for equation `g': and residual"
*sum `t`g'' `u`g'' ${y`g'}

local listu "`listu' `u`g''" /* this is the list of all residuals - use to calculate u'u */
local listup "`listup' u`g'"  /* list of permanent variables with residuals */
local nameu "`nameu' eq`g'"  /*  list of equation numbers - to name rows and columns of uu */
  local g=`g'+1
}
if $G>1 {  /* correlation matrix is only calculated if more then 1 equation is specified */
**** generating U'U matrix of variance-covariance of resuduals *****
*di "these are summary  of residuals:"
*sum `listu' $names

qui mat accum uu=`listup' , nocons    /* this is u'u for all equations */
mat uu=(1/`nused')*uu
if "$names"~="" {       local nameu "$names"        } 
  * for VAR only - $names is the list of all Y variable names 
  * for others - make the list later using Y variables, for now will use eq1 eq2 ...
mat colnames uu = `nameu'     
mat rownames uu = `nameu'
*mat list uu

di " "
di in g "Residuals correlation matrix "
pwcorr `listup' , sig
}

*di "b2Mata = "
*matrix list b2Mata
matrix list b2M
di ""
di "bgmm ist the same, but vectorized, litrly vec(b) *on mata"
di " "
di "GMM finished : $S_TIME"
di " "

end



















*******************************************
*** Auxilary programs needed to run GMM ***
*******************************************

************************************************************************************
************************************************************************************
************************************************************************************
************************************************************************************
************************************************************************************
************************************************************************************
************************************************************************************
**************************************************** Mata Code beginnt. ;)
version 6.0
mata:
function getb2sls(varname)
{
	/*  st_view(matTest0 =.,.,("l.h_gc l.h_gdp l.h_ca l.h_reer"))		*/
	/* Va a haber que modificar esto para que sea generalizable...		*/
	/* get parameters  */
pp = strtoreal(st_global("P"))
numvars = strtoreal(st_global("G"))
        /* numPais = strtoreal(st_global("numCountries"))                 */
numPais = strtoreal(st_global("nc"))
G = numvars
/* matrixData = st_data(.,("h_gc h_gdp h_ca h_reer")) */

matrixData = st_data(.,(varname))			
		/* past line for when function in stata */
				/* and next line only for notepad */ 
		/*   (when building Z 4 first time) */
	/*matrixData = st_data(.,("`varlist'"))   */

Y = matrixData
	/*  Y = st_data(.,("h_gc h_gdp h_ca h_reer"))				*/
	/* A_z = st_data(.,("l.h_gc l.h_gdp l.h_ca l.h_reer"))			*/
A_z = matrixData
	/*   matTest1[101..500,1..4]      para ver la matriz en un rango	*/
	/*   matTest1[135,1]      para ver un dato                      	*/
	/*  Limpiamos la matriz	... la peluqueamos para que sólo tenga valores	*/
	/*	en renglones COMPLETOS						*/
m = rows(matrixData)
n = cols(matrixData)
numtrim = m/numPais
	/* OJO: Assuming T is equal for each country.				*/
A = A_z
	/*  A será como una A_z pelona, y también pelamos a Y			*/
vec10000 = J(1,n,100000000)
vecpoints = J(1,n,.)
for(ii =  1; ii <= m; ii++){
	/* printf(" ii = %g \n", ii) */
	vecAz = A_z[ii,.]
	vecY = Y[ii,.]
	if( vecAz < vec10000 ){
		/* si todos son menores, son números, en otro caso s vacío...	*/
	}else {
		A[ii,.] = vecpoints
	}
	if( vecY < vec10000 ){
		/* si todos son menores, son números, en otro caso s vacío...	*/
	}else {
		Y[ii,.] = vecpoints
	}
}
	/*    A[8000..m,.]    ó    A[1..600,.]  Sólo para ver la matriz :)	*/
	/*    La idea ahora es crear la matriz Z a la manera de Arellano - Bond */
	/*    Vemos los límites de donde hay datos para cada país en A		*/

empiezanDatosPaisj = 0
terminanDatosPaisj = 0
/* indica si encontramos el primer dato */
flag = 0
/* indica si encontramos el último dato */
flag2 = 0
	/*    OJO: MUY importante, en IMV hay 44 países, pero puede cambiar	*/
	/* 	Lo mismo que en general la estructura de trimestres/años	*/
vec_EmpTermPaises = J(numPais,2,0)
vec10000 = J(1,cols(A_z),100000000)
for(j=1; j <= numPais; j++){
	flag = 0
	flag2 = 0
	counter = 1
	empiezanDatosPaisj = -1
	terminanDatosPaisj = -1
	while(flag2 == 0){
		if(A[numtrim*(j-1) + counter,.] < vec10000 && flag == 0){
			empiezanDatosPaisj = counter
			flag = 1
		}else if( !(A[numtrim*(j-1) + counter,.] < vec10000) && flag2 == 0 && flag == 1 ){
			terminanDatosPaisj =counter-1
			flag2 = 1
		}else if (counter == numtrim && flag2 ==0){
			flag2 = 1
			terminanDatosPaisj =counter
		}
		counter++
	}
	vec_EmpTermPaises[j,.] = (empiezanDatosPaisj , terminanDatosPaisj ) 
}
	/* hay que buscar mínimos comunes de años (o maximos).			*/
l= J(0,0,.)
	/* l guarda el lugar donde está el minmax(de hecho es vector y los ordena)*/
w = J(0,0,.)
	/* w no se usa */
maxindex(vec_EmpTermPaises[.,1],1,l ,w)
maxYearInit = vec_EmpTermPaises[l[1],1] 
l= J(0,0,.)
	/* l guarda el lugar donde está el minmax(de hecho es vector y los ordena)*/
w = J(0,0,.)
	/* w no se usa */
minindex(vec_EmpTermPaises[.,2],1,l ,w)
minYearEnd = vec_EmpTermPaises[l[1],2] 
	

	/*    Vemos los límites de donde hay datos para cada país en Y		*/
	/*    EncontramosminYearEnd  y de Y 					*/
empiezanDatosPaisj = 0
terminanDatosPaisj = 0
flag = 0
flag2 = 0
	/*    OJO: MUY importante, en IMV hay 44 países, pero puede cambiar	*/
	/* 	Lo mismo que en general la estructura de trimestres/años	*/
vec_EmpTermPaises = J(numPais,2,0)
vec10000 = J(1,cols(Y),100000000)
for(j=1; j <= numPais; j++){
	flag = 0
	flag2 = 0
	counter = 1
	empiezanDatosPaisj = -1
	terminanDatosPaisj = -1
	while(flag2 == 0){
		if(Y[numtrim*(j-1) + counter,.] < vec10000 && flag == 0){
			empiezanDatosPaisj = counter
			flag = 1
		}else if( !(Y[numtrim*(j-1) + counter,.] < vec10000)  && flag2 == 0 && flag == 1 ){
			terminanDatosPaisj = counter-1
			flag2 = 1
		}else if (counter == numtrim && flag2 ==0){
			flag2 = 1
			terminanDatosPaisj =counter
		}
		counter++
	}
	vec_EmpTermPaises[j,.] = (empiezanDatosPaisj , terminanDatosPaisj ) 
}
	/* hay que buscar mínimos comunes de años (o maximos).			*/
l= J(0,0,.)
	/* l guarda el lugar donde está el minmax(de hecho es vector y los ordena)*/
w = J(0,0,.)
	/* w no se usa */
maxindex(vec_EmpTermPaises[.,1],1,l ,w)
maxYearInit_Y = vec_EmpTermPaises[l[1],1] 
l= J(0,0,.)
	/* l guarda el lugar donde está el minmax(de hecho es vector y los ordena)*/
w = J(0,0,.)
	/* w no se usa */
minindex(vec_EmpTermPaises[.,2],1,l ,w)
minYearEnd_Y = vec_EmpTermPaises[l[1],2] 
T = minYearEnd_Y - maxYearInit_Y + 1
pT = T-pp



	/*   maxYearInit vec:	*/
mxYI = vec_EmpTermPaises[.,1]
	/*   minYearEnd vec:	*/
mnYE = vec_EmpTermPaises[.,2]
TT = mnYE - mxYI  + J(rows(vec_EmpTermPaises),1,1)
pTT = TT - pp*J(rows(TT),1,1)
timesMatrix = ( vec_EmpTermPaises, TT, pTT)


	/*Con minYearEnd y maxYearInit, creamos la matriz newY, a partir de Y.	*/
newY = J(0,G,.)
for(i=1;i<=numPais;i++){
	matTemp = Y[numtrim*(i-1) + mxYI[i] + pp ..numtrim*(i-1) + mnYE[i],.]
	newY = (newY \ matTemp)
}
newYfat = newY
	/*   newYfat contains the newY in the old fashion written */



/* First we build newYvec, that is Y in a vectorized fashion... */
sumpT = 0
newYvec = J(0,1,0)
for(i=1;i<=numPais;i++){
	betai = vec(newY[sumpT+1..sumpT + pTT[i],.])
	newYvec = (newYvec \ betai)
	sumpT = sumpT + pTT[i]
}

/*Creamos X, que es la matriz regresora de Y (newY)			*/
X = J(0,pp*G,.)
for(i=1;i<=numPais; i++){
	matr = J(0,pp*G,.)
	for(k=mxYI[i]+pp;k<=mnYE[i];k++){
		mk = J(1,0,.)
		for(p=1;p<= pp; p++){
			vec = Y[numtrim*(i-1)+k-p,.]
			mk = (mk , vec)
		}
		matr = (matr \ mk)
	}
	X = (X \ matr)
}



/* y después con X, creamos una nueva X, que es la matriz regresora de (newYvec) */
oldX = X
Xvec = J(0, G*G*pp, 0)
sumpT = 0
for(i=1;i<=numPais;i++){
	Xii = oldX[sumpT+1..sumpT + pTT[i],.]
	Xi = J(G*pTT[i],G*G*pp,0)
	for(g=1;g<=G;g++){
		Xi[pTT[i]*(g-1)+1..pTT[i]*g, G*pp*(g-1)+1..G*pp*g ] = Xii
	}
	
	Xvec = (Xvec \ Xi)
	sumpT = sumpT + pTT[i]
}

X = Xvec



/*
	/*    Creamos la matriz Z a la manera de Arellano - Bond  		*/
	/*      ... instrumentando correctamente...				*/
	/*    Ya tenemos maxYearInit, minYearEnd, y con p se puede hacer todo	*/
/*Z= J(m,n,.)*/
numforsZ = minYearEnd_Y - maxYearInit_Y - pp
Zold = J(numPais*(numforsZ+1),0,.)
for(i = 1; i <= numforsZ; i++){
	BigMatr = J(0,numvars*i,.)
	for(j=1;j<=numPais;j++){
		matr = J(numforsZ+1,numvars*i,0)
		/*matr[maxYearInit_Y + pp + 1..minYearEnd_Y ,.] = J(numforsZ,numvars*i,0)*/
		/*matr[1.. ,.] = J(numforsZ,numvars*i,0)*/
		/* construimos vecj */
		vecj = J(1,0,.)
		for(k=1;k<=i;k++){
			/*vecj = (vecj, A[(j-1)*200+maxYearInit+k-1,.])*/
			vecj = (vecj, Y[(j-1)*200+maxYearInit_Y+k-1,.])
		}
		matr[1 + i,.] = vecj
		BigMatr = (BigMatr \ matr)
	}
	/* Agregamos la matriz a Zold  = [Zold , BigMatr] 					*/
	Zold  = (Zold , BigMatr)
}

*/




l= J(0,0,.)
	/* l guarda el lugar donde está el minmax(de hecho es vector y los ordena)*/
w = J(0,0,.)
	/* w no se usa */
maxindex(pTT,1,l ,w)
maxpTT = pTT[l[1],1] 

	/* We build the Z with the correct instruments...	*/
	/* numforsZ = mnYE- mxYI - pp*J(rows(TT),1,1)  */
	/* 2 alternative ways of building the same matrix	*/
numforsZ = pTT - J(rows(TT),1,1)
gauss = (maxpTT *(maxpTT - 1) ) / 2
Z = J(0,G*(gauss + pp*(G-1)),0)
for(i=1;i<= numPais;i++){
	i
	Zi = J(G*pTT[i],G*( gauss + pp*(G-1) ), 0)
	for(g=1; g<= G; g++){
		sumlenpast = 0
		sumlen = 0
		yk = J(1,0,0)
		Zgg = J(pTT[i],gauss,0)
		for(k=2;k<= numforsZ[i]+1;k++){
			ygik = A[numtrim*(i-1)+mxYI[i]+k-2,g]
			if(k<=3){
				yk = (yk,ygik)
			}else{
				yk = (yk[2..cols(yk)],ygik)
			}
			lenk = cols(yk)
			sumlen = lenk + sumlenpast
			Zgg[k,sumlenpast+1..sumlen] = yk
			sumlenpast = sumlen
		}
		cInd = 1..G
		vec = cInd:!=g
		Yg = J(pTT[i],0,0)
		for(p=1;p<= pp; p++){
			Ygp =A[numtrim*(i-1)+mxYI[i]+pp-p..numtrim*(i-1)+mnYE[i]-p,.]
			Ygp2 = select(Ygp,vec)
			Yg = (Yg,Ygp2)
		}
		Zg = (Zgg,Yg)
		Zi[pTT[i]*(g-1)+1..pTT[i]*g,( gauss + (G-1)*(pp) )*(g-1) + 1..( gauss + (G-1)*(pp) )*g] = Zg
	}
	Z = (Z \ Zi)
}

printf("det(Z'*Z) = ")
det(Z'*Z)


/*  A[1..300,.] */
/* Z[1..20,.] */
/* A[378..384,.]  */

/* The last step in the beta building is gettting y, we */
numobs = rows(newY)
newY = newYvec
st_local("numobs", strofreal(numobs))
		/* numObs												*/
		/* newY													*/
printf("ZtX processing... ")
ZtX = Z'*X
printf("ZtY processing... ")
ZtY = Z'*newY
	/* ojo, ZZ ya está a la menos uno...						*/
printf("invsym(cross(Z,Z)) processing... ")
ZZ = invsym(cross(Z,Z))
printf("processing... ZtXtZZZtX = invsym(ZtX'* ZZ * ZtX ) processing... ")
ZtXtZZZtX = invsym(ZtX'* ZZ * ZtX )
printf("processing... b = ZtXtZZZtX * ZtX' * ZZ * ZtY processing... ")
b = ZtXtZZZtX * ZtX' * ZZ * ZtY
printf("passing b to b2Mata... ")
st_matrix("b2Mata", b)
	/*  b = vec(b)		*/
printf("passing b to b2M... ")
st_matrix("b2M", b)

printf("matvar processing = invsym(ZtX'* ZZ * ZtX )... ")
matvar=numobs*ZtXtZZZtX
/* mat var=`nused'*invsym(`zx'' * `W' * `zx') ... variance-covariance matrix */
		/*       st_global("numobs", numObs)					*/


		/*     agregamos perturbación a matvar para que sea pdf	*/
nn = rows(matvar)
epsilon = 0.001
perturb = J(nn,1,epsilon)
perturb = diag(perturb)
matvar = matvar + perturb 

printf("passing matvar to matvar... ")
st_matrix("matvar", matvar)
printf("passing b to b2slsMata... ")
st_matrix("b2slsMata", b)

	/* Xvec[1..10,.]		*/
	/* X[1..100,1..4] 		*/
	/* newY[1..100,.] 		*/
	/* Y[1..600,.] 			*/


}

end
**************************************************** Mata Code endet.. ;)
************************************************************************************
************************************************************************************
************************************************************************************
************************************************************************************
************************************************************************************
************************************************************************************
************************************************************************************








capture program drop diag
program define diag
** this will construct block diagonal matrix from the matrices supplied as list of parameters
** TO USE:  diag outname = zx1 zx2 zx3 
** outname will be the name of big block diagonal matrix 
local out `1'  /* get name of the output matrix into `out' local */
mac shift      /* skip next parameter - it is equal sign */ 
mac shift
tempname a temp1 temp2
mat `a'=`1'
mac shift
  local i 1      /* i is counter of equations - to use in naming columns and rows*/
*   matrix roweq `a' = eq1_
*   matrix coleq `a' = eq1_
while "`1'"~="" {
mat `temp1'=`a', J(rowsof(`a'),colsof(`1'),0)
mat `temp2'=J(rowsof(`1'),colsof(`a'),0), `1'
  local i=`i'+1                /* all idented lines  are just for*/
*   matrix roweq `1' = eq`i'_    /* corret labeling of variables and equations */
*   matrix coleq `1' = eq`i'_
  local cola : colnames(`a')
  local rowa : rownames(`a')
  local col1 : colnames(`1')
  local row1 : rownames(`1')
   local reqa : roweq(`a')
   local req1 : roweq(`1') 
   local ceqa : coleq(`a')
   local ceq1 : coleq(`1')   
mat `a'= `temp1' \ `temp2'
  mat colnames `a' = `cola' `col1' 
  mat rownames `a' = `rowa' `row1'
   mat roweq `a'=`reqa' `req1'
   mat coleq `a'=`ceqa' `ceq1'
mac shift
}
mat `out'=`a'
end


capture program drop matroot
program define matroot
 * this program calculates matrix with square roots of each element in the
 * incoming matrix
 * to use :  matroot input output
local input="`1'"    /* unload input vector into name input */
local output="`2'"
local rows=rowsof(`1')
local cols=colsof(`1')
mat `output'=`input'   /* create "dummy" matrix which then will be replaced */
local i=1
while `i'<=`rows'{
  local j=1
  while `j'<=`cols'{
     if sqrt(`input'[`i',`j'])==. { di in white "negative values in input matrix - cannot take sqrt" }
     mat `output'[`i',`j']=sqrt(`input'[`i',`j'])
     local j=`j'+1
    }
  local i=`i'+1
}
end


*capture program drop element
program define element
 * this program will perform element by element operation on matricies
 * possible functions will be +, -, *, /
 * the format is
 * element  output = input1 * input2     - for multiplication
 * element  output = input1 / input2     - for division
 * where * is in place of the operation that is needed to perform
 /* parameters 1 - output matrix
                2 - equal sign, for better readability
                3 - input matrix 1
                4 - operation to perform 
                5 - input matrix 2  */
  * both input patricies must be same size - for now

local input1="`3'"    /* unload input vector into name input */
local input2="`5'"
local o="`4'"           /* operation */
local output="`1'"
local rows=rowsof(`input1')
local cols=colsof(`input1')
mat `output'=`input1'   /* create "dummy" matrix which then will be replaced */
local i=1
while `i'<=`rows'{
  local j=1
  while `j'<=`cols'{
        if "`o'"=="/" & `input2'[`i',`j']==0 { di in white "division by zero in `i' row `j' col" }
     mat `output'[`i',`j']=`input1'[`i',`j']`o'`input2'[`i',`j']
     local j=`j'+1
    }
  local i=`i'+1
}
end


capture program drop correl
program define correl
 * this program will generate correlation matrix from the covariance matrix


local input1="`3'"    /* unload input vector into name input */
local input2="`5'"
local o="`4'"           /* operation */
local output="`1'"
local rows=rowsof(`input1')
local cols=colsof(`input1')
mat `output'=`input1'   /* create "dummy" matrix which then will be replaced */
local i=1
while `i'<=`rows'{
  local j=1
  while `j'<=`cols'{
        if "`o'"=="/" & `input2'[`i',`j']==0 { di in white "division by zero in `i' row `j' col" }
     mat `output'[`i',`j']=`input1'[`i',`j']`o'`input2'[`i',`j']
     local j=`j'+1
    }
  local i=`i'+1
}
end


