#ifndef MYUTILS_H_
#define MYUTILS_H_


#include <R.h>
/* FlushConsole*/
#include <Rinternals.h>
/* rng's, log, etc..*/
#include <Rmath.h>
/*  R_alloc*/
#include <R_ext/RS.h>
/* BLAS */
#include <R_ext/BLAS.h>
/* Lapack*/
#include <R_ext/Lapack.h>
/* Linpack : */
#include <R_ext/Linpack.h>
/* CheckUserInterrupt*/
#include <R_ext/Utils.h>
/*Fortran routines DQRLS etc...*/
#include <R_ext/Applic.h>

//#define DEBUG_MSG

/* ####################################### */
/*Define struct for blockwise QR-updates */
typedef struct {
	int qA; //length of active param-vector
	int qI; //length of inactive param-vector
	int indA1; // 1st element of active param
	int indA2; // last element of active param
	int reject;

	double *ya; // blockwise residuals (=y if only 1 block), rescaled with 1/sqrt(sigma2), with qA zeros attached
	double *Xa; // active cols of X, with qA x qA sqrt(Precision) attached
	double *Xi; // inactive X(colmns of X not in this block)
	double *coefI; //inactive coefficients

	double *qraux; //helper containers (qa x 1) for dqrdfc, dqrcf
	double *work;
	int *pivots;

	double *m; //containers for mean and drift of proposal
	double *err;
} XBlockQR;


void R_FreeXBlockQR(XBlockQR *blocks, int n){
	for (int i =0; i < n; i++) {
		R_Free(blocks[i].ya);
		R_Free(blocks[i].Xa);
		R_Free(blocks[i].Xi);
		R_Free(blocks[i].coefI);

		R_Free(blocks[i].qraux);
		R_Free(blocks[i].work);
		R_Free(blocks[i].pivots);

		R_Free(blocks[i].m);
		R_Free(blocks[i].err);

	}
	R_Free(blocks);
}


/*Aux. function:  add diag to diagonal of qxq-Matrix res */
static R_INLINE void addDiag (double *res, double *diag, int q)
{
	for(int j=0; j < q; j++){
		res[(j*q)+j] += diag[j];
	}
	return;
}

/*Aux. function:  divide vector by a scalar */
static R_INLINE void divVecByScalar (double *res, double *vec, double scalar, int q)
{
	for(int j=0; j < q; j++)
		res[j] = vec[j]/scalar;
	return;
}

/*Aux. function:  divide scalar by vector*/
static R_INLINE void divScalarByVec (double *res, double *vec, double scalar, int q)
{
	for(int j=0; j < q; j++)
		res[j] = scalar/vec[j];
	return;
}

/*Aux. function:  square elements of a vector*/
static R_INLINE void squareVec(double *res, double *vec, int q)
{
	for(int j=0; j < q; j++)
		res[j] = vec[j]* vec[j];
	return;
}

/*Aux. function: replicate elements in vector each times
 * res has length sum(each)
 * lenVec is length of vec and each
 * */
static R_INLINE void repEach (double *res, double *vec, int *each, int lenVec)
{
	int i, j, begin=0;
	for(i=0; i < lenVec; i++){
		for(j=0; j < each[i]; j++){
			res[begin+j] = vec[i];
		}
		begin += each[i];
	}
	return;
}

///*Aux. function:  set an array to zero */
static R_INLINE void setToZero (double *res, int len)
{
	for(int j=0; j < len; j++)
		res[j] = 0.0;
	return;
}

/*Aux. function:  sum vector by (sorted!) index */
/* index runs from 1 to length(vec)
 * res : the vector to be returned, length = no. of levels in index
 * vec: the vector to be summed up
 * index: same length as vec, containing numbers 1 to length(res), mapping entries in vec to entries in res
 * lengthVec: length of vector and index
 * */
static R_INLINE void sumVecByIndex (double *res, double *vec, int *index, int lengthVec, int lengthRes)
{
	int j, indexNow;
	for(j=0; j < lengthRes; j++){
		res[j] = 0;
	}
	for(j=0; j < lengthVec; j++){
		indexNow = index[j]-1;
		res[indexNow] += vec[j];
	}
	return;
}

/*Aux. function:  get standard deviation of a vector */
static R_INLINE double sdVec (double *vec, int lengthVec)
{
	int i;
	double mean=0.0, sqDev = 0.0, res;
	for(i=0; i < lengthVec; i++) mean += vec[i];
	mean = mean / lengthVec;
	for(i=0; i < lengthVec; i++) sqDev += R_pow_di((vec[i]-mean),2);
	res = sqrt(sqDev / (lengthVec-1));
	return(res);
}

// remove submatrix with upper left corner (index(r1, c1)), lower right corner (index(r2,c2)) from oldMat with dimensions rOld, cOld
// and return the remaining matrices (r1,0),(r2,c1-1); (r1,c2+1),(r2,cOld)
// uses zero based coordinates for corners but number of rows and columns for rOld, cOld!!
static void removeSubmatrix(double *newMat, double *oldMat, int rOld, int cOld, int r1, int c1, int r2, int c2)
{
	int j, i, k=0;
	if(c1 > 0){
		//left submatrix
		for(j=0; j < c1; j++){
			for(i=r1; i < (r2+1); i++){
				newMat[k] = oldMat[j*rOld + i];
				k++;
			}
		}
	}
	if(c2+1 < cOld){
		//right submatrix
		for(j = (c2 + 1); j < cOld; j++){
			for(i = r1; i < (r2+1); i++){
				newMat[k] = oldMat[j*rOld + i];
				k++;
			}
		}
	}
}

static double logLik(double* y, double* eta, int family, double* scale, int n){
	double ret = 0.0;
	double sd = 1/scale[0];
	if(family == 0){ //gaussian
		for(int i=0; i < n; i++){
			ret += dnorm(y[i], eta[i], sd, 1);
		}
		return ret;
	}
	if(family == 1){ //binomial
		for(int i=0; i < n; i++){
			ret += dbinom(y[i]*scale[i], scale[i], 1/(1+exp(-eta[i])), 1);
		}
		return ret;
	}
	if(family == 2){ //poisson
		for(int i=0; i < n; i++){
			ret += dpois(y[i], exp(eta[i]), 1);
		}
		return ret;
	}
	return 0.0;
}

//check for NAs in vec
static int anyNA(double *vec, int length){
	int ret = 0;
	for(int i = 0; i<length; i++){
		if(ISNAN(vec[i])){
			ret = 1;
			return(ret);
		}
	}
	return ret;
}

void debugMsg(char *format, ...){
#ifdef DEBUG_MSG
	va_list(ap);

	va_start(ap, format);
	Rvprintf(format, ap);
	va_end(ap);
#endif
}

static void printDouble(char *name, double *vec, int len)
{
	debugMsg("%s: (", name);
	if(len >0){//catch attempts to print empty vectors
		for(int j=0; j<(len-1); j++) debugMsg("%f, ", vec[j]);
		debugMsg("%f", vec[(len-1)]);
	}
	debugMsg(")\n");
	return;
}

/* old functions / debugging stuff starts here
static void printInt(char *name, int *vec, int len)
{
	debugMsg("%s: (", name);
	if(len >0){//catch attempts to print empty vectors
		for(int j=0; j<(len-1); j++) debugMsg("%d, ", vec[j]);
		debugMsg("%d", vec[(len-1)]);
	}
	debugMsg(")\n");
	return;
}
static void printArray(char *name, double *arr, int rows, int cols)
{
	debugMsg("%s: \n	", name);
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			debugMsg("%f, ", arr[j*rows + i]);
		}
		debugMsg("\n	");
	}
	return;
}



//extract submatrix with upper left corner (index(r1, c1)), lower right corner (index(r2,c2)) from oldMat with dimensions rOld, cOld
//uses zero based coordinates for corners but number of rows and columns for rOld, cOld!!
static void getSubmatrix(double *newMat, double *oldMat, int rOld, int cOld, int r1, int c1, int r2, int c2)
{
	int j, i, k=0;

	for(j=c1; j < (c2+1); j++){
		for(i=r1; i < (r2+1); i++){
			newMat[k] = oldMat[j*rOld + i];
			k++;
		}
	}
}

void  printType(SEXP x){
	switch(TYPEOF(x)){
	case LGLSXP:
		debugMsg("LGLSXP\n"); break;
	case INTSXP:
		debugMsg("INTSXP\n"); break;
	case REALSXP:
		debugMsg("REALSXP\n"); break;
	case CPLXSXP:
		debugMsg("CPLSXP\n"); break;
	case STRSXP:
		debugMsg("STRSXP\n"); break;
	case VECSXP:
		debugMsg("VECSXP\n"); break;
	case LISTSXP:
		debugMsg("LISTSXP\n"); break;
	case DOTSXP:
		debugMsg("DOTSXP\n"); break;
	case NILSXP:
		debugMsg("NILSXP\n"); break;
	case SYMSXP:
		debugMsg("SYMSXP\n"); break;
	case CLOSXP:
		debugMsg("CLOSXP\n"); break;
	case ENVSXP:
		debugMsg("ENVSXP\n"); break;
	default:
		debugMsg("unknown\n"); break;
	};
}
*/

#endif /* MYUTILS_H_ */
