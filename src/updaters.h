/*
 *
 *  Created on: 24.04.2009
 *      Author: fabians
 *
  */

#ifndef UPDATERS_H_
#define UPDATERS_H_

#include "utils.h"


/* #############################################*/
/* define update function for rescaled/pre-multiplied Design Matrices */
static void updateXAlpha(double *XAlpha, double *X, double *G, double *ksi, int *q, int *qKsiUpdate, int *p, int *n){
	//calculates X*diag(ksi)*G
	int i,j, nqNoUpdate=(*q - *qKsiUpdate) * (*n), nrv=1;
	double one=1.0, zero=0.0;
	double *temp = R_Calloc(*q * (*n), double);

	//copy first qKsiNoUpdate columns of X to temp
	if(nqNoUpdate > 0) F77_CALL(dcopy)(&nqNoUpdate, X, &nrv, temp, &nrv);
	//scale remaining columns of X by ksi,write to temp
	for(j=*q-*qKsiUpdate; j < *q; j++){//cols
		for(i=0; i < *n; i++){ //rows
			temp[j*(*n) + i] = ksi[j] * X[j*(*n) + i];
		}
	}
	//XAlpha = temp %*%G = (X*ksi) %*% G
	F77_CALL(dgemm)("N","N", n, p, q, &one, temp, n, G, q, &zero, XAlpha, n FCONE FCONE);
	R_Free(temp);
}

static void updateXKsi(double *XKsiUpdate, double *X, double *alphaLong, int *q, int *qKsiNoUpdate, int *n){
	//calculates X*diag(G * alpha) for columns with variable ksi
	int i,j;
	for(j=*qKsiNoUpdate; j < *q; j++){
		for(i=0; i < *n; i++){
			XKsiUpdate[(j-*qKsiNoUpdate)*(*n) + i] = alphaLong[j] * X[j*(*n) + i];
		}
	}
}

static void initializeBlocksQR(XBlockQR *XBlocks, double *X,
		int n, int blocks, int lengthCoef, double *varCoef, double* scale
		){
	int j, i, l, k, rows;
	//fill up Xa, Xi
	for(j=0; j < blocks; j++){
		l = 0;
		rows= n + (XBlocks[j]).qA;
		for(i=(XBlocks[j]).indA1; i < (XBlocks[j]).indA2 + 1; i++){ //cols
			for(k=0; k < n; k++){	//rows
				(XBlocks[j]).Xa[l] = X[i*n + k] * scale[0]; //upper n x p is X/sqrt(sigma)
				l++;
				if(k == n-1){ //lower p x p is diag(1/sqrt(varCoef))
					(XBlocks[j]).Xa[l + (i - (XBlocks[j]).indA1)] = R_pow(varCoef[i], -0.5);
					l = (i+1-(XBlocks[j]).indA1)*rows; //set l to first entry of next column
				}
			} //end for k
		} //end for i
		if(blocks > 1){
			removeSubmatrix((XBlocks[j]).Xi, X, n, lengthCoef, 0, (XBlocks[j]).indA1, n-1, (XBlocks[j]).indA2);
		}
		(XBlocks[j]).reject=0;
	}
}

/* #############################################*/
/* define update functions for parameter blocks */

static void updateVarKsi(double *ksiUpdate, double *varKsi, double *ksiDF, double *mKsi, int qKsiNoUpdate, int q)
/* update varKsi ~ 1/Ga(df/2, df/2) (--> Ksi ~ t(df)) */
{
	double shape, rate;
	for(int j = 0; j < q-qKsiNoUpdate; j++){
		shape = *ksiDF/2 + 0.5;
		rate = 0.5 * ( (*ksiDF) + R_pow(ksiUpdate[j] - mKsi[j], 2));
		varKsi[j] = 1/rgamma(shape, rate);
	}
	return;
}

static void updateTau(double *penAlphaSq, double *gamma, double *tau2, double a1, double a2, int pPen)
{
	int j=0;
	for(j=0; j < pPen; j++){
		tau2[j] = 1/rgamma((a1 + 0.5), 1/(a2 + (penAlphaSq[j] / (2*gamma[j]))));
	}

	return;
}

static void updateP1Gamma(double *penAlphaSq, double *tau2, double *p1, double *gamma, double v0, double w, int pPen)
{

	double A;
	double c1 = log(w/(1-w)) + log(v0)* 1/2;
	double c2 = (1 - v0)/(2*v0);
	for(int j=0; j < pPen; j++){
		A = exp( c1  +  c2 * penAlphaSq[j]/tau2[j]);
		if(!ISNAN(A)){
			if(A > 10000){
				p1[j] = 1;
				gamma[j] = 1;
			}  else {
				if(A < 0.00001) {
					p1[j] = 0;
					gamma[j] = v0;
				} else {
					p1[j] =  A / (1 + A);
					gamma[j] = (runif(0,1) < p1[j]) ? 1.0 : v0;
				}
			}
		} // end if(A==NA)
	}
	return;
}

static void  updateCoefBlockQR(int j, int n, double *X, double *y,
		double *coef, XBlockQR *Xblocks, double sqrtInvSigma2,
		double *varCoef, int lengthCoef, int blocks, int nrv, int changeMean, double *mCoef, int scaleMode)
{
	int i=0, k=0, l=0, info, oneInt = 1;
	int qA = (Xblocks[j]).qA, qI = (Xblocks[j]).qI, rows = n + qA;
	double minusOne = -1.0, one=1.0;

	//fill up Xa, Xi
	for(i=(Xblocks[j]).indA1; i < (Xblocks[j]).indA2 + 1; i++){ //cols
		for(k=0; k < n; k++){	//rows
			(Xblocks[j]).Xa[l] = X[i*n + k] * sqrtInvSigma2; //upper n x p is X/sqrt(sigma)
			l++;
			if(k == n-1){ //lower p x p is diag(1/sqrt(varCoef))
				(Xblocks[j]).Xa[l + (i - (Xblocks[j]).indA1)] = R_pow(varCoef[i], -0.5);
				l = (i+1-(Xblocks[j]).indA1)*rows; //set l to first entry of next column
			}
		} //end for k
	} //end for i
	if(blocks > 1){
		removeSubmatrix((Xblocks[j]).Xi, X, n, lengthCoef, 0, (Xblocks[j]).indA1, n-1, (Xblocks[j]).indA2);
	}

	/*ya = c((y - Xi*coefI)/sqrt(sigma), rep(0,qA))*/
	F77_CALL(dcopy)(&n, y, &nrv, (Xblocks[j]).ya, &nrv);
	if(blocks > 1){
	/*set (Xblocks[j]).ya = y - (Xblocks[j]).Xi *coefI*/
		l = 0;
		for(k = 0; k < (Xblocks[j]).indA1; k++){
			(Xblocks[j]).coefI[l] = coef[k];
			l++;
		}
		for(k = (Xblocks[j]).indA2 + 1; k < lengthCoef; k++){
			(Xblocks[j]).coefI[l] = coef[k];
			l++;
		}
		//ya = (y - X_i * beta_i)/sqrt(sigma2)
		F77_CALL(dgemm)("N", "N", &n, &oneInt, &qI, &minusOne,  (Xblocks[j]).Xi, &n, (Xblocks[j]).coefI, &qI, &one, (Xblocks[j]).ya, &n FCONE FCONE);
	}
	F77_CALL(dscal)(&n, &sqrtInvSigma2, (Xblocks[j]).ya, &oneInt);

	/*do QR*/
	F77_CALL(dqrdc)((Xblocks[j]).Xa, &rows, &rows, &qA, (Xblocks[j]).qraux, (Xblocks[j]).pivots, (Xblocks[j]).work, &oneInt);

	/*set prior mean to zero if changeMean=0 else to -1 or 1 depending on previous sign of coef*/
	for(k = 0; k < qA; k++){
		(Xblocks[j]).ya[n+k] = (changeMean) ?  mCoef[(Xblocks[j]).indA1 + (Xblocks[j]).pivots[k]-1]  : 0.0;
	}

	/*get coefs, write into (Xblocks[j]).m - these are permuted , see pivots !*/
	F77_CALL(dqrcf)((Xblocks[j]).Xa, &rows, &qA, (Xblocks[j]).qraux, (Xblocks[j]).ya, &oneInt, (Xblocks[j]).m, &info);

	/*get eps*/
	 for(k = 0; k < qA; k++)
	        (Xblocks[j]).err[k] = rnorm(0, 1);
  	//solve R' * err = (Xblocks[j]).err for err (<=> err  = R'^-1 (Xblocks[j]).err => (Xblocks[j]).err ~ N(0, R^-1 R'^-1) = N(0, (R' R)^-1 )
	F77_CALL(dtrtrs)("U", "T", "N", &qA, &oneInt, (Xblocks[j]).Xa, &rows, (Xblocks[j]).err, &qA, &info FCONE FCONE FCONE);

	// un-permute draw and write to coef
	for(k = 0; k < qA; k++){
	    coef[(Xblocks[j]).indA1 + (Xblocks[j]).pivots[k]-1] = (Xblocks[j]).m[k] + (Xblocks[j]).err[k];
	}
}

static void  updateCoefBlockQR_IWLS(int j, int n, double *X, double *y,
		double *coef, XBlockQR *Xblocks, double* scale,
		double *varCoef, int lengthCoef, int blocks, int nrv, int changeMean, double *mCoef,
		int family, // 0=normal, 1=binom, 2=poisson, ....
		double *modeCoef,
		double *eta,
		double *accept,
		double *offset,
		double modeSwitching,
		int scaleMode)
{
	//PhD-Thesis A. Brezger: p.57,  sampling scheme 1
	int i, k, l, info=0, oneInt = 1, switchMode = 0;
	int qA = (Xblocks[j]).qA, qI = (Xblocks[j]).qI, rows = n + qA;
	double minusOne = -1.0, one=1.0, likC, logAccProb, proposalDens;
	double *coefC =  R_Calloc(qA, double);
	double *coefChange = R_Calloc(qA, double);
	double *var = R_Calloc(n, double);
	double *w = R_Calloc(n, double);
	double *mu = R_Calloc(n, double);
	double *etaOffset = R_Calloc(n, double);
	double *etaMode= R_Calloc(n, double);

	//0. update Xa, Xi (see updateBlocks)
	l = 0;
	for(i=(Xblocks[j]).indA1; i < (Xblocks[j]).indA2 + 1; i++){ //cols
		for(k=0; k < n; k++){	//rows
			(Xblocks[j]).Xa[l] = X[i*n + k]; //upper n x p is X
			l++;
			if(k == n-1){ //lower p x p is diag(1/sqrt(varCoef))
				(Xblocks[j]).Xa[l + (i - (Xblocks[j]).indA1)] = R_pow(varCoef[i], -0.5);
				l = (i+1-(Xblocks[j]).indA1)*rows; //set l to first entry of next column
			}
		} //end for k
	} //end for i
	if(blocks > 1){
		removeSubmatrix((Xblocks[j]).Xi, X, n, lengthCoef, 0, (Xblocks[j]).indA1, n-1, (Xblocks[j]).indA2);
	}

	//1. compute lik for current coef
	//update etaOffset = offset + X*coefP:
	for(i=0; i<n; i++) etaOffset[i] = eta[i] + offset[i];
	likC = logLik(y, etaOffset, family, scale, n);
	//printDouble("likC  ", &likC, 1);
	//printDouble("eta   ", eta, n);
	//printDouble("offset", offset, n);

	//2. replace coef with previous mode
	//set mode of proposal to current coef every now and then to avoid getting stuck
	if(runif(0,1) < modeSwitching){
		switchMode = 1;
		//debugMsg("switchMode!");
	}
	for(k = 0; k < qA; k++){
			coefC[k] = 	coef[(Xblocks[j]).indA1 + k];
		    coef[(Xblocks[j]).indA1 + k] = modeCoef[(Xblocks[j]).indA1 + k];
		    coefChange[k] = coefC[k] - modeCoef[(Xblocks[j]).indA1 + k];
	}
	//printDouble("curr", coefC, qA);
	//printDouble("mode", &modeCoef[(Xblocks[j]).indA1], qA);

	//3. draw prop:
	//3 a) compute working w and working y based on previous mode
	//update etaMode= eta -X_a(coefC - modeCoef)
	F77_CALL(dcopy)(&n, eta, &oneInt, etaMode, &oneInt);
	F77_CALL(dgemm)("N", "N", &n, &oneInt, &qA, &minusOne,
			&(X[(Xblocks[j]).indA1*n]), &n, coefChange, &qA, &one, etaMode, &n FCONE FCONE);
	//printDouble("etaMode", etaMode, n);
	//compute mu, var, w = sqrt(var)
	if(family==1) { //binomial
		 for(k=0; k<n; k++){
			 mu[k] = 1 / (1 + exp(-(etaMode[k] + offset[k])));
			 var[k] = fmax2(0.00001, mu[k] * (1 - mu[k]) * scale[k]);
			 w[k] = R_pow(var[k], 0.5);
		 }
	}
	if(family==2) { //poisson
			 for(k=0; k<n; k++){
				 mu[k] = fmin2(100000, exp(etaMode[k] + offset[k]));
				 w[k] = R_pow(mu[k], 0.5);
				 var[k] = fmax2(0.000001, mu[k]);
			 }
	}

	//compute working response ya = w(k)* working reponse
	if(blocks > 1){
		for(k=0; k<n; k++){
				(Xblocks[j]).ya[k] = etaMode[k] + (y[k]-mu[k])/var[k];
		}

		/*set (Xblocks[j]).ya = (Xblocks[j]).ya - (Xblocks[j]).Xi *coefI*/
			l = 0;
			for(k = 0; k < (Xblocks[j]).indA1; k++){
				(Xblocks[j]).coefI[l] = coef[k];
				l++;
			}
			for(k = (Xblocks[j]).indA2 + 1; k < lengthCoef; k++){
				(Xblocks[j]).coefI[l] = coef[k];
				l++;
			}
			//ya = (ya - X_i * beta_i)
			F77_CALL(dgemm)("N", "N", &n, &oneInt, &qI, &minusOne,  (Xblocks[j]).Xi, &n, (Xblocks[j]).coefI, &qI, &one, (Xblocks[j]).ya, &n FCONE FCONE);
			//ya = (ya - X_i * beta_i)*w
			for(k=0; k<n; k++){
					(Xblocks[j]).ya[k] = (Xblocks[j]).ya[k] * w[k];
			}
	} else {
		for(k=0; k<n; k++){
			(Xblocks[j]).ya[k] = (etaMode[k] + (y[k]-mu[k])/var[k]) * w[k];
		}
	}
	//printDouble("ya*w   ", (Xblocks[j]).ya, n);

	//3 b) do QR
	//Xa=Xa*w
	l = 0;
	for(i=0; i < qA; i++){ //cols
			for(k=0; k < n; k++){	//rows
				(Xblocks[j]).Xa[i*rows + k] = (Xblocks[j]).Xa[i*rows + k]* w[k]; //upper n x qA is Xa / sqrt(weight)
			} //end for k
	} //end for i

	/*do QR (with pivoting) : upper left qA*qA contains triangular R*/
	F77_CALL(dqrdc)((Xblocks[j]).Xa, &rows, &rows, &qA, (Xblocks[j]).qraux, (Xblocks[j]).pivots, (Xblocks[j]).work, &oneInt);
	/*set prior mean to zero if changeMean=0 else to mCoef*/
	for(k = 0; k < qA; k++){
				(Xblocks[j]).ya[n+k] = (changeMean) ?  mCoef[(Xblocks[j]).indA1 + (Xblocks[j]).pivots[k]-1]*R_pow(varCoef[(Xblocks[j]).indA1 + (Xblocks[j]).pivots[k]-1], -0.5)  : 0.0;
	}

	//3 c) draw prop
	/*get mode, write into (Xblocks[j]).m - these are permuted , see pivots !*/
	F77_CALL(dqrcf)((Xblocks[j]).Xa, &rows, &qA, (Xblocks[j]).qraux, (Xblocks[j]).ya, &oneInt, (Xblocks[j]).m, &info);
	//if(info != 0) debugMsg("		dqrcf error!\n");

	if(anyNA((Xblocks[j]).m, qA)){
		//	Error handling??
	}
	//printDouble("mode   ", (Xblocks[j]).m, qA);

	if( (switchMode > 0) || anyNA((Xblocks[j]).m, qA)>0){
		switchMode = 1;
		for(k = 0; k < qA; k++){
				(Xblocks[j]).m[k]  = coefC[(Xblocks[j]).pivots[k]-1];
			}
	}

	/*get eps*/
	for(k = 0; k < qA; k++) (Xblocks[j]).err[k] = rnorm(0, 1);

	//solve R * err = (Xblocks[j]).err for err (<=> err  = R^-1 (Xblocks[j]).err => (Xblocks[j]).err ~ N(0, R^-1 (R^-1)') = N(0, (R' R)^-1 )
	F77_CALL(dtrtrs)("U", "N", "N", &qA, &oneInt, (Xblocks[j]).Xa, &rows, (Xblocks[j]).err, &qA, &info FCONE FCONE FCONE);
	if(info != 0) //FIXME

	//4. exchange coef with prop
	// un-permute draw and write to coef
	for(k = 0; k < qA; k++){
		    coef[(Xblocks[j]).indA1 + (Xblocks[j]).pivots[k]-1] = (Xblocks[j]).m[k] + (Xblocks[j]).err[k];
	}
	//printDouble("prop   ", &coef[(Xblocks[j]).indA1], qA);

	//rescale coef, (Xblocks[j]).m, (Xblocks[j]).err accordingly
	// R*(prop - (Xblocks[j]).m)
	if(switchMode<1){
		for(k = 0; k < qA; k++){
			l = qA-k;
			coefChange[k] = F77_CALL(ddot)(&l, &((Xblocks[j]).Xa[k*rows + k]), &rows, &((Xblocks[j]).err[k]), &oneInt);
		}
		proposalDens = -.5 * F77_CALL(ddot)(&qA, coefChange, &oneInt, coefChange, &oneInt);
	} else {
		//if modeSwitch was done, the ratio of proposals is \phi(current; mu=proposal, Sigma=R'R)/ \phi(proposal; current, R'R),
		// which is 1
		proposalDens = 0;
	}
	//printDouble("propDens:", &proposalDens, 1);


	//5. compute lik-ratio,  prior ratio, proposal ratio
	//update etaOffset = offset + X*coefP:
	F77_CALL(dcopy)(&n, offset, &oneInt, etaOffset, &oneInt);
	F77_CALL(dgemm)("N", "N", &n, &oneInt, &lengthCoef, &one,
						X, &n, coef, &lengthCoef, &one, etaOffset, &n FCONE FCONE);
	// log(acc prob) = likP - likC  ...
	logAccProb = logLik(y, etaOffset, family, scale, n) - likC;
	//printDouble("likP - likC", &logAccProb, 1);
	//printDouble("offset   ", offset, n);
	//printDouble("etaOffset", etaOffset, n);
	//printDouble("eta      ", eta, n);

	for(k = 0; k < qA; k++){
		logAccProb += 0.5 * (((R_pow(coefC[k]-mCoef[(Xblocks[j]).indA1 + k], 2.0) -
				R_pow(coef[(Xblocks[j]).indA1 + k] - mCoef[(Xblocks[j]).indA1 + k], 2.0))/
				varCoef[(Xblocks[j]).indA1 + k]));
	}
	//printDouble("likP-likC+logPriorRatio", &logAccProb, 1);

	if(switchMode<1){
		//+ propC - propP: propC ~ N((Xblocks[j]).m, (R' R)^-1); propP ~ N((Xblocks[j]).m, (R' R)^-1)
		// -.5( (current - (Xblocks[j]).m)'(R' R)(current - (Xblocks[j]).m) -  (prop - (Xblocks[j]).m)'(R' R)(prop - (Xblocks[j]).m)):
		//	.5 * (prop - (Xblocks[j]).m)'(R' R)(prop - (Xblocks[j]).m) = .5 * sum(R%*%(Xblocks[j]).err)^2 (added to logAccProb above)
		// calculate diff to current permuted like R from dqrdc
		for(k = 0; k < qA; k++){
			coefChange[k] = coefC[(Xblocks[j]).pivots[k]-1] - (Xblocks[j]).m[k];
		}
		// R*(current - (Xblocks[j]).m)
		for(k = 0; k < qA; k++){
			l = qA-k;
			coefChange[k] = F77_CALL(ddot)(&l, &((Xblocks[j]).Xa[k*rows + k]), &rows, &(coefChange[k]), &oneInt);
		}

		logAccProb += - .5 * F77_CALL(ddot)(&qA, coefChange, &oneInt, coefChange, &oneInt) - proposalDens;
	}// if Mode was switched to current, the proposal ratio is 1: \phi(current; proposal, R'R) = \phi(proposal; current, R'R)
	//printDouble("likP-likC+logPriorRatio+propRatio", &logAccProb, 1);

	//6. accept(+ update eta)/ reject( + replace prop with current coef if reject), update mode
	if(log(runif(0,1)) <  logAccProb){
		//accept
		accept[j] = accept[j] + 1.0;
		for(i=0; i<n; i++) eta[i] = etaOffset[i] - offset[i];
		printDouble("eta      ", eta, n);
		(Xblocks[j]).reject = 0;

	} else {
		//reject
		for(k = 0; k < qA; k++){
			coef[(Xblocks[j]).indA1 + k] = coefC[k];
		}
		(Xblocks[j]).reject ++;
	}
	if(anyNA((Xblocks[j]).m, qA) == 0){
		for(k = 0; k < qA; k++){
			modeCoef[(Xblocks[j]).indA1 + (Xblocks[j]).pivots[k]-1] = (Xblocks[j]).m[k];
		}
	}

	R_Free(etaMode);
	R_Free(etaOffset);
	R_Free(mu);
	R_Free(w);
	R_Free(var);
	R_Free(coefChange);
	R_Free(coefC);
}

static  void updateCoefQR(
		double *y, double *XCoef,
		XBlockQR *CoefBlocks,
		int blocksCoef,
		double *coef,
		double *varCoef, int lengthCoef,
		double *scale,
		int n, int nrv,int oneInt, int info, double minusOne, double zero, double one, int changeMean, double *priorMeanCoef,
		int family, // 0=normal, 1=binom, 2=poisson, ....
		double *modeCoef,
		double *eta,
		double *accept,
		double *offset,
		double modeSwitching,
		int scaleMode
){
// FIXME: DIRTYDIRTY- scale[0] is 1/sqrt(sigma2) for gaussian, else the real scale parameter (e.g. n_i for binomial)
	int j, k, length=blocksCoef;
	double *yOffset = R_Calloc(n, double);
	int *shuffle = R_Calloc(blocksCoef, int);
	int *iterator = R_Calloc(blocksCoef, int);

	if(family == 0) { // offset for IWLS incorporated into determining working obs, only need this for non-IWLS
		for(j = 0; j < n; j++) yOffset[j] = y[j] - offset[j];
	}

	for(k=0;  k < blocksCoef; k++) iterator[k] = k;
	for(k=0;  k < blocksCoef; k++){
		j = length * runif(0,1);
		shuffle[k] = iterator[j];
		iterator[j] = iterator[--length];
	}

	for(k=0; k < blocksCoef; k++){
		j = shuffle[k];

		if(family == 0) {
			updateCoefBlockQR(j, n, XCoef, yOffset, coef, CoefBlocks,
					scale[0],  varCoef, lengthCoef, blocksCoef, nrv, changeMean, priorMeanCoef, scaleMode);
		} else {
			updateCoefBlockQR_IWLS(j, n, XCoef, y, coef, CoefBlocks,
								scale,  varCoef, lengthCoef, blocksCoef, nrv, changeMean, priorMeanCoef,
								family, modeCoef, eta, accept, offset, modeSwitching, scaleMode);
		}

	}
	R_Free(iterator);
	R_Free(shuffle);
	R_Free(yOffset);
}

static  void  rescaleKsiAlpha(double *ksi, double *alpha, double *varKsi, double *tau2, double *G, int *d, int p, int q,
		int qKsiNoUpdate, int pPen, int scaleMode, double *modeAlpha, double *modeKsi, int family){
	/*rescale ksi groupwise so that sum(abs(ksi_g))/dim(g)=1, change alphas accordingly*/
	/* also tried rescaling by 1/max(abs(ksi_g)) --> stronger regularization*/
	int  i ,j, oneInt =1;
	double *scale = R_Calloc(p, double);
		setToZero(scale, p);
	double *scaleLong = R_Calloc(q, double);
		setToZero(scaleLong, q);
	double one = 1.0, zero=0.0;

	if(scaleMode>0){
		for(i = 0; i < p; i++){
				//step through (lower part of) i^th column of G, compute scale factor
				for(j = i * q + qKsiNoUpdate; j < (i + 1) * q; j++){
					if(G[j] > 0.0){
						if(scaleMode == 1){
							//rescale s.t. mean(abs(ksi_g)) = 1
							scale[i] += sign(ksi[j - i * q]) * ksi[j - i * q];
						} else {
							if(scaleMode==2){
								//rescale s.t. max(abs(ksi_g)) = 1
								scale[i] = fmax2(scale[i], sign(ksi[j - i * q]) * ksi[j - i * q]);
							} else {
								scale[i] = 1.0;
							}
						}
					}
				}
				if(scale[i] == 0.0) scale[i] = 1.0; //set scale for ksis that aren't updated to 1
				scale[i] = d[i]/scale[i];
				//add change to mode s.t. "mode" stays reasonably close to current coefficent otherwise chain gets stuck
				if(family != 0) modeAlpha[i] = modeAlpha[i] + (alpha[i]/scale[i] - alpha[i]);
				alpha[i] = alpha[i]/scale[i];
				//tau2[i-pPen] = tau2[i-pPen]*R_pow(scale[i], -2.0);
		}


		F77_CALL(dgemm)("N","N", &q, &oneInt, &p, &one, G, &q, scale, &p, &zero, scaleLong, &q FCONE FCONE);

		for(j = qKsiNoUpdate; j < q; j++) {
			//add change to mode s.t. "mode" stays reasonably close to current coefficent otherwise chain gets stuck
			if(family != 0) modeKsi[j - qKsiNoUpdate] = modeKsi[j - qKsiNoUpdate] + (ksi[j]*scaleLong[j] - ksi[j]);
			ksi[j] *= scaleLong[j];
			//varKsi[j-qKsiNoUpdate] *= R_pow(scaleLong[j], 2.0);
		}
	}
	R_Free(scaleLong);
	R_Free(scale);
}

static double updateLogPost(
		double *y,
		double *alpha, double *varAlpha,
		double *ksi, double *varKsi,
		double *scale, double b1, double b2,
		double *gamma, double w, double alphaW, double betaW,
		double *tau2, double a1, double a2,
		int n, int q, int p, int pPen, int pIncluded,
		int qKsiNoUpdate,
		double* priorMeanKsi,
		int family,
		double lik)
{
	double lp = 0.0,
		gammaPost = 0.0,
		tau2Post = 0.0,
		alphaPost=0.0,
		ksiPost=0.0,
		sigma2Post=0.0;
	int j=0;
	//likelihood
	//sum(dbinom(gamma == 1, size=1, prob=w, log=TRUE)) + # gamma
	gammaPost = pIncluded * log(w) + (pPen-pIncluded) * log(1-w);
	//sum(- (a1 + 1) * log(tau2) - (a2 / tau2)) + # tau2
	for(j=0; j < pPen; j++) tau2Post += -(a1 + 1) * log(tau2[j]) - (a2 / tau2[j]);
	tau2Post += pPen * (a1 * log(a2) -lgamma(a1));
	for(j=0; j < p; j++)  alphaPost += dnorm(alpha[j], 0, R_pow(varAlpha[j], 0.5), 1);
	for(j=qKsiNoUpdate; j < q; j++)  ksiPost += dnorm(ksi[j], priorMeanKsi[j-qKsiNoUpdate], R_pow(varKsi[j-qKsiNoUpdate], 0.5), 1);
	if(family==0){
		sigma2Post = (- (b1 + 1) * log(R_pow(scale[0], -2.0)) - (b2 / R_pow(scale[0], -2.0)) + b1*log(b2) - lgamma(b1));
	}
	lp = lik + gammaPost + tau2Post + alphaPost + ksiPost + sigma2Post + (dbeta(w, alphaW, betaW, 1));
	return(lp);
}

#endif /* UPDATERS_H_ */
