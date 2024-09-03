/*
 * sampler.c
 *
 *  Created on: 27.04.2009
 *      Author: fabians
 *
 *
 */

/*   ****************************************
 *
 *
 *    SAMPLER
 *
 *
 * ****************************************** */
#define USE_FC_LEN_T
#include "utils.h"
#include "updaters.h"

//extern "C"{
SEXP sampler(
		/*prior params*/
		double *a1, double *a2, /* prior for tau2*/
		double *b1, double *b2, /* prior for sigma2 */
		double *alphaW, double *betaW, /* prior for w */
		double *v0, /* gamma */
		double *varKsi, /*vector length qKsiUpdate!!*/

		/*model dimensions*/
		int *q, /*length of ksi*/
		int *qKsiUpdate, /*length of updated ksi*/
		int *p,   /*length alpha*/
		int *pPen,   /*length penalized alpha/ tau2 / gamma*/
		int *n,   /* no. of  obs.*/
		int *d,   /*vector (length p): group sizes*/

		/*parameter vectors*/
		double *beta,
		double *alpha,
		double *ksi,
		double *tau2,
		double *gamma,
		double *sigma2,
		double *w,

		/* (precomputed) constants */
		double *y,
		double *X,
		double *G,
		double *scale,
		double *offset,

		/*info about updateBlocks*/
		int *blocksAlpha,
		int *indA1Alpha,
		int *indA2Alpha,

		int *blocksKsi,
		int *indA1Ksi,
		int *indA2Ksi,

		/*MCMC parameters*/
		int *pcts,
		int *burnin,
		int *thin,
		int *totalLength,
		int *verbose,
		double *ksiDF,
		int *scaleMode,
		double *modeSwitching,
		int *family,
		double *acceptKsi,
		double *acceptAlpha,

		/*return matrices*/
		double *betaMatR,
		double *alphaMatR,
		double *ksiMatR,
		double *gammaMatR,
		double *probV1MatR,
		double *tau2MatR,
		double *sigma2MatR,
		double *wMatR,
		double *likMatR,
		double *logPostMatR
)
{
	// ############################################### //
	// ######## unwrap/initialize args ############### //
	// ############################################### //
	int pIncluded=0, i=0, j=0, startPen = *p-*pPen, qKsiNoUpdate = *q - *qKsiUpdate,
			save = 0, keep = *burnin,  nrv =1,  info=0,
			nSamp=(*totalLength-*burnin)/(*thin), oneInt = 1, zeroInt = 0;

	double *p1 =R_Calloc(*pPen, double);

	double infV  = 100000, oneV = 1.0, zeroV = 0.0, minusOneV =-1.0;
	double *one=&oneV, *zero=&zeroV, *minusOne=&minusOneV, *inf=&infV, acceptance=0;
	double invSigma2 = 1 / *sigma2, sqrtInvSigma2 = R_pow(invSigma2, 0.5);
	double  *penAlphaSq, *alphaLong, *varAlpha, *priorMeanAlpha, *modeAlpha, *offsetAlpha;;
	penAlphaSq	= R_Calloc(*pPen, double);
	for(int i=*p-*pPen; i<*p; i++) penAlphaSq[i- *p + *pPen] = R_pow(alpha[i], 2.0);
	alphaLong = R_Calloc(*q, double);
	F77_CALL(dgemm)("N","N", q, &oneInt, p, one, G, q, alpha, p, zero, alphaLong, q FCONE FCONE);
	varAlpha = R_Calloc(*p, double);
	for(int i=0; i<startPen; i++) varAlpha[i] = *inf; /*unpenalized*/
	for(int i=startPen; i<*p; i++) varAlpha[i] = tau2[i-startPen]*gamma[i-startPen]; /*penalized*/
	priorMeanAlpha	= R_Calloc(*p, double);
	setToZero(priorMeanAlpha, *p);
	modeAlpha = R_Calloc(*p, double);
	F77_CALL(dcopy)(p, alpha, &oneInt, modeAlpha, &oneInt);
	offsetAlpha = R_Calloc(*n, double);
	F77_CALL(dcopy)(n, offset, &oneInt, offsetAlpha, &oneInt);


	double *ksiUpdate, *priorMeanKsi, *modeKsi,  *offsetKsi;
	int safeQKsiUpdate = imax2(1, *qKsiUpdate);
	//ksiUpdate contains the last qKsiUpdate elements in ksi
	ksiUpdate = R_Calloc(safeQKsiUpdate, double);
	F77_CALL(dcopy)(&safeQKsiUpdate, &ksi[*q-safeQKsiUpdate], &oneInt, ksiUpdate, &oneInt);
	priorMeanKsi = R_Calloc(safeQKsiUpdate, double);
	setToZero(priorMeanKsi, safeQKsiUpdate);
	for(int i=0; i<*qKsiUpdate; i++) priorMeanKsi[i] = 1.0;
	modeKsi = R_Calloc(safeQKsiUpdate, double);
	setToZero(modeKsi, safeQKsiUpdate);
	for(int i=0; i<*qKsiUpdate; i++) modeKsi[i] = ksi[i+qKsiNoUpdate];
	// offsetKsi = offset + X_d=1*alpha : use lin.predictor of grps with ksi==1 as offset
	offsetKsi = R_Calloc(*n, double);
	F77_CALL(dcopy)(n, offset, &oneInt, offsetKsi, &oneInt);
	if(qKsiNoUpdate < *q){
		if(qKsiNoUpdate > 0){
			F77_CALL(dgemm)("N","N", n, &oneInt, &qKsiNoUpdate, one, X, n, alpha, &qKsiNoUpdate, one, offsetKsi, n FCONE FCONE);
		}
	}

	double	*eta, *resid, rss, *XAlpha, *XKsiUpdate, *etaOffset;
	eta	= R_Calloc(*n, double);
	F77_CALL(dgemm)("N","N", n, &oneInt, q, one, X, n, beta, q, zero, eta, n FCONE FCONE);
	resid = R_Calloc(*n, double);
	rss = 0;
	for(int i=0; i<*n; i++) {
		resid[i] = y[i]-eta[i] - offset[i];
		rss += R_pow(resid[i], 2.0);
	}
	XAlpha = R_Calloc(*p * (*n), double);
	updateXAlpha(XAlpha, X, G, ksi, q, qKsiUpdate, p, n);
	XKsiUpdate = R_Calloc( *n * safeQKsiUpdate, double);
	setToZero(XKsiUpdate, *n * safeQKsiUpdate);
	if(qKsiNoUpdate < *q){
		updateXKsi(XKsiUpdate, X, alphaLong, q, &qKsiNoUpdate, n);
	}
	etaOffset	= R_Calloc(*n, double);
	for(int i=0; i<*n; i++) etaOffset[i] = eta[i]+offset[i];


	// ############################################################ //
	// ######## set up blocks for blockwise updates ############### //
	// ############################################################ //

	XBlockQR *AlphaBlocks = R_Calloc(*blocksAlpha, XBlockQR);
	XBlockQR *KsiBlocks = R_Calloc(*blocksKsi, XBlockQR);


	for(int i=0; i < *blocksAlpha; i++){
		(AlphaBlocks[i]).indA1 = indA1Alpha[i];
		(AlphaBlocks[i]).indA2 = indA2Alpha[i];

		(AlphaBlocks[i]).qA = (AlphaBlocks[i]).indA2 - (AlphaBlocks[i]).indA1 + 1;
		(AlphaBlocks[i]).qI = *p - (AlphaBlocks[i]).qA;

		(AlphaBlocks[i]).qraux = R_Calloc((AlphaBlocks[i]).qA, double);
		setToZero((AlphaBlocks[i]).qraux, (AlphaBlocks[i]).qA);
		(AlphaBlocks[i]).work = R_Calloc((AlphaBlocks[i]).qA, double);
		setToZero((AlphaBlocks[i]).work, (AlphaBlocks[i]).qA);
		(AlphaBlocks[i]).pivots = R_Calloc((AlphaBlocks[i]).qA, int);
		for(int j=0; j < (AlphaBlocks[i]).qA; j++) (AlphaBlocks[i]).pivots[j] = 0;

		(AlphaBlocks[i]).coefI = R_Calloc((AlphaBlocks[i]).qI, double);
		setToZero((AlphaBlocks[i]).coefI, (AlphaBlocks[i]).qI);

		(AlphaBlocks[i]).Xa = R_Calloc(((AlphaBlocks[i]).qA + *n) * (AlphaBlocks[i]).qA, double);
		setToZero((AlphaBlocks[i]).Xa, ((AlphaBlocks[i]).qA + *n) * (AlphaBlocks[i]).qA);
		(AlphaBlocks[i]).Xi = R_Calloc(*n * (AlphaBlocks[i]).qI, double);
		setToZero((AlphaBlocks[i]).Xi, *n * (AlphaBlocks[i]).qI );
		(AlphaBlocks[i]).ya = R_Calloc(((AlphaBlocks[i]).qA + *n), double);
		F77_CALL(dcopy)(n, y, &nrv, (AlphaBlocks[i]).ya, &nrv);
		setToZero((AlphaBlocks[i]).ya + *n, (AlphaBlocks[i]).qA);

		(AlphaBlocks[i]).m = R_Calloc((AlphaBlocks[i]).qA, double);
			setToZero((AlphaBlocks[i]).m, (AlphaBlocks[i]).qA);
		(AlphaBlocks[i]).err = R_Calloc((AlphaBlocks[i]).qA, double);
			setToZero((AlphaBlocks[i]).err, (AlphaBlocks[i]).qA);

	}
	initializeBlocksQR(AlphaBlocks, XAlpha, *n, *blocksAlpha, *p, varAlpha, scale);


	if(*qKsiUpdate > 0){
		for(int i=0; i < *blocksKsi; i++){
			(KsiBlocks[i]).indA1 = indA1Ksi[i];
			(KsiBlocks[i]).indA2 = indA2Ksi[i];

			(KsiBlocks[i]).qA = (KsiBlocks[i]).indA2 - (KsiBlocks[i]).indA1 + 1;
			(KsiBlocks[i]).qI = *qKsiUpdate - (KsiBlocks[i]).qA;

			(KsiBlocks[i]).qraux = R_Calloc((KsiBlocks[i]).qA, double);
			setToZero((KsiBlocks[i]).qraux, (KsiBlocks[i]).qA);
			(KsiBlocks[i]).work = R_Calloc((KsiBlocks[i]).qA, double);
			setToZero((KsiBlocks[i]).work, (KsiBlocks[i]).qA);
			(KsiBlocks[i]).pivots = R_Calloc((KsiBlocks[i]).qA, int);
			for(int j=0; j < (KsiBlocks[i]).qA; j++) (KsiBlocks[i]).pivots[j] = 0;

			(KsiBlocks[i]).coefI = R_Calloc((KsiBlocks[i]).qI, double);
			setToZero((KsiBlocks[i]).coefI, (KsiBlocks[i]).qI);

			(KsiBlocks[i]).Xa = R_Calloc(((KsiBlocks[i]).qA + *n) * (KsiBlocks[i]).qA, double);
			setToZero((KsiBlocks[i]).Xa, ((KsiBlocks[i]).qA + *n) * (KsiBlocks[i]).qA);
			(KsiBlocks[i]).Xi = R_Calloc(*n * (KsiBlocks[i]).qI, double);
			setToZero((KsiBlocks[i]).Xi, *n * (KsiBlocks[i]).qI );
			(KsiBlocks[i]).ya = R_Calloc(((KsiBlocks[i]).qA + *n), double);
			F77_CALL(dcopy)(n, y, &nrv, (KsiBlocks[i]).ya, &nrv);
			setToZero((KsiBlocks[i]).ya + *n, (KsiBlocks[i]).qA);

			(KsiBlocks[i]).m = R_Calloc((KsiBlocks[i]).qA, double);
			setToZero((KsiBlocks[i]).m, (KsiBlocks[i]).qA);
			(KsiBlocks[i]).err = R_Calloc((KsiBlocks[i]).qA, double);
			setToZero((KsiBlocks[i]).err, (KsiBlocks[i]).qA);
		}
		initializeBlocksQR(KsiBlocks, XKsiUpdate, *n, *blocksKsi, *qKsiUpdate, varKsi, scale);
	}

	// ############################################### //
	// ########     start sampling     ############### //
	// ############################################### //


#ifdef Win32
	R_FlushConsole();
#endif
	/* sampling */
	GetRNGstate();
	for(i = 0; i < *totalLength; i++)
	{
		debugMsg("\n###########################################\n\n");
		//update alpha
		{
			//update varAlpha
			for(j=startPen; j<*p; j++) varAlpha[j] = tau2[j-startPen] * gamma[j-startPen];
			//update alpha
			updateCoefQR(y, XAlpha, AlphaBlocks,
					*blocksAlpha,
					alpha,
					varAlpha, *p,
					scale,
					*n, nrv, oneInt, info, *minusOne, *zero, *one, 1, priorMeanAlpha,
					*family, modeAlpha, eta, acceptAlpha, offsetAlpha, *modeSwitching, zeroInt);
		}


		//update ksi
		if(qKsiNoUpdate < *q){

			//update alphaLong = G %*% alpha
			F77_CALL(dgemm)("N","N", q, &oneInt, p, one, G, q, alpha, p, zero, alphaLong, q FCONE FCONE);

			//update design for ksi
			updateXKsi(XKsiUpdate, X, alphaLong, q, &qKsiNoUpdate, n);

			//update offsetKsi
			if(qKsiNoUpdate > 0){
				F77_CALL(dcopy)(n, offset, &oneInt, offsetKsi, &oneInt);
				F77_CALL(dgemm)("N","N", n, &oneInt, &qKsiNoUpdate, one, X, n, alpha, &qKsiNoUpdate, one, offsetKsi, n FCONE FCONE);
			}

			for(j = 0; j < *qKsiUpdate; j++){
				priorMeanKsi[j] = sign(  1/(1 + exp(-2*ksiUpdate[j]/varKsi[j])) - runif(0,1) );
			}


			if(*ksiDF>0){
				updateVarKsi(ksiUpdate, varKsi, ksiDF, priorMeanKsi, qKsiNoUpdate, *q);
			}


			updateCoefQR(y, XKsiUpdate, KsiBlocks,
					*blocksKsi,
					ksiUpdate, varKsi, *qKsiUpdate,
					scale,
					*n, nrv, oneInt, info, *minusOne, *zero, *one, 1, priorMeanKsi,
					*family, modeKsi, eta, acceptKsi, offsetKsi, *modeSwitching, *scaleMode);
			//write back to ksi
			F77_CALL(dcopy)(qKsiUpdate, ksiUpdate, &oneInt, &ksi[*q-*qKsiUpdate], &oneInt);


			//rescale ksi, alpha & put back in ksiUpdate
			if(*scaleMode > 0){
				rescaleKsiAlpha(ksi, alpha, varKsi, tau2, G, d, *p, *q, qKsiNoUpdate, *pPen, *scaleMode, modeAlpha, modeKsi, *family);
				F77_CALL(dcopy)(qKsiUpdate, &ksi[*q-*qKsiUpdate], &oneInt, ksiUpdate, &oneInt);
			}

			//update XAlpha
			updateXAlpha(XAlpha, X, G, ksi, q, qKsiUpdate, p, n);

			//update alphaLong = G %*% alpha
			F77_CALL(dgemm)("N","N", q, &oneInt, p, one, G, q, alpha, p, zero, alphaLong, q FCONE FCONE);

		} else {
			F77_CALL(dcopy)(q, alpha, &oneInt, alphaLong, &oneInt);
		}

		for(int i = *p-*pPen; i < *p; i++) penAlphaSq[i - *p + *pPen] = R_pow(alpha[i], 2.0);
		updateTau(penAlphaSq, gamma, tau2, *a1, *a2, *pPen);

		updateP1Gamma(penAlphaSq, tau2, p1, gamma, *v0, *w, *pPen);
		pIncluded = 0;
		for(j=0; j<*p - startPen; j++) pIncluded += (gamma[j] == 1.0);

		*w = rbeta( *alphaW + pIncluded, *betaW + *p - pIncluded );

		// update beta
		for(j = 0; j < *q; j++){
			beta[j] = alphaLong[j]*ksi[j];
		}

		//update eta, eta+offset
		F77_CALL(dgemm)("N", "N", n, &oneInt, q, one, X, n, beta, q, zero, eta, n FCONE FCONE);
		for(int i=0; i<*n; i++) etaOffset[i] = eta[i] + offset[i];

		//update sigma_eps
		if(*family == 0){
			//resid = y - eta - offset
			F77_CALL(dcopy)(n, y, &nrv, resid, &nrv);  //resid <- y
			F77_CALL(daxpy)(n, minusOne, etaOffset, &nrv, resid, &nrv); //resid <- resid - eta - offset

			//rss = resid'resid
			rss = F77_CALL(ddot)(n, resid, &oneInt, resid, &oneInt);

			//update sigma2
			invSigma2 = rgamma(*n/2 + *b1, 1/(rss/2 + *b2));
			sqrtInvSigma2 = R_pow(invSigma2, 0.5);
			scale[0] = sqrtInvSigma2;
			*sigma2 = 1 / invSigma2;
		}


		if(i >= *burnin){
			/* report progress */
			if(*verbose){
				for(j=0; j<9; j++){
					if(i == pcts[j]){
						Rprintf(".");
						#ifdef Win32
							R_FlushConsole();
						#endif
						break;
					}
				}
			}
			/* save samples*/
			if(i == keep){
				for(j = 0; j < *q; j++){
					(betaMatR)[save + j*nSamp] = beta[j];
					(ksiMatR)[save + j*nSamp] = ksi[j];
				}
				for(j=0; j < *p; j++){
					(alphaMatR)[save + j*nSamp] = alpha[j];
				}
				for(j=0; j < *pPen; j++){
					(tau2MatR)[save + j*nSamp] = tau2[j];
					(gammaMatR)[save + j*nSamp] = gamma[j];
					(probV1MatR)[save + j*nSamp] = p1[j];
				}
				(wMatR)[save] = *w;
				(sigma2MatR)[save] = *sigma2;
				likMatR[save] = logLik(y, etaOffset, *family, scale, *n);
				(logPostMatR)[save] = updateLogPost(y, 	alpha, varAlpha,
						ksi, varKsi, scale, *b1, *b2, gamma, *w, *alphaW, *betaW,
						tau2, *a1, *a2,	*n, *q, *p, *pPen, pIncluded, qKsiNoUpdate, priorMeanKsi, *family, likMatR[save]);
				keep = keep + *thin;
				save ++;
				R_CheckUserInterrupt();
			}
		} else {
			if(*verbose){
				if(i == (*burnin-1)){
					Rprintf("b");
					#ifdef Win32
						R_FlushConsole();
					#endif
				}
			}
		}
	} /* end for i*/

	PutRNGstate();

	if(*verbose) Rprintf(".");
	if(*family > 0) {
		acceptance = 0.0;
		for(j=0; j<*blocksAlpha; j++) acceptance += acceptAlpha[j];
		acceptance = 0.0;
		if(qKsiNoUpdate < *q){
			for(j=0; j<*blocksKsi; j++) acceptance += acceptKsi[j];
		}
	}

	R_Free(etaOffset); R_Free(XKsiUpdate); R_Free(XAlpha);  R_Free(resid); R_Free(eta);
	R_Free(offsetKsi); R_Free(modeKsi); R_Free(priorMeanKsi);	R_Free(ksiUpdate);
	R_Free(offsetAlpha);
	R_Free(modeAlpha);
	R_Free(priorMeanAlpha);
	R_Free(varAlpha);
	R_Free(alphaLong);
	R_Free(penAlphaSq);
	R_FreeXBlockQR(AlphaBlocks, *blocksAlpha);
	if(qKsiNoUpdate < *q) R_FreeXBlockQR(KsiBlocks, *blocksKsi);
	R_Free(p1);
	return(R_NilValue);
}/* end sampler ()*/

