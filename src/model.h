//
//  model.h
//  
//
//  Created by Charley Chen on 7/31/15.
//
//

#ifndef _model_h
#define _model_h

#include "stdafx.h"
#include "MersenneTwister.h"

void Init(int B, int G, double *para, double ga, double gb, double sigma_square,MTRand & mtr );
void ReformatData(double *data, int N, int &B, int &G, double *new_data);
int Sample(double *p, int n, double d);
void UpperTruncatedPoissonPDF(double *lambda, double *p, int N, int n);
void CalculateLambda(double *data, double *Beta, double *Gamma, double *Lambda, int n);
int SampleBG(double *data, double *BETA, double *GAMMA, double MU, double PHI, int which, double scale, int upperlimit, int Acc, int n, int B, int G, double *BG, MTRand & mtr);
void SampleMuPhi(double *Para, int N, double PHI, double sigma_square, double ga, double gb, double *pMu, double *pPhi, MTRand & mtr);
double* DoModel(double *data, double *new_data, int n, int niters, int burnin, int stride,int& betaAcc, int& gammaAcc, int& numberofmodels, int& B, int& G, MTRand& mt, int in_upper);
void ImportanceSampling(double *lambdas, double *data, double *yobs, double *result,int L, int m, int n, int upperlimit);
void SampleLambdas(double *lambdas, double *result, int m, int n, int upperlimit, MTRand & mtr);
void CalculateRiskSummary(double *data, int n, double *riskmeasure, double *result, int in_upper);

#endif
