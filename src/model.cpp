#include "stdafx.h"
#include "MersenneTwister.h"
//#include "SpecialFunctions.h"
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;

//SpecialFunctions sf;

//The upper truncated limmit for Poisson
//#define UPPERLIMIT 8
//commented out since this is made an input

//declare the method
double gammarand(double a, double b, MTRand& mt);


/* MCMC initilization
	Input:
		B:				the unique number of source counties
		G:				the unique number of destination counties
		ga,ga:			gamma hyper-parameters
		sigma_square:	sigma_square hyper0paramtere
		mtr:			the reandom number generator
	Output:
		para:			the paramater vector to be initilized, at least B+G+4 in size, preallicated by caller
*/

void Init(int B, int G, double *para, double ga, 
        double gb, double sigma_square,MTRand & mtr ) {
    int i;
    para[0] = sqrt(sigma_square)*mtr.randNorm(0.0,1.0); //mu_beta
    para[1] = gammarand(ga, 1.0/gb,mtr);  //phi_beta
    para[2] = sqrt(sigma_square)*mtr.randNorm(0.0,1.0); //mu_gamma
    para[3] = gammarand(ga, 1.0/gb,mtr); //phi_gamma
    for (i = 4; i < 4+B; i++) {
        para[i] = para[0] + sqrt(para[1])*mtr.randNorm(0.0,1.0); //BETA
    }
    for (i = 4+B; i < 4+B+G; i++) {
        para[i] = para[2] + sqrt(para[3])*mtr.randNorm(0.0,1.0);//GAMMA
    }
}

/*	parse input data, and add unique indeix to each row for source/dest countied
	Input:
		data:	a N by 3 matrix for input data. Each row is a record, and each column is (source,dest,y)
		N:		the number of recrods
	Output:
		B:		the unique number of source counties
		G:		the unique number of dest counties
		new_data: the resulting N by 5 matrix, with two extra columns for unique source/dest county indexes.
				Caller is responsible for memory allocation.
*/
void ReformatData(double *data, int N, int &B, int &G, double *new_data) { //one based columns
    set<int> UniqueSource;
    set<int> UniqueDest; 
    map<int,int> SourceMap;
	map<int,int> DestMap;
	set<int>::iterator it;
	int i;
	int p = 3; //number of er o
    for (i = 0; i < N; i++) {
		//cout << data[i*p] << " " << data[i*p+1] << endl;
        UniqueSource.insert((int)data[i*p]);
        UniqueDest.insert((int)data[i*p+1]);
    }
	B = 0;
	for (it=UniqueSource.begin(); it!=UniqueSource.end(); ++it) {
		B++;
		SourceMap[*it] = B;
	}
	G = 0;
	for (it=UniqueDest.begin(); it!=UniqueDest.end(); ++it) {
		G++;
		DestMap[*it] = G;
	}
	int newcol = p + 2;
	for (i = 0; i < N; i++) {
		int base = i * newcol;
		new_data[base] = data[i *p];
		new_data[base+1] = data[i*p+1];
		new_data[base+2] = data[i*p+2];
		new_data[base+3] = SourceMap[(int)data[i *p]];
		new_data[base+4] = DestMap[(int)data[i *p+1]];
	}
}

/* draw from a multinomial distribution
	Input:
		p: a n by 1 vector of probabilities, (can be non-normalized)
		n: the length of p
		d: a uniform random number drawing outside of the function
	Output: the random sample (one based)
*/
int Sample(double *p, int n, double d) {
    double dsum;
    int i,k;
    dsum = 0;
    double *myw;
    myw = new double[n];
    for (i = 0; i < n;i++) {
        dsum+=p[i];
    }
    myw[0] = p[0] / dsum;
    for (i = 1; i < n;i++) {
        myw[i] = p[i] / dsum + myw[i-1];
    }
    
    for(k=0;k < n && d>myw[k];k++)
        ;
    delete [] myw;
    return k+1;
}

/*caclualte UpperTruncatedPoissonPDF for a list of lambda's
	Intput: 
		lambda:		the a N by 1 vector of lambdas
		N:			the number of lambda's in the vector
		n:			the upper truncated limit + 1 
	Output:			
		p:			a N by (n+1) matrix of normilized probabilities. 
					Caller is responsible for memory allocation/deallocation.
	Note: This function is vectorized for performance.

*/

void UpperTruncatedPoissonPDF(double *lambda, double *p, int N, int n) {
    for (int i = 0; i < N;i++) {
        int base = i * n;
        p[base] = 1.0;
        double dsum = p[base];
        int j;
        for (j = 1; j < n; j++) {
            p[base + j] = p[base + j-1] * (lambda[i] / double(j)); 
            dsum += p[base + j];
        }
        for (j = 0; j < n; j++) {
            p[base + j] = p[base + j] / dsum; 
        }
    }
} 

/* Caculate lambdas from data, beta and gamma,
	Input:
		data:	the n by 5 matrix from call to ReformatData 
		Beta:	the pointer to beta vector 
		Gamma:	the pointer to gamma vector
		n:		numer of observations
	Output:
		Lambda:	a n by 1 vector of lambdas. Caller is responsible for memory allocation.
*/
void CalculateLambda(double *data, double *Beta, double *Gamma, double *Lambda, int n) {
    int Ncol = 5;
    for (int i = 0; i < n; i++) {
        int bindex = (int)data[i * Ncol + 3] - 1; //from one based to zero based
        int gindex = (int)data[i * Ncol + 4] - 1; 
        Lambda[i] = exp(Beta[bindex] + Gamma[gindex]);
    }
}

/* sample beta and gamma in MCMC iteration
	Input:
		data:		the n by 5 matrix for formated data
		BETA:		the pointer to beta vector
		GAMMA:		the pointer to gamma vector
		Mu, PHI:	the mu,phi parameters
		which:		indicator for sampling beta or gamma.0 for beta, non-zero for gamma
		scale:		the sigma scale hyper-parameter
		upperlimt:  the truncated Poisson upperlimit
		Acc:		the previous acceptance number
		n:			the number of observations
		B:			the number of unique source counties
		G:			the number of unique dest counties
		mtr:		the random number generator
	Output:
		BG:			the sampled values for beta's or gamma's
		return:		the updated Acc value
*/
int SampleBG(double *data, double *BETA, double *GAMMA, double MU, double PHI, int which, double scale, int upperlimit, int Acc,
        int n, int B, int G, double *BG, MTRand & mtr) {
    
    double *lambda = new double[n];
    double *lambda_new = new double[n];
    double *p = new double[n * (upperlimit+1)];
    double *p_new = new double[n * (upperlimit+1)];
    double *logdiff = NULL;
    int i;
    int nCol = 5;
    
    CalculateLambda(data, BETA, GAMMA, lambda, n);
    UpperTruncatedPoissonPDF(lambda, p, n, upperlimit+1);
   
    int dataindex;
    int nBG;
    double *PARA = NULL;
    if (which == 0) {
        dataindex = 3;
        nBG = B;
        PARA = BETA;
    } else {
        dataindex = 4;
        nBG = G;
        PARA = GAMMA;
    }
    double *PARA_new = new double[nBG];
    for (i =0; i < nBG; i++) {
        PARA_new[i] = PARA[i] + scale * mtr.randNorm(0.0,1.0);
    }
    if (which == 0) {
        CalculateLambda(data, PARA_new, GAMMA, lambda_new, n);
    } else {
        CalculateLambda(data, BETA, PARA_new, lambda_new, n);
    }
    UpperTruncatedPoissonPDF(lambda_new, p_new, n, upperlimit+1);
    
    logdiff = new double[nBG];
    for (i =0; i < nBG; i++) {
        logdiff[i] =  -PHI*(PARA_new[i] - MU)* (PARA_new[i] - MU) / 2.0;
        logdiff[i] +=  PHI*(PARA[i] - MU)* (PARA[i] - MU) / 2.0; 
    }
    
    for (i = 0; i < n; i++) {
        int BGindex = (int)data[i * nCol + dataindex] - 1;
        int y  = (int)data[i * nCol + 2] - 1;
        logdiff[BGindex] +=  log(p_new[i*(upperlimit+1)+y]);
        logdiff[BGindex] -=  log(p[i*(upperlimit+1)+y]);
    }
    
    for (i =0; i < nBG; i++) {
        if (logdiff[i] > log(mtr.rand())) {
            BG[i] = PARA_new[i];
            Acc = Acc + 1;
        } else {
            BG[i] = PARA[i];
        }
    }
    
    delete [] logdiff;
    delete [] PARA_new;
    delete [] p;
    delete [] lambda_new;
    delete [] lambda;
    return Acc;
}

/*sample mu, and phi in MCMC iterations
	Input:		
		Para:			the beta/gamma values
		N:				number of beta/gamma;s
		PHI:			phi
		sigma_square:   the sigma_square hyper-parameter
		ga,gb:			the ga,gb hyper-parameters
		mtr:			the random number generator 
	Output:
		pMu:			the newly sampled mu
		pPhi:			the newly sampled phi
*/
void SampleMuPhi(double *Para, int N, double PHI, double sigma_square, double ga, double gb, double *pMu, double *pPhi, MTRand & mtr) {
    
	double mean_para = 0.0;
	int i;
	for (i=0; i < N; i++) {
		mean_para += Para[i];
	}
    double phi_star = double(N) * PHI + 1.0 / sigma_square;
    double mu_star  = double(N) * PHI * mean_para / phi_star;
    double Mu = mu_star +  mtr.randNorm(0.0,1.0) / sqrt(phi_star);
    
    double a = ga + double(N)/2.0;
	double dsum = 0.0;
	for (i =0; i < N; i++) {
		dsum += (Para[i] - mean_para) * (Para[i] - mean_para);
	}
    double b = gb + dsum / 2.0;
    double Phi = gammarand(a, 1.0/b,mtr);
	*pMu = Mu;
	*pPhi = Phi;
}

/* model fitting 
	Input:
		data:		the raw n by 3 input data 
		new_data:   the pointer to the formatted n by 5 data
		n:			number of observations
		mt:			the random number generator
		niters:		the total number of iterations
		burnin:		the total number of burnins
		stride:		save result for every stride number of iterations
	Output:
		betaAcc:		beta acceptance count
		gammaAcc:		gamma acceptance count
		numberofmodels: the number of models to be saved
		B:				the unique number of source counties
		G:				the unique number of dest counties
		return:			the saved model parameters from each saved iterations. caller is reponable for memory deallocation

*/
double* DoModel(double *data, double *new_data, int n,  
        int niters, int burnin, int stride, 
        int& betaAcc, int& gammaAcc, int& numberofmodels, int& B, int& G, MTRand& mt, int in_upper) {
    int i;
    
    
    double beta_c = 0.1; 
    double gamma_d = 0.1;;
    
    double sigma_square = 25.0;
    double ga = 0.01; 
    double gb = 0.01;
   
    int nsave = (niters - burnin) / stride + 2;

    ReformatData(data,n ,B,G,new_data);
    int npara = B + G + 4;

	//cout <<"npara = " << npara << endl;
    double *result = new double[npara * nsave];
    double  *current_para = result;
    Init(B, G, current_para, ga, gb,sigma_square,mt);
    
    int count = 0;
    numberofmodels = 1;
    
    betaAcc  = 0;
    gammaAcc = 0;
    
    for (int iter = 0; iter < niters; iter++) {
        count++;
        double *BETA = current_para + 4;
        double *GAMMA = BETA + B;
        betaAcc  += SampleBG(new_data, BETA, GAMMA, current_para[0], current_para[1], 0, beta_c,  in_upper, 0, n, B, G,BETA,mt);
        gammaAcc += SampleBG(new_data, BETA, GAMMA, current_para[2], current_para[3], 1, gamma_d, in_upper, 0, n, B, G,GAMMA,mt);
        SampleMuPhi(BETA, B, current_para[1], sigma_square, ga, gb, current_para+0, current_para+1, mt);
        SampleMuPhi(GAMMA,G, current_para[3], sigma_square, ga, gb, current_para+2, current_para+3, mt);
        
        if (count == stride) {
            if (iter >= burnin) {
                numberofmodels++;
                for (i =0; i< npara; i++) {
                    current_para[i+npara] = current_para[i];   
                }
                current_para += npara;
            }
            count = 0;
        }
    }
	numberofmodels--; //the last model was duplicated but not updated.
    return result;
}
/*
Input:
	lambdas:	A L by n matrix for Lambdas. Note: each row a lambda, easier for GPU later 
	data:		A m by n matrix for Synthetic data. 
	yobs:		A 1 by m vector for observed y values
	L:			Number of Lambdas used in importance sampling
	m:			Number fo synthetic data sets
	n:			Number of observations
	upperlimit:	always 8 in our case
Output:
	result:		A n by (upperlimit+1) matrix for risk measures. Caller is responable for allocating memory for result

*/
void ImportanceSampling(double *lambdas, double *data, double *yobs, 
        double *result,int L, int m, int n, int upperlimit) {
    int t,y,h,l,s,k;
    
    double *P1 = new double[L * n * (upperlimit+1)];
    UpperTruncatedPoissonPDF(lambdas, P1, L*n, upperlimit+1);
    
    double *P2 = new double[L* n * (upperlimit+1)];
    int SizeYL = (upperlimit+1) * L; 
    for (t = 0; t < n; t++) {
        for (y=1;y <= upperlimit+1; y++) {
            double yminusyt = y - yobs[t];
            double dsum = 0;
            for (h =0; h < L; h++) {
                int index = t * SizeYL + (y-1) * L + h;
                P2[index] = pow(lambdas[h*n+t],yminusyt);
                dsum += P2[index];
            }
            for (h =0; h < L; h++) {
                int index = t * SizeYL + (y-1) * L + h;
                P2[index] /= dsum;
            }
        }
    }
    
    //calculate the left-hand side of 1.4.5, which has 3 indexes (y,l,t) for
    //possible y values, saved lambda draws and each cell(obs)
    double *P_YMN = new double[(upperlimit+1)* m * n];
    for (y=1;y <= upperlimit+1; y++) {
        for (l=0; l < m; l++) {
            for (t=0; t< n; t++) { //outer loops
                
                //start of inner loop
                double dylt = 0.0;
                for (s=0; s < L; s++) {
                    double dprod = 1.0;
                    for (k=0; k < n; k++) {
                        dprod *= P1[(s*n+k)*(upperlimit+1) + (int)data[l*n+k] - 1];
                    }
                    dprod *= P2[t * SizeYL + (y-1) * L + s];
                    dylt += dprod;
                }
                P_YMN[(y-1)* m * n + l * n + t] = dylt / (double)L;
                //end of inner loop
                
            }
        }
    }
    
    //calculate the left-hand side of 1.4.3, which has 2 indexes (y,t) 
    for (t=0; t< n; t++) {
        double dsum = 0.0;
        for (y=1;y <= upperlimit+1; y++) {
            double dprod = 1.0;
            for (l=0; l < m; l++) {
                dprod *= P_YMN[(y-1)* m * n + l * n + t];
            }
            result[t * (upperlimit+1) + y -1] = dprod;
            dsum += dprod;
        }
        for (y=1;y <= upperlimit+1; y++) {
            result[t * (upperlimit+1) + y -1] /= dsum;
        }
    }
                 
    delete [] P1;
    delete [] P2;
    delete [] P_YMN;
                    
}

/*
 *draw samples from a list of lambdas
Input:
    lambds:     a 1 by m vector of lambdas
    m:          Number of lambdas
    n:          Number of samples per lambda
    upperlimt:  the upper truncated limit
    mtr:        the random number generator
 Output:        
    result: the random samples. Caller is responable for allocating memory, (samples are in rows)
 */
void SampleLambdas(double *lambdas, double *result, int m, int n, int upperlimit, MTRand & mtr) {
    int i,j;
    double *p = new double[m*(upperlimit+1)];
    UpperTruncatedPoissonPDF(lambdas, p, m, upperlimit+1);

    for (i=0; i< m; i++) {
        for (j =0; j < n; j++) {
            double d = mtr.rand();
            result[j*m+i] = Sample(p + i*(upperlimit+1), (upperlimit+1), d);
        }
    }
    delete [] p;               
}

double *ReadData(int n, int ncol, string infile, int verbose) {
	double *data = new double[n * ncol]; 
	int i,j;
	ifstream theFile(infile.c_str());
	if (theFile.fail()) {
		std::cout << "Failed to open the file " << infile.c_str() << endl;
		exit(1);
	}
	
	double dCurrent = 0;
	for (i = 0; i < n*ncol; i++) {
		if (!theFile.eof()) {
			theFile >> dCurrent;
			data[i] = dCurrent;
		} else {
			std::cout << "Not enough numbers to be read" << endl;
			exit(1);
		}
	}
	theFile.close();

	//debug output, to make sure the data is read in correctly
	if (verbose) {
		for (i = 0; i<n; i++){
			for (j = 0; j < ncol; j++) {
				fprintf(stdout,"%d ",(int)data[i*ncol+j]);
			}
			fprintf(stdout,"\n");
		}
	}	
	return data;
}

void CalculateRiskSummary(double *data, int n, double *riskmeasure, double *result, int in_upper) {
	int *rt = new int[n];
	int *at = new int[n];
	int *arg_y = new int[n];
	int i,j;

	//calculate rt
	int ncol = 5;
	int nsum = 0;
	for (i =0; i < n;i++) {
		//find the maximum value
		double dmax = -1;
		for (j=0; j < in_upper+1; j++) {
			if (riskmeasure[i *(in_upper+1) + j] > dmax) { dmax = riskmeasure[i *(in_upper+1) + j];}
		}

		//
		rt[i] = 0;
		arg_y[i] = -1;
		for (j=0; j < in_upper+1; j++) {
			if (fabs(riskmeasure[i *(in_upper+1) + j] - dmax) < 1e-20) {
				int y = (int)data[i*	ncol + 2];
				if (y == j+1) {
					rt[i]++;
				}
				arg_y[i] = j+1;
			}
		}
		if (rt[i] > 1) { //ties, then count as zero
			rt[i] = 0;
		}

		nsum += rt[i];
		result[i] = fabs(data[i*	ncol + 2] - arg_y[i]);
	}
	double rall = (double)nsum / (double)n;

	int nsum_at = 0;
	int nsum_at_rt = 0;
	for (i = 0; i < n; i++) {
		int counts = 0; 
		int countt = 0;
		at[i] = 0;
		for (j = 0; j < n; j++) {
			if (data[i * ncol + 3] == data[j * ncol + 3]) {
				counts++;
			}
			if (data[i * ncol + 4] == data[j * ncol + 4]) {
				countt++;
			}
	
		}
		if (counts == 1 && countt == 1) {
			at[i]++;
		}
		//cout << (i+1) << "\t" << at[i] << endl;
		nsum_at += at[i];
		nsum_at_rt += at[i] * rt[i];
		
	}
	double runq = 0.0;
	if (nsum_at > 0) {
		runq = (double)nsum_at_rt / (double) nsum_at;
	} 

	delete [] rt;
	delete [] at;
	delete [] arg_y;
	result[n] = rall;
	result[n+1] = runq;

}
void Run(string infile, int n, unsigned int seed, int niters, int burnin, int stride, string outfile,
		 string syndatafile, string riskfile, int m, bool verbose, int in_upper) {
	
	int i,j;
	int ncol = 3; 
	int betaAcc,gammaAcc,L,npara;

	// random number intiliization
	MTRand mt;
	if (seed > 0) {
		mt.seed(seed);
	} else {
		mt.seed((unsigned)time( NULL ));
	}

	//read data 
	double *data = ReadData(n, ncol, infile, verbose);
    
	//yobs
	double *yobs = new double[n];
	for (i = 0; i < n; i++) {
		yobs[i] = data[i*ncol+2]; // y is in column 2
	}

	// model fitting
	int B,G;
	double *new_data = new double[n*(ncol+2)];
	double *models = DoModel(data,new_data, n, niters, burnin, stride, betaAcc, gammaAcc, L,B,G,mt,in_upper);
	
	npara = B + G + 4; 

	if (verbose) {
		//write the model file	
		FILE *OUT; OUT = fopen(outfile.c_str(),"w");
		for ( i = 0; i <L;i++) {
			for (j = 0; j < npara; j++) {
				fprintf(OUT, "%10.8g ",models[i*npara+j]);
			}
			fprintf(OUT, "\n");
		}
		fclose(OUT);
	}

	//calculate Lambdas from saved mcmc interations and sample synthetic data from all Lambdas
	double * Lambdas = new double[L * n]; 
	double *syndata = new double[L*n];
	for (i =0; i < L; i++) {
		double *MODEL = models + i*npara;
		double *BETA = MODEL + 4;
		double *GAMMA = BETA + B;
		CalculateLambda(new_data, BETA, GAMMA, Lambdas + i * n, n);
		SampleLambdas(Lambdas + i * n, syndata + i * n , n, 1, in_upper, mt);
	}
	
	double *CI = new double[n*2];
	//caculate confidence interval
	int lindex = ceil(L * 0.05);
	int uindex = floor(L * 0.95);
	for (i = 0; i < n; i++) {
		vector<double> v;
		for (j =0; j < L; j++) {
			v.push_back(syndata[j*n+i]);
		}
		std::sort(v.begin(),v.end());
		CI[i*2] = v[lindex];
		CI[i*2+1] = v[uindex];
	}
	
	double *riskmeasure = new double[n*(in_upper+1)];
	double *risksummary = new double[m*(n+2)];
	FILE *RISKMEASURE; 
	if (verbose) {
		RISKMEASURE = fopen(riskfile.c_str(),"w");
	}

	int *synIndex = new int[m];
	for (int s= 0; s< m; s++) {
		int index = (int)(mt.rand() * (double)L);
		synIndex[s] = index;
		ImportanceSampling(Lambdas, syndata + index*n, yobs, riskmeasure, L, 1, n, in_upper);
		if (verbose) {
			//output risk measure matrix to a output file
			for ( i = 0; i <n;i++) {
				for (j = 0; j < (in_upper+1); j++) {
					fprintf(RISKMEASURE, "%10.8g ",riskmeasure[i*(in_upper+1)+j]);
				}
				fprintf(RISKMEASURE, "\n");
			}
			//fprintf(RISKMEASURE, "\n");
		}
		CalculateRiskSummary(new_data, n, riskmeasure, risksummary+s*(n+2), in_upper);
	}
	if (verbose) {
		fclose(RISKMEASURE);
	}
	
	
	//output synthetic data
	FILE *SDATA; SDATA = fopen(syndatafile.c_str(),"w"); //transpose
	fprintf(SDATA, "sc dc y ");
	
	for (j = 0; j < m; j++) {
		fprintf(SDATA, "sd%d ",j+1);
	}
	
	fprintf(SDATA, "cl cu ");
	for (j = 0; j < m; j++) {
		fprintf(SDATA, "dt%d Rall%d Runq%d ",j+1, j+1, j+1);
	}
	fprintf(SDATA, "\n");

	for (i = 0; i < n; i++) {
		//output original data
		for (j=0; j < 3; j++) {
			fprintf(SDATA, "%d ",(int)data[i*3+j]);
		}
		for (j = 0; j < m; j++) {
			fprintf(SDATA, "%d ",(int)syndata[synIndex[j]*n+i]);
		}
		//Confidence interval
		fprintf(SDATA, "%d ",(int)CI[i*2]);
		fprintf(SDATA, "%d ",(int)CI[i*2+1]);
		
		for (j = 0; j < m; j++) {
			fprintf(SDATA, "%d ",(int)risksummary[j*(n+2)+i]);
			fprintf(SDATA, "%10.8g ",risksummary[j*(n+2)+ n]);
			fprintf(SDATA, "%10.8g ",risksummary[j*(n+2)+ n+1]);
		}
		
		fprintf(SDATA, "\n");
	}
	fclose(SDATA);


	delete [] data;
	delete [] models;
	delete [] new_data;
	delete [] syndata;
	delete [] yobs;
	delete [] Lambdas;
	delete [] riskmeasure;
	delete [] risksummary;
	delete [] CI;
	delete [] synIndex;

}

// START OF MIGRATED FUNCTIONS
//Functions are migrated since it's troublesome to add dependency
//and it's easier to just copy paste when the object has only
//one function used

double gammarand(double a, double b, MTRand& mt)
{
    
    //Return -1 if a or b is not positive.
    if ((a <= 0) || (b <= 0)) return -1.0;
    
    //If a == 1, then gamma is exponential. (Devroye, page 405).
    if (a == 1) {
        return -b * log(1 - mt.randExc());
    }
    
    //Devroye, page 418 Johnk's generator
    if ((a < 1) && (a > 0)) {
        double c = 1.0 / a;
        double d = 1.0 / (1 - a);
        while (1) {
            double dX = pow(1-mt.randExc(), c);
            double dY = pow(1-mt.randExc(), d);
            if (dX + dY <= 1.0) {
                double dE = -log(1 - mt.randExc());
                return b * dE * dX / (dX + dY );
            }
        }
    }
    
    //Use a rejection method for a > 1
    //Devroye, page 410 Best's algorithm
    if (a > 1) {
        double bb = a - 1;
        double c = 3 * a - 0.75;
        
        while (1) {
            double dU = 1-mt.randExc();
            double dV = 1-mt.randExc();
            
            double dW = dU * (1 - dU);
            double dY = sqrt(c / dW) * (dU - 0.5);
            double dX = bb + dY;
            if (dX >=0) {
                double dZ = 64 * dW * dW * dW * dV * dV;
                if (dZ <= (1.0 - 2.0 * dY * dY / dX)) {
                    return b * dX;
                } else {
                    if (log(dZ)<= 2 * (bb * log(dX / bb) - dY)) {
                        return b * dX;
                    }
                }
            }
        }
    }
    return -1.0;
}

//START OF CONVERSION HELPERS FOR EXPORTING FUNCTIONS
//conversion helper get numeric vector from array and the reverse also

NumericVector vectorFromArray(double *arr, int length) {
    NumericVector v = NumericVector(length);
    for (int i = 0; i<length; i++) {
        v[i] = arr[i];
    }
    return v;
}

double *arrayFromVector(NumericVector v) {
    
    int n = v.size();
    double *arr = new double[n];
    for (int i = 0; i<n; i++) {
        arr[i] = v[i];
    }
    return arr;
}

List listFromArray(double *arr, int ncol, int length) {
    
    List lst;
    
    int n = length/ncol;
    for (int j = 0; j<ncol; j++) {
      NumericVector v = NumericVector(n);
      for (int i = 0; i < n; i++) {
          v[i] = arr[i*ncol+j];
        }
      lst.push_back(v);
    }
    
    return lst;
}

List listForSynData(double *syndata, int*synIndex, int m, int n) {
    
    List lst;
    for (int j = 0; j < m; j++) {
        NumericVector v = NumericVector(n);
        for (int i = 0; i < n; i++) {
            v[i] = (int)syndata[synIndex[j]*n+i];
        }
        lst.push_back(v);
    }
    return lst;
    
}


//START OF EXPORTING FUNCTIONS
//EXPORTING FUNCTION links to Rcpp to provide R interface for all
//the C++ functions here


// [[Rcpp::export]]
List getSynData(NumericVector vdata, int n, int seed, int niters, int burnin, int stride, int m, bool verbose, int in_upper) {
    
    int i,j;
    int ncol = 3;
    int betaAcc,gammaAcc,L,npara;
    
    // random number intiliization
    MTRand mt;
    if (seed > 0) {
        mt.seed(seed);
    } else {
        mt.seed((unsigned)time( NULL ));
    }
    
    double *data = arrayFromVector(vdata);
    
    //yobs
    double *yobs = new double[n];
    for (i = 0; i < n; i++) {
        yobs[i] = data[i*ncol+2]; // y is in column 2
    }
    
    // model fitting
    int B,G;
    double *new_data = new double[n*(ncol+2)];
    double *models = DoModel(data,new_data, n, niters, burnin, stride, betaAcc, gammaAcc, L,B,G,mt,in_upper);
    
    npara = B + G + 4;
    
    //calculate Lambdas from saved mcmc interations and sample synthetic data from all Lambdas
    double * Lambdas = new double[L * n];
    double *syndata = new double[L*n];
    for (i =0; i < L; i++) {
        double *MODEL = models + i*npara;
        double *BETA = MODEL + 4;
        double *GAMMA = BETA + B;
        CalculateLambda(new_data, BETA, GAMMA, Lambdas + i * n, n);
        SampleLambdas(Lambdas + i * n, syndata + i * n , n, 1, in_upper, mt);
    }
    
    double *CI = new double[n*2];
    //caculate confidence interval
    int lindex = ceil(L * 0.05);
    int uindex = floor(L * 0.95);
    for (i = 0; i < n; i++) {
        vector<double> v;
        for (j =0; j < L; j++) {
            v.push_back(syndata[j*n+i]);
        }
        std::sort(v.begin(),v.end());
        CI[i*2] = v[lindex];
        CI[i*2+1] = v[uindex];
    }
    
    double *riskmeasure = new double[n*(in_upper+1)];
    double *risksummary = new double[m*(n+2)];
    
    
    int *synIndex = new int[m];
    for (int s= 0; s< m; s++) {
        int index = (int)(mt.rand() * (double)L);
        synIndex[s] = index;
        ImportanceSampling(Lambdas, syndata + index*n, yobs, riskmeasure, L, 1, n, in_upper);
        CalculateRiskSummary(new_data, n, riskmeasure, risksummary+s*(n+2), in_upper);
    }
    
    List lst1 = listForSynData(syndata,synIndex,m,n);
    List lst2 = listFromArray(CI,2,2*n);

    lst1.push_back(lst2[0]);  lst1.push_back(lst2[1]);

    return lst1;
    
}
