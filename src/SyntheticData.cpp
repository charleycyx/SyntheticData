//
//  SyntheticData.cpp
//  
//
//  Created by Charley Chen on 7/31/15.
//
//

#include <Rcpp.h>
#include "model.h"

using namespace Rcpp;
using namespace std;

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

List listFromRisk(double *risksummary, int n, int m) {
    
    List lst;
    for (int j = 0; j < m; j++) {
        NumericVector v1 = NumericVector(n);
        NumericVector v2 = NumericVector(n);
        NumericVector v3 = NumericVector(n);
        for (int i = 0; i < n; i++) {
            v1[i] = (int)risksummary[j*(n+2)+i];
            v2[i] = risksummary[j*(n+2)+ n];
            v3[i] = risksummary[j*(n+2)+ n+1];
        }
        lst.push_back(v1);  lst.push_back(v2);  lst.push_back(v3);
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
    List lst3 = listFromRisk(risksummary, n, m);
    
    //push list 2 and 3 into list 1 to form an entire list
    lst1.push_back(lst2[0]);  lst1.push_back(lst2[1]);
    for (i=0; i<lst3.size(); i++) {
        lst1.push_back(lst3[i]);
    }
    
    return lst1;
    
}