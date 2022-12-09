#include<vector>
#include<Rcpp.h>
#include<math.h>
#include<iostream>

using namespace std;
using namespace Rcpp;

vector<double> groundTruth;
vector<vector<double> > inputValues;
vector<double> weight;


long epochs = 10;
double lr = 0.001;
double e = 2.71828;


void updateWeight(double predicted, double expected, vector<double> inputs);
void getAcc();
double sigmoid(double z);


NumericVector logistic_regression(NumericMatrix explain, NumericVector response);

//' @title A logistic regression method using Rcpp
//' @description A logistic regression method using Rcpp
//' @param explain the design matrix 
//' @param response the response vector y
//' @return the estimated coefficients of variables
//' @export
// [[Rcpp::export]]
NumericVector logistic_regression(NumericMatrix explain, NumericVector response) {
  int nr = explain.nrow();
  int nc = explain.ncol();
  for (int i = 0; i < nc; i++) {
    weight.push_back(0.2);
  }
  long i, j;
  
  // fill every sample attribution of dataset
  for (int k = 0; k < nr; k++) {
    vector<double> inputRow;
    for (int kk = 0; kk < nc; kk++) {
      inputRow.push_back(explain(k, kk));
    }
    inputValues.push_back(inputRow);
  }
  
  // fill ground truth
  int n_true = response.size();
  for (int k = 0; k < n_true; k++) {
    groundTruth.push_back(response[k]);
  }
  
  while (epochs--) {
    // cout << "********************----" << epochs << endl;
    for (i = 0; i < inputValues.size(); i++) {
      double yhat, z = 0;
      for (j = 0; j < inputValues[0].size(); j++) {
        z += weight[j] * inputValues[i][j];
      }
      yhat = sigmoid(z);
      
      updateWeight(yhat, groundTruth[i], inputValues[i]);
    }
  }
  
  
  getAcc();
  NumericVector out(nc);
  for (int k = 0; k < nc; ++k) {
    out[k] = weight[k];
  }
  return out;
}

void updateWeight(double yhat, double groundTruth, vector<double> inputValue) {
  for (int i = 0; i < inputValue.size(); i++) {
    double gradientDescent;
    gradientDescent = (yhat - groundTruth) * inputValue[i];
    weight[i] = weight[i] - (lr * gradientDescent);
  }
}

double sigmoid(double z) {
  return 1/(1 + pow(e, (-1 * z)));
}

void getAcc() {
  
  long correct = 0, totalCases = inputValues.size();
  
  for (int i = 0; i < totalCases; i++) {
    double yhat, z = 0;
    for (int j = 0; j < inputValues[0].size(); j++) {
      z += weight[j] * inputValues[i][j];
    }
    
    
    yhat = (sigmoid(z) < 0.5) ? 0 : 1;
    
    correct += (yhat == groundTruth[i]);
    
  }
  //cout << "Accuracy is: " << (correct * 100) / (totalCases) << "%" << endl;
  
}