#include "memwrapper.h"
#include "nn.h"
#include "numeric.h" /* Using:  if(FLOAT_EQ(NumOne, NumTwo));*/
#include "metricspace.h"
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>


typedef struct{
  //weights matrix
  matrix *wi;
  matrix *wo;
  //last change in weights for momentum 
  matrix *ci;
  matrix *co;
  
  // activations for nodes
  dvector *ai;
  dvector *ah;
  dvector *ao;
  size_t ni, no, nh;
} NNModel;

// sigmoid function
double sigmoid(double x){
  return tanh(x);
}

// derivative of sigmoid function, in terms of the output (i.e. y)
double dsigmoid(double y){
  return 1.0 - pow(y, 2);
}

/*
 * Simple Backpropagation Neural Network
 * input:
 * ni = number of input nodes (generally are the number of variables of the dataset)
 * nh = number of hidden nodes (variable whic could be optimized manually)
 * no = number of output nodes (generally are the number of predictors)
 * lr = learning rate
 * mf = momentum factor
 * 
 * regression = 1 if is a regression problem and 0 if is a classification problem
 */
void NN(matrix *mx, matrix *my, size_t no, double lr, double mf){
  
}

void update(dvector *xrow, NNModel *nn, int regression){
  if(xrow->size != ni-1){
    fprintf(stderr, "wrong number of variables\n");
  }
  else{
    size_t i, j, k;
    // input activations
    for(i = 0; i < nn->ni-1; i++)
      nn->ai->data[i] = xrow->data[i];

    // hidden activations
    for(j = 0; j < nn->nh - 1; j++){
      double total = 0.f;
      for(i = 0; i < nn->ni; i++)
        total += nn->ai->data[i] * nn->wi->data[i][j];
      nn->ah[j] = sigmoid(total);
    }
    
    // output activations
    for(k = 0; k < nn->no; k++){
      double total = 0.f;
      for(j = 0 j < nn->nh; j++)
        total += nn->ah->data[j] * nn->wo->data[j][k];
      nn->ao->data[k] = total;
      
      if(regression == 0)
        nn->ao[k] = sigmoid(total);
    }
//   return self.ao[:]
  }
}


double BackPropagate(dvector *yrow, NNModel *nn, double lr, double mf, int regression)
{
  if(yrow->size != nn->no){
    fprintf(seterr, "wrong number of y\n");
    return 99999999;
  }
  else{
    size_t i, j, k;
    // calculate error terms for output
    dvector *odelta;
    NewDVector(&odelta, nn->no);
    DVectorSet(odelta, 0.f);
 
    for(k = 0; k < nn->no; k++){
      odelta->data[k] = yrow->data[k] - nn->ao->data[k];
      if(regression == 0)
        odelta->data[k] = dsigmoid(nn->ao->data[k]) * odelta->data[k];
    }

    // calculate error terms for hidden
    dvector *hdelta;
    NewDVector(&hdelta, nn->nh);
    DVectorSet(hdelta, 0.f);
    for(j = 0; j < nn->nh; j++){
      double error = 0.f;
      for(k = 0; k < nn->no; k++)
        error += odelta->data[k] * nn->wo->data[j][k];
      hdelta->data[j] = dsigmoid(nn->ah->data[j]) * error;
    }
    
    // update output weights
    for(j = 0; j < nn->nh; j++){
      for(k = 0; k < nn->no; k++){
        double change = odelta->data[k] * nn->ah->data[j];
        nn->wo->data[j][k] = nn->wo->data[j][k] + lr * change + mf*nn->co->data[j][k];
        nn->co->data[j][k] = change;
      }
    }

    // update input weights
    for(i = 0; i < nn->ni; i++){
      for(j = 0; j < nn->nh; j++){
        double change = hdelta->data[j]*nn->ai->data[i];
        nn->wi->data[i][j] = nn->wi->data[i][j] + lr * change + mf * nn->ci->data[i][j];
        nn->ci->data[i][j] = change;
      }
    }
    
    // calculate error
    double error = 0.f;
    for(k = 0; k < yrow->size; k++){
      error += 0.5*pow((yrow->data[k]-nn->ao->data[k]),2);
    }
    return error;
  }
}
