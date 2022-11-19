#include <stdio.h>
#include <scientific.h>

int main(void)
{
    matrix *x, *y; // Define the feature matrix x and the target to predict y
    dvector *betas; // Define the beta coefficients
    PLSMODEL *m;

    // Allocate the matrix
    NewMatrix(&x, 14, 6);
    NewMatrix(&y, 14, 1);

    // Fill the matrix with values
    // This is a manual filling.
    // Of course we can read a csv file and fill it automatically

    x->data[0][0] = 4.0000;  x->data[0][1] = 4.0000;  x->data[0][2] = 1.0000;  x->data[0][3] = 84.1400;  x->data[0][4] = 1.0500;  x->data[0][5] = 235.1500;
    x->data[1][0] = 5.0000;  x->data[1][1] = 5.0000;  x->data[1][2] = 1.0000;  x->data[1][3] = 79.1000;  x->data[1][4] = 0.9780;  x->data[1][5] = 231;
    x->data[2][0] = 4.0000;  x->data[2][1] = 5.0000;  x->data[2][2] = 1.0000;  x->data[2][3] = 67.0900;  x->data[2][4] = 0.9700;  x->data[2][5] = 249.0000;
    x->data[3][0] = 4.0000;  x->data[3][1] = 4.0000;  x->data[3][2] = 1.0000;  x->data[3][3] = 68.0700;  x->data[3][4] = 0.9360;  x->data[3][5] = 187.3500;
    x->data[4][0] = 3.0000;  x->data[4][1] = 4.0000;  x->data[4][2] = 2.0000;  x->data[4][3] = 68.0800;  x->data[4][4] = 1.0300;  x->data[4][5] = 363.0000;
    x->data[5][0] = 9.0000;  x->data[5][1] = 7.0000;  x->data[5][2] = 1.0000;  x->data[5][3] = 129.1600;  x->data[5][4] = 1.0900;  x->data[5][5] = 258.0000;
    x->data[6][0] = 10.0000;  x->data[6][1] = 8.0000;  x->data[6][2] = 0.0000;  x->data[6][3] = 128.1600;  x->data[6][4] = 1.1500;  x->data[6][5] = 352.0000;
    x->data[7][0] = 6.0000;  x->data[7][1] = 6.0000;  x->data[7][2] = 0.0000;  x->data[7][3] = 78.1118;  x->data[7][4] = 0.8765;  x->data[7][5] = 278.6400;
    x->data[8][0] = 16.0000;  x->data[8][1] = 10.0000;  x->data[8][2] = 0.0000;  x->data[8][3] = 202.2550;  x->data[8][4] = 1.2710;  x->data[8][5] = 429.1500;
    x->data[9][0] = 6.0000;  x->data[9][1] = 12.0000;  x->data[9][2] = 0.0000;  x->data[9][3] = 84.1600;  x->data[9][4] = 0.7800;  x->data[9][5] = 279.0000;
    x->data[10][0] = 4.0000;  x->data[10][1] = 8.0000;  x->data[10][2] = 1.0000;  x->data[10][3] = 72.1100;  x->data[10][4] = 0.8900;  x->data[10][5] = 164.5000;
    x->data[11][0] = 4.0000;  x->data[11][1] = 9.0000;  x->data[11][2] = 1.0000;  x->data[11][3] = 71.1100;  x->data[11][4] = 0.8660;  x->data[11][5] = 210.0000;
    x->data[12][0] = 5.0000;  x->data[12][1] = 11.0000;  x->data[12][2] = 1.0000;  x->data[12][3] = 85.1500;  x->data[12][4] = 0.8620;  x->data[12][5] = 266.0000;
    x->data[13][0] = 5.0000;  x->data[13][1] = 10.0000;  x->data[13][2] = 1.0000;  x->data[13][3] = 86.1300;  x->data[13][4] = 0.8800;  x->data[13][5] = 228.0000;

    y->data[0][0] = 357.1500;
    y->data[1][0] = 388.0000;
    y->data[2][0] = 403.0000;
    y->data[3][0] = 304.5500;
    y->data[4][0] = 529.0000;
    y->data[5][0] = 510.0000;
    y->data[6][0] = 491.0000;
    y->data[7][0] = 353.3000;
    y->data[8][0] = 666.6500;
    y->data[9][0] = 354.0000;
    y->data[10][0] = 339.0000;
    y->data[11][0] = 360.0000;
    y->data[12][0] = 379.0000;
    y->data[13][0] = 361.0000;

    // Allocate the PLS model
    NewPLSModel(&m);

    /* Calculate the partial least squares algorithm taking as input:
     * x: the feature matrix x
     * y: the target matrix y
     * nlv: the number of latent variable nlv
     * xautoscaling: the autoscaling type for the x matrix
     * yautoscaling: the autoscaling type for the y Matrix
     * model: the PLSMODEL previously allocated
     * ssignal: a scientific signal to stop the calculation if requested by the user
     *
     * more information in the pls.h header file
     * void PLS(matrix *mx, matrix *my, size_t nlv, size_t xautoscaling, size_t yautoscaling, PLSMODEL *model, ssignal *s);
     */
    PLS(x, y, 3, 1, 0, m, NULL);

    PrintPLSModel(m); // Print to video the PLS model

    /*Validate the model using the internal validation method*/
    MODELINPUT minpt; // Define the model input for the validation method
    minpt.mx = &x;
    minpt.my = &y;
    minpt.nlv = 3;
    minpt.xautoscaling = 1;
    minpt.yautoscaling = 0;

    // Use the boot strap random group cross validation.
    BootstrapRandomGroupsCV(&minpt, 3, 100, _PLS_, m->predicted_y, m->pred_residuals, 4, NULL, 0);
    // We can also compute leave one out in case...
    // LeaveOneOut(&minpt, _PLS_, m->predicted_y, m->pred_residuals, 4, NULL, 0);

    // Calculate the model validation statistics
    PLSRegressionStatistics(y, m->predicted_y, m->q2y, m->sdep, m->bias);
    //Print to video the results of the validation and the predicted values

    puts("Q2 Cross Validation");
    PrintMatrix(m->q2y);
    puts("SDEP Cross Validation");
    PrintMatrix(m->sdep);
    puts("BIAS Cross Validation");
    PrintMatrix(m->bias);;

    // Calculate the beta coefficients to see the importance of each feature
    puts("Beta coefficients");
    initDVector(&betas);
    PLSBetasCoeff(m, GetLVCCutoff(m->q2y), betas); // GetLVCCutoff select the best Q2 value from all the possibilities
    PrintDVector(betas);

    puts("PREDICTED VALUES");
    PrintMatrix(m->predicted_y);

    puts("PREDICTED RESIDUALS");
    PrintMatrix(m->pred_residuals);

    puts("REAL Y");
    PrintMatrix(y);

    // Free the memory spaces
    DelDVector(&betas);
    DelPLSModel(&m);
    DelMatrix(&x);
    DelMatrix(&y);
}
