
/*
 * Description: Calculate the PLS Bootstrap Random group cross validation.
 */
void PLSRandomGroupsCV(matrix *mx, matrix *my, size_t xautoscaling, /*Inputs*/
                       size_t yautoscaling, size_t nlv, size_t group, size_t iterations, /*Inputs*/
                      matrix **q2y, matrix **sdep, matrix **bias, matrix **predicted_y, /*Ouputs*/
                      matrix **pred_residuals, size_t nthreads, ssignal *s); /*Ouputs*/

/*
 * Description: Calculate the PLS leave one out cross validation.
 */
void PLSLOOCV(matrix *mx, matrix *my, size_t xautoscaling, size_t yautoscaling, size_t nlv,/*Inputs*/
                        matrix **q2y, matrix **sdep, matrix **bias, /*Ouputs*/
                        matrix **predicted_y, matrix **pred_residuals, /*Ouputs*/
                        size_t ntreads, ssignal *s); /*Inputs*/


/*
 * Description: Validate the model stability by reducing and fixing the sample size
*/
void PLSStaticSampleValidator(matrix *mx, matrix *my, uivector *obj_class,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t sample_size, size_t niters,
                        size_t rgcv_group, size_t rgcv_iterations, size_t nthreads,
                        matrix **q2_distr, matrix **sdep_distr, uivector **bestid, ssignal *s);

/*
 * Description: Validate the model stability using a dinamic incremental object sampling
 */
void PLSDynamicSampleValidator(matrix *mx, matrix *my,
                        size_t xautoscaling, size_t yautoscaling,
                        size_t nlv, size_t niters,
                        uivector *obj_class, size_t deltaobj, size_t maxobj,
                        size_t rgcv_group, size_t rgcv_iterations, size_t nthreads,
                        matrix **q2_surface, matrix **sdep_surface, uivector **bestid, ssignal *s);
