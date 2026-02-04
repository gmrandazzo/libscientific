/* Implements variable selection algorithms.
 * Copyright (C) 2018-2026 designed, written and maintained by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */


/* Variable Selection using the Particle Swarm Optimization algorithm
 *
 * Input:
 *   mx and my matrix for make model; optionally px and py for testset validation.
 *   randomness = a number betwee 0 and 1 that represent the randomness to add on the first model generation.
 *   nrandomnessiter = for compute random models and evaluate also the model consistency...
 *
 * Output:
 * - varselected vector: vector of varialble selected
 * - map = a matrix of 4 column with the best r2, q2, f-test and number of variables.
 * - vardistribution = a vector of variables distribution accumulated in all model. Each row is a variable and the value of this row is the number of time that this variable was used.
 */
void PSLPSOVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                       size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                       size_t population_size, double randomness,
                       uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s);

/*
 * Variable selection using the Genetic Algorithm
 *
 * Input
 *   mx and my matrix for make model; optionally px and py for testset validation.
 *
 * population_size = population of models to calculate; is a number > 2
 * fraction_of_population = fraction of population; is a number between 0 and 1
 * mutation_rate = rate of mutation; is a number between 0 and 1
 * crossovertype = 0 for single point crossover, 1 for two point crossover and 2 for uniform crossover
 * nswapping = used by the uniform crossover (crossovertype = 2) and is a number between 0 and 1. It means the % of variable to swapping
 *
 * Output
 *   varselected vector: vector of varialble selected
 *   rq2 = coefficient of validation q2  or r2
 *  f_distribution = coefficient of distribution
 *  map = a matrix of 4 column with the best r2, q2, f-test and number of variables.
 */
void PLSGAVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                       size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                       size_t population_size, double fraction_of_population, double mutation_rate, size_t crossovertype, double nswapping, double populationconvergence,
                       uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s);

/* Variable selection using the Spearmann correlation coefficient
 * Input
 * mx and my matrix for make the model; optionally px and py for testset validation
 *
 * nlv = number of principal component
 * validation_type = validation type: 0 = LOO; 1 RGCV
 * ngroup = number of group for RGCV
 * niter = number of iterations for RGCV
 * nrandomnessiter = number of random iteration for Model consistency
 *
 * threshold = the step size to check if variable is in or out
 *
 * Output:
 *  varselected vector: vector of varialble selected
  *  map = a matrix of 4 column with the best r2, q2, treshold and number of variables.
 */
void PLSSpearmannVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                                   size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type,
                                   size_t ngroup, size_t niter,
                                   double threshold,
                                   uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s);
