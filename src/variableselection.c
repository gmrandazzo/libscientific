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


/*
 * Variable Selection using the Particle Swarm Optimization algorithm
 *
 * Input:
 *   mx and my matrix for make model; optionally px and py for testset validation.
 *
 * Output:
 * - varselected vector: vector of varialble selected
 * - validation vector:  r2 if you have a test set or q2 if the variable selection is made into de model
 *
 *
 * Each particle is a model. Each model have is velocity and his position.
 *
 * Algoritm:
 *
 * 1) Initialize the particle
 *
 * Population = The model with all varible setted
 * P_g_best = Q^2 of the model as is with all variable or R^2 if you have a prediction!
 *
 * 2) Initialize the models particle.
 * for(i = 1 to i = Population Size):
 *   P_velocity = randomvelocity for the model
 *   P_position = randomposition. This means that the variable could be setted on or off in the model.
 *
 *   P_best = Cost of model with the actual position. So calculation of R^2 or Q^2.
 *
 *   if the cost of actual model P_best is <= to the global cost P_g_best:
 *     P_g_best = P_best
 * end
 *
 * While(R^2 or Q^2 is not good enough...):
 *   For Model P in population:
 *     P_new_velocity of the model = UpdateVelocity(P_old_velocity, P_g_best, P_best)
 *     P_new_position of the model = UpdatePosition(P_old_position, P_velocity)
 *
 *     P_actual_model = R^2 or Q^2 of the actual model setted..
 *
 *     if P_actual_model < = P_best:
 *       P_best = P_actual_model
 *       if P_best <= P_g_best:
 *         P_g_best = P_best
 *
 *
 *
 * N.B.: 0 < rand() < 1
 *
 * UpdateVelocity(P_old_velocity, P_g_best, P_best):
 *   P_new_velocity = P_old_velocity + (c1 *rand() * (P_best - P_old_position)) +
 *     (c2 * rand() * (P_g_best + P_old_position))
 *
 * UpdatePosition(P_old_position, P_new_velocity):
 *   P_new_position = P_old_position + P_new_velocity
 *
 * P_new_position is the position to set on or off so is a number between 1 and the number of variables
 *
 */

void PLSCostPopulation(matrix *mx, matrix *my,
                      matrix *px, matrix *py,
                      uivector *popvector,
                      size_t xautoscaling, size_t yautoscaling, size_t nlv,
                      int validation_type, size_t ngroup, size_t niter,
                      double *r2, double *q2, double *bias_, size_t nthreads, ssignal *s)
{
  size_t i, j, k, subsize, cutoff;
  double r2_, q2_, _bias_;
  matrix *submx, *subpx, *r2y, *sdec, *q2y, *sdep, *bias;
  PLSMODEL *m;


  subsize = 0;
  if(popvector != NULL){
    for(i = 0; i < popvector->size; i++){
      if(getUIVectorValue(popvector, i) == 1){
        subsize++;
      }
      else{
        continue;
      }
    }
    NewMatrix(&submx, mx->row, subsize);
    for(i = 0; i < mx->row; i++){
      k = 0;
      for(j = 0; j < mx->col; j++){
        if(getUIVectorValue(popvector, j) == 1){
          setMatrixValue(submx, i, k, getMatrixValue(mx, i, j));
          k++;
        }
        else{
          continue;
        }
      }
    }
  }
  else{
    initMatrix(&submx);
    MatrixCopy(mx, &submx);
    subsize = mx->col;
  }

  if(px != NULL && py != NULL && validation_type == 0){ /* cost = R^2*/
    NewPLSModel(&m);
    PLS(submx, my, nlv, xautoscaling, yautoscaling, m, s);

    NewMatrix(&subpx, px->row, subsize);

    k = 0;
    for(i = 0; i < px->row; i++){
      for(j = 0; j < px->col; j++){
        if(getUIVectorValue(popvector, j) == 1){
          setMatrixValue(subpx, i, k, getMatrixValue(px, i, j));
          k++;
        }
        else{
          continue;
        }
      }
    }

    initMatrix(&r2y);
    initMatrix(&sdec);

    PLSRegressionStatistics(subpx, py, m, nlv, &r2y, &sdec, NULL);

    cutoff = GetLVCCutoff(r2y);
//     MatrixGetMaxValue(r2y, &cutoff, NULL);
    r2_ = 0.f;
    for(j = 0; j < r2y->col; j++){
      r2_ += getMatrixValue(r2y, cutoff, j);
    }

    if(r2_ != 0){
      r2_ /= (double)r2y->col;
    }
    else{
      r2_ = 0.f;
    }

    if(r2!= NULL)
      (*r2) = (*q2) = r2_;
    else
      (*q2) = r2_;
    if(bias_ != NULL)
      (*bias_) = 0.f;

    DelMatrix(&sdec);
    DelMatrix(&r2y);
    DelMatrix(&subpx);
    DelPLSModel(&m);
  }
  else if(validation_type == 1){ /*LEAVE ONE OUT*/
    initMatrix(&q2y);
    initMatrix(&sdep);
    initMatrix(&bias);
    PLSLOOCV(submx, my, xautoscaling, yautoscaling, nlv,
                      &q2y,
                      &sdep,
                      &bias,
                      NULL, NULL, nthreads, s);

    NewPLSModel(&m);
    PLS(submx, my, nlv, xautoscaling, yautoscaling, m, s);

    PLSRegressionStatistics(submx, my,  m, nlv, &m->r2y_model, &m->sdec, NULL);

    cutoff = GetLVCCutoff(q2y);
//     MatrixGetMaxValue(q2y, &cutoff, NULL);

    q2_ = 0.f;
    for(j = 0; j < q2y->col; j++){
      q2_ += getMatrixValue(q2y, cutoff, j);
    }

    _bias_ = 0.f;
    for(j = 0; j < bias->col; j++){
      _bias_ += getMatrixValue(bias, cutoff, j);
    }

    r2_ = 0.f;
    for(j = 0; j < m->r2y_model->col; j++){
      r2_ += getMatrixValue(m->r2y_model, cutoff, j);
    }

    if(q2_ != 0){
      q2_ /= (double)q2y->col;
      r2_ /= (double)m->r2y_model->col;
      _bias_ /= (double)bias->col;
    }
    else{
      q2_ = 0.f;
      r2_ = 0.f;
      _bias_ = 99;
    }

    (*q2) = q2_;
    if(r2 != NULL)
      (*r2) = r2_;
    if(bias_ != NULL)
      (*bias_) = _bias_;

    DelPLSModel(&m);
    DelMatrix(&q2y);
    DelMatrix(&sdep);
    DelMatrix(&bias);
  }
  else{ /*CROSS VALIDATION*/
    initMatrix(&q2y);
    initMatrix(&sdep);
    initMatrix(&bias);
    PLSRandomGroupsCV(submx, my, xautoscaling, yautoscaling, nlv, ngroup, niter,
                  &q2y,
                  &sdep,
                  &bias,
                  NULL, NULL, nthreads, s);

    NewPLSModel(&m);
    PLS(submx, my, nlv, xautoscaling, yautoscaling, m, s);
    PLSRegressionStatistics(submx, my,  m, nlv, &m->r2y_model, &m->sdec, NULL);

    cutoff = GetLVCCutoff(q2y);
//     MatrixGetMaxValue(q2y, &cutoff, NULL);

    q2_ = 0.f;
    for(j = 0; j < q2y->col; j++){
      q2_ += getMatrixValue(q2y, cutoff, j);
    }

    _bias_ = 0.f;
    for(j = 0; j < bias->col; j++){
      _bias_ += getMatrixValue(bias, cutoff, j);
    }

    r2_ = 0.f;
    for(j = 0; j < m->r2y_model->col; j++){
      r2_ += getMatrixValue(m->r2y_model, cutoff, j);
    }

    if(q2_ != 0){
      q2_ /= (double)q2y->col;
      r2_ /= (double)m->r2y_model->col;
      _bias_ /= (double)bias->col;
    }
    else{
      q2_ = 0.f;
      r2_ = 0.f;
      _bias_ = 99;
    }

    (*q2) = q2_;
    if(r2 != NULL)
      (*r2) = r2_;

    if(bias_ != NULL)
      (*bias_) = _bias_;

    DelPLSModel(&m);
    DelMatrix(&q2y);
    DelMatrix(&sdep);
    DelMatrix(&bias);
  }

  DelMatrix(&submx);
}


void PSLPSOVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                       size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                       size_t population_size, double randomness,
                       uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s)
{
  size_t i, j, iter, new_position, population_size_, tmp_var_id, nvarON;
  double model_r2, cost_g_best, cost_best, tmp_model_r2, tmp_cost_best, tmp_g_best, new_v, new_pos, rand1, rand2;
  matrix *models, *position, *velocity;
  uivector *popvector, *best_model;
  dvector *b_position, *g_position, *rowmap;

  /*1) Initialize the population*/

  population_size_ = population_size + 1; /* + 1 is the model with all variables */

  NewMatrix(&models, population_size_, mx->col);
  NewMatrix(&position, population_size_, mx->col); /*position chagned stored if are -1 then the position is not stored.... else is  a number >= 0 and < mx->col*/
  NewMatrix(&velocity, population_size_, mx->col);
  NewDVector(&b_position, mx->col);
  NewDVector(&g_position, mx->col);
  NewUIVector(&popvector, mx->col);
  NewUIVector(&best_model, mx->col);
  NewDVector(&rowmap, 4);

  UIVectorSet(popvector, 1); /*get all variables*/

  PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &model_r2, &cost_g_best, NULL, nthreads, s);

  setDVectorValue(rowmap, 0, model_r2);
  setDVectorValue(rowmap, 1, cost_g_best);
  setDVectorValue(rowmap, 2, (cost_g_best  * (population_size - mx->col -1)) /  (mx->col*(1-cost_g_best )));
  setDVectorValue(rowmap, 3, mx->col);

  MatrixAppendRow(map, rowmap);

  for(i = 0; i < mx->col; i++){
    UIVectorAppend(vardistribution, 1);
  }

  /*the best and global position are the same in the first model that have all the variable setted ON*/
  for(i = 0; i < mx->col; i++){
    /*All variable are global and best...*/
    setDVectorValue(g_position, i, i);
    setDVectorValue(b_position, i, i);
    setMatrixValue(velocity, 0, i, 0.f);
  }

  /*
  printf("model 0 with all variable: %f\n", cost_g_best);
  */

  /*2) Initialize the models particle.*/
  MatrixSet(models, 1.0); // all variables are setted ON!
  MatrixSet(position, -1.0); // all variables on model are unvariate...
  MatrixSet(velocity, 0); // all variables have no speed...

  srand(mx->row + mx->col + my->col + population_size_);

  for(i = 1; i < population_size_; i++){

    for(j = 0; j < (size_t)floor(mx->col * randomness); j++){ /* % of object to randomize in this case is the 20 % of variable*/
      setMatrixValue(velocity, i, j, (double)rand()/RAND_MAX); /*velocity is a number between 0 an 1*/

      new_position = rand() % mx->col;

      setMatrixValue(position, i, j, new_position);

      if((size_t)getMatrixValue(models, i, new_position) == 1){ /* set to OFF*/
        setMatrixValue(models, i, new_position, 0.0);
      }
      else{ /* Set to ON*/
        setMatrixValue(models, i, new_position, 1.0);
      }
    }

    /*COMPUTE THE COST OF ACTUAL MODEL RANDOM MODEL*/
    for(j = 0; j < models->col; j++){
      setUIVectorValue(popvector, j, (size_t)getMatrixValue(models, i, j));
    }

    nvarON = 0;
    for(j = 0; j < models->col; j++){
      if(getUIVectorValue(popvector, j) == 1){
        nvarON++;
      }
    }

    if(nvarON > 0){
      PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &tmp_model_r2, &cost_best, NULL, nthreads, s);
    }
    else{
      tmp_model_r2 = cost_best = 0;
    }

    for(j = 0; j < models->col; j++){
      setDVectorValue(b_position, j, getMatrixValue(position, i, j));
    }

    if(cost_best > cost_g_best){
      cost_g_best = cost_best;
      model_r2 = tmp_model_r2;

      for(j = 0; j < models->col; j++){
        setDVectorValue(g_position, j, getMatrixValue(position, i, j));
        setUIVectorValue(best_model, j, getUIVectorValue(popvector, j));
      }

      setDVectorValue(rowmap, 0, model_r2);
      setDVectorValue(rowmap, 1, cost_g_best);
      setDVectorValue(rowmap, 2, (cost_g_best  * (population_size - mx->col -1)) /  (mx->col*(1-cost_g_best )));
      setDVectorValue(rowmap, 3, nvarON);

      MatrixAppendRow(map, rowmap);
      DVectorSet(rowmap, 0.f);
    }

    /*
    printf("Modello %u Generato\n", (size_t)i);
    PrintUIVector(popvector);
    */
  }


  /*
  puts("modelli");
  PrintMatrix(models);

  puts("Best Position");
  PrintDVector(b_position);
  puts("Great Best Position");
  PrintDVector(g_position);

  puts("Second step!");
  */


/* While(R^2 or Q^2 is not good enough...):
 *   For Model P in population:
 *     P_new_velocity of the model = UpdateVelocity(P_old_velocity, P_g_best, P_best)
 *     P_new_position of the model = UpdatePosition(P_old_position, P_velocity)
 *
 *     P_actual_model = R^2 or Q^2 of the actual model setted..
 *
 *     if P_actual_model < = P_best:
 *       P_best = P_actual_model
 *       if P_best <= P_g_best:
 *         P_g_best = P_best
 */

  iter = 0;
  do{
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      /*
      puts("Velocity");
      PrintMatrix(velocity);
      puts("Positions");
      PrintMatrix(position);
      while (!getchar()=='\n');
      */

      tmp_g_best = cost_g_best;

      for(i = 0; i < models->row; i++){
        /*
        * update velocity and position and set variable on/off in model
        */
        for(j = 0; j < velocity->col; j++){
            rand1 = (double)rand()/(double)RAND_MAX;
            rand2 = (double)rand()/(double)RAND_MAX;
            new_v = getMatrixValue(velocity, i, j);
            new_v += (1 * rand1 * ((double)getDVectorValue(b_position, j) - getMatrixValue(position, i, j))); /*1 is a constant weight and could be variable */
            new_v += (1 * rand2 * ((double)getDVectorValue(g_position, j) - getMatrixValue(position, i, j))); /*1 is a constant weight and could be variable */
            setMatrixValue(velocity, i, j, new_v);

            /*
            * p(t+1) = p(t) + v(t); this value must be positive...
            */
            new_pos = fabs(floor(getMatrixValue(position, i, j) + new_v));
            setMatrixValue(position, i, j, new_pos);
            tmp_var_id = (size_t)getMatrixValue(position, i, j);

            if(new_pos > models->col-1 || new_pos < 0){
              /*
              printf("Error on setting the variable. Out of range: %f\n", getMatrixValue(position, i, j));
              */
              continue;
            }
            else{
              if((size_t)getMatrixValue(models, i, tmp_var_id) == 1){ /* set to OFF*/
                setMatrixValue(models, i, tmp_var_id, 0.0);
              }
              else{ /* Set to ON*/
                setMatrixValue(models, i, tmp_var_id, 1.0);
              }
            }
        }

        for(j = 0; j < models->col; j++){
          setUIVectorValue(popvector, j, (size_t)getMatrixValue(models, i, j));
          setUIVectorValue((*vardistribution), j, (getUIVectorValue((*vardistribution), j)+(size_t)getMatrixValue(models, i, j)));
        }


        nvarON = 0;
        for(j = 0; j < position->col; j++){
          if(getUIVectorValue(popvector, j) == 1){
            nvarON++;
          }
        }

        if(nvarON > 0){
          PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &tmp_model_r2, &tmp_cost_best, NULL, nthreads, s);
        }
        else{
          tmp_model_r2 = tmp_cost_best = 0;
        }

        if(tmp_cost_best > cost_best){
          cost_best = tmp_cost_best;

          for(j = 0; j < position->col; j++){
            setDVectorValue(b_position, j, fabs(floor(getMatrixValue(position, i, j))));
          }

          if(cost_best > cost_g_best){
            cost_g_best = cost_best;
            model_r2 = tmp_model_r2;

            for(j = 0; j < position->col; j++){
              setDVectorValue(g_position, j, fabs(floor(getMatrixValue(position, i, j))));
              setUIVectorValue(best_model, j, getUIVectorValue(popvector, j));
            }

            setDVectorValue(rowmap, 0, model_r2);
            setDVectorValue(rowmap, 1, cost_g_best);
            setDVectorValue(rowmap, 2, (cost_g_best  * (population_size - mx->col -1)) /  (mx->col*(1-cost_g_best )));
            setDVectorValue(rowmap, 3, nvarON);

            MatrixAppendRow(map, rowmap);

            DVectorSet(rowmap, 0.f);

            iter = 0;
          }
        }

        /*
        printf("iteration %u best_val %f global_best_value %f\n", (size_t) iter, cost_best, cost_g_best);
        */
        iter++;
      }
    }
  }while(FLOAT_EQ(cost_g_best, tmp_g_best, PLSCONVERGENCE) && iter < 100);

  /*
  PrintDVector(g_position);

  printf("Selected Model: %u\n", (size_t)mod_id);
  PrintMatrix(models);
  */

  for(i = 0; i < best_model->size; i++){
    UIVectorAppend(varselected, getUIVectorValue(best_model, i));
  }

  DelUIVector(&best_model);
  DelDVector(&rowmap);
  DelDVector(&g_position);
  DelDVector(&b_position);
  DelUIVector(&popvector);
  DelMatrix(&velocity);
  DelMatrix(&position);
  DelMatrix(&models);
}


/* Variable selection using Genetic algorithm
 *
 *
 */

int FitnessMaxsimized(dvector* modelsfitness, dvector* modelsfitness_old, double populationconvergence)
{
  size_t i, popsize = 0;

  /*
  puts("OLD");
  PrintDVector(modelsfitness_old);
  puts("NEW");
  PrintDVector(modelsfitness);

  sleep(1);
  */

  for(i = 0; i < modelsfitness->size; i++){
    if(FLOAT_EQ(getDVectorValue(modelsfitness, i), getDVectorValue(modelsfitness_old, i), 1e-3)){
      continue;
    }
    else{
      popsize++;
    }
  }

  if((double)popsize/(double)modelsfitness->size < populationconvergence){
    return 1;
  }
  else{
    return 0;
  }
}

void PLSGAVariableSelection(matrix *mx, matrix *my, matrix *px, matrix *py,
                       size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                       size_t population_size, double fraction_of_population, double mutation_rate, size_t crossovertype, double nswapping, double populationconvergence,
                      uivector **varselected, matrix **map, uivector **vardistribution, size_t nthreads, ssignal *s)
{
  size_t i, j, k, position, mutatebit, bitswap, totbitswapp, copy_selection, crossover_selection, mutation_selection, init, a, b, nvarON, blockitercount;
  double best_fitness = 0, best_bias = 99, model_r2, tmp_fitness, tmp_bias, tmp_model_r2;
  matrix *models1, *models2;
  dvector *modelsfitness, *best_modelsfitness, *model_f_distribution, *best_model_f_distribution, *rowmap;
  uivector *popvector, *best_model, *crosspoints, *selectedid;


  /* Initialize with rando models..*/
  NewMatrix(&models1, population_size, mx->col); /*models of the generation k */
  NewMatrix(&models2, population_size, mx->col); /*models of the generation k + 1*/
  NewDVector(&modelsfitness, population_size); /* Stored the Q^2 or R^2 of models*/
  NewDVector(&best_modelsfitness, population_size); /* Stored the Q^2 or R^2 of models*/
  NewDVector(&model_f_distribution, population_size); /* Stored the Q^2 or R^2 of models*/
  NewDVector(&best_model_f_distribution, population_size); /* Stored the Q^2 or R^2 of models*/
  NewUIVector(&best_model, mx->col);
  NewUIVector(&popvector, mx->col);

  for(i = 0; i < mx->col; i++){
    UIVectorAppend(vardistribution, 0);
  }

  NewDVector(&rowmap, 4); /* 1 = R2, 2 = Q2, 3 = f-test, 4 = number of variables ON */
  best_fitness = 0.f;

  copy_selection = (size_t)floor((1-fraction_of_population) * population_size);
  crossover_selection = (size_t)floor(fraction_of_population * population_size);
  mutation_selection = (size_t)floor(mutation_rate*population_size);

  if(crossovertype == 0){ /*single point crossover */
    totbitswapp = 1;
  }
  else if(crossovertype == 1){ /* two point crossover */
    totbitswapp = 2;
  }
  else{ /*uniform crossover*/
    if(nswapping > 1){
      totbitswapp = (size_t)floor(((double)models1->col * nswapping)/ 100.f);
    }
    else{
      totbitswapp = (size_t)floor((double)models1->col * nswapping);
    }

    if(totbitswapp > 0){
      if(totbitswapp % 2 != 0){
        totbitswapp +=1;
      }
    }
    else{
      totbitswapp = 1;
    }
  }

  NewUIVector(&crosspoints, totbitswapp);

  /*printf("copy %u cross %u muta %u bit to swap: %u\n", (size_t)copy_selection, (size_t)crossover_selection, (size_t)mutation_selection, (size_t)totbitswapp);*/
  init = (size_t)(mx->row + mx->col + my->col + population_size + copy_selection + mutation_selection);
  srand(init);

  DVectorSet(modelsfitness, 0.f);
  DVectorSet(best_modelsfitness, 0.f);
  DVectorSet(model_f_distribution, 0.f);
  DVectorSet(best_model_f_distribution, 0.f);
  MatrixSet(models1, 0.f);
  MatrixSet(models2, 0.f);

  for(i = 0; i < population_size; i++){
    for(j = 0; j < mx->col; j++){
      position = rand() % mx->col;
      if((size_t)getMatrixValue(models1, i, position) == 1){ /* set to OFF*/
        setMatrixValue(models1, i, position, 0.0);
      }
      else{ /* Set to ON*/
        setMatrixValue(models1, i, position, 1.0);
      }
    }

    /*COMPUTE THE COST OF ACTUAL MODEL RANDOM MODEL*/
    for(j = 0; j < models1->col; j++){
      setUIVectorValue(popvector, j, (size_t)getMatrixValue(models1, i, j));
      setUIVectorValue((*vardistribution), j, getUIVectorValue((*vardistribution),j)+(size_t)getMatrixValue(models1, i, j));
    }

    nvarON = 0;
    for(j = 0; j < popvector->size; j++){
      if(getUIVectorValue(popvector, j) == 1){
        nvarON++;
      }
    }


    if(nvarON > 0){
      PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &tmp_model_r2, &tmp_fitness, &tmp_bias, nthreads, s);
    }
    else{
      tmp_model_r2 = tmp_fitness = 0;
    }

    setDVectorValue(modelsfitness, i, tmp_fitness);
    setDVectorValue(model_f_distribution, i, (tmp_fitness * (population_size - mx->col -1)) /  (mx->col*(1-tmp_fitness)));

    if(tmp_fitness > best_fitness && tmp_bias < best_bias){
      best_fitness = tmp_fitness;
      model_r2 = tmp_model_r2;
      best_bias = tmp_bias;

      for(j = 0; j < models1->col; j++){
        setUIVectorValue(best_model, j, getUIVectorValue(popvector, j));
      }

      setDVectorValue(rowmap, 0, model_r2);
      setDVectorValue(rowmap, 1, best_fitness);
      setDVectorValue(rowmap, 2, getDVectorValue(model_f_distribution, i));
      setDVectorValue(rowmap, 3, (double)nvarON);

      MatrixAppendRow(map, rowmap);
      DVectorSet(rowmap, 0.f);
    }


    /*
    printf("Modello %u Generato\n", (unsigned int)i);
    PrintUIVector(popvector);
    */
  }

  /*
  puts("first best model");
  PrintUIVector(best_model);

  PrintMatrix((*map));

  puts("Models1");
  PrintMatrix(models1);
  puts("Models2");
  PrintMatrix(models2);
  */

  blockitercount = 0;
  do{
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      /*printf("blockitercount %u\n", (unsigned int) blockitercount);*/

      tmp_fitness = tmp_model_r2 = 0.f;
      tmp_bias = 99;
      /*Create generation k+1*/
      /*1. Copy
        *
        * Select (1-fraction_of_population) * population_size members of models1 and insert into the new population models2
        * through the roulette whell selection
        */
      if(crossover_selection % 2 != 0){
        crossover_selection -=1;
        copy_selection +=1;
      }

      initUIVector(&selectedid);
      StochasticUniversalSample(modelsfitness, copy_selection, init, &selectedid);
      /*RouletteWheelselection(modelsfitness, copy_selection, init, &selectedid);*/

      k = 0;
      for(i = 0; i < selectedid->size; i++, k++){
        position = getUIVectorValue(selectedid, i);
        for(j = 0; j < models1->col; j++){
          setMatrixValue(models2, i, j, getMatrixValue(models1, position, j));
        }
      }
      DelUIVector(&selectedid);


      /*
      * 2. Crossover
      *
      * Select fraction_of_population * population_size members of models1; pair them up; produce offspring; inser the offspring into models2
      *
      * The selection is made with the roulette whell selection
      *
      *
      */

      initUIVector(&selectedid);

      StochasticUniversalSample(modelsfitness, crossover_selection, init, &selectedid);
      /*RouletteWheelselection(modelsfitness, crossover_selection, init, &selectedid);*/

      for(i = 0; i < selectedid->size; i+=2){

        position = getUIVectorValue(selectedid, i);

        for(bitswap = 0; bitswap < totbitswapp; bitswap++){
          setUIVectorValue(crosspoints, bitswap, rand() % models1->col);
        }

        SortUIVector(crosspoints);

        if(totbitswapp >=2){
          /*offspring 1*/
          for(j = 0; j < models1->col; j++){
            for(bitswap = 0; bitswap < totbitswapp; bitswap += 2){
              if(j > getUIVectorValue(crosspoints, bitswap) && j < getUIVectorValue(crosspoints, bitswap+1)){ // set pos1
                setMatrixValue(models2, k, j, getMatrixValue(models1, position, j));
              }
              else{ // set pos 2
                setMatrixValue(models2, k, j, getMatrixValue(models1, position+1, j));
              }
            }
          }

          /*offspring 2*/
          k++;
          for(j = 0; j < models1->col; j++){
            for(bitswap = 0; bitswap < totbitswapp; bitswap += 2){
              if(j > getUIVectorValue(crosspoints, bitswap) && j < getUIVectorValue(crosspoints, bitswap+1)){ // set pos1
                setMatrixValue(models2, k, j, getMatrixValue(models1, position+1, j));
              }
              else{ // set pos 2
                setMatrixValue(models2, k, j, getMatrixValue(models1, position, j));
              }
            }
          }
        }
        else{
          /*offspring 1*/
          for(j = 0; j < models1->col; j++){
            if(j < getUIVectorValue(crosspoints, 0)){
              setMatrixValue(models2, k, j, getMatrixValue(models1, position, j));
            }
            else{
              setMatrixValue(models2, k, j, getMatrixValue(models1, position+1, j));
            }
          }

          /*offspring 2*/
          k++;
          for(j = 0; j < models1->col; j++){
            if(j < getUIVectorValue(crosspoints, 0)){
              setMatrixValue(models2, k, j, getMatrixValue(models1, position+1, j));
            }
            else{
              setMatrixValue(models2, k, j, getMatrixValue(models1, position, j));
            }
          }
        }
        k++;

      }
      DelUIVector(&selectedid);

    /*3. Mutate
      * Select a mutation_rate * population_size elements of models2 and invert a randomly selected bit in each;
      * The selection is always made by the roulette Whell Selection algorithm
      */

      initUIVector(&selectedid);

      StochasticUniversalSample(modelsfitness, mutation_selection, init, &selectedid);
  /*      RouletteWheelselection(modelsfitness, mutation_selection, init, &selectedid);*/

      for(i = 0; i < selectedid->size; i++){

        position = getUIVectorValue(selectedid, i);
        mutatebit = rand() % models1->col;

        if((size_t)getMatrixValue(models2, position, mutatebit) == 1){ /*set OFF*/
          setMatrixValue(models2, position, mutatebit, 0.0);
        }
        else{ /*set ON*/
          setMatrixValue(models2, position, mutatebit, 1.0);
        }
      }
      DelUIVector(&selectedid);


    /* Evaluate fittness for the new generation.....
      *
      * Copy the old fitness to
      *
      */

      for(i = 0; i < modelsfitness->size; i++){
        if(getDVectorValue(modelsfitness, i) > getDVectorValue(best_modelsfitness, i)){
          setDVectorValue(best_modelsfitness, i, getDVectorValue(modelsfitness, i));
        }

        if(getDVectorValue(model_f_distribution, i) > getDVectorValue(best_model_f_distribution, i)){
          setDVectorValue(best_model_f_distribution, i, getDVectorValue(model_f_distribution, i));
        }
      }

      for(i = 0; i < models2->row; i++){
        for(j = 0; j < models2->col; j++){
          setUIVectorValue(popvector, j, (size_t)getMatrixValue(models2, i, j));
          setUIVectorValue((*vardistribution), j, getUIVectorValue((*vardistribution),j)+(size_t)getMatrixValue(models1, i, j));
        }

        nvarON = 0;
        for(j = 0; j < models2->col; j++){
          if(getUIVectorValue(popvector, j) == 1){
            nvarON++;
          }
        }

        if(nvarON > 0){
          PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &tmp_model_r2, &tmp_fitness, &tmp_bias, nthreads, s);
        }
        else{
          tmp_model_r2 = tmp_fitness = 0;
          tmp_bias = 99;
        }

        setDVectorValue(modelsfitness, i, tmp_fitness); /* Stored the Q^2 or R^2 of models*/
        setDVectorValue(model_f_distribution, i, (tmp_fitness * (population_size - mx->col -1)) /  (mx->col*(1-tmp_fitness)));

        /*printf("[LIBSCIENTIFIC DEBUG] >> tmpfitness: %f bestfitness: %f tmpbias: %f bestbias: %f\n", tmp_fitness, best_fitness, tmp_bias, best_bias);*/

        if(tmp_fitness > best_fitness && tmp_bias < best_bias){

          best_bias = tmp_bias;
          best_fitness = tmp_fitness;
          model_r2 = tmp_model_r2;

          for(j = 0; j < models2->col; j++){
            setUIVectorValue(best_model, j, getUIVectorValue(popvector, j));
          }

          /*
          puts("POPVECTOR");
          PrintUIVector(popvector);

          puts("Models1");
          PrintMatrix(models1);
          puts("Models2");
          PrintMatrix(models2);

          sleep(1);
          */
          setDVectorValue(rowmap, 0, model_r2);
          setDVectorValue(rowmap, 1, best_fitness);
          setDVectorValue(rowmap, 2, getDVectorValue(model_f_distribution, i));
          setDVectorValue(rowmap, 3, nvarON);

          MatrixAppendRow(map, rowmap);
          DVectorSet(rowmap, 0.f);
        }
      }

      /* copy the new generation... */
      for(i = 0; i < models2->row; i++){
        for(j = 0; j < models2->col; j++){
          setMatrixValue(models1, i, j, getMatrixValue(models2, i, j));
        }
      }

      a = FitnessMaxsimized(modelsfitness, best_modelsfitness, populationconvergence);
      b = FitnessMaxsimized(model_f_distribution, best_model_f_distribution, populationconvergence);

      /*
      printf("a %d  b %d\n", (unsigned int)a, (unsigned int)b);
      */

      if(a == b && a == 1){
        if(blockitercount > 1000){
          /*No more good solutions found...
          * decrease the population size and the population convergence
          */
          break;
        }
        else{
          blockitercount++;
        }
      }
      else{
        blockitercount = 0;
      }
    }
  }while((a != 0 && b != 0));


  /*
  puts("BEST_MODEL");
  PrintUIVector(best_model);
  */

  for(i = 0; i < best_model->size; i++){
    UIVectorAppend(varselected, getUIVectorValue(best_model,  i));
  }

  /*
  puts("VARSELEVTED");
  PrintUIVector((*varselected));
  puts("----------");
  */

  DelDVector(&rowmap);
  DelUIVector(&crosspoints);
  DelDVector(&best_modelsfitness);
  DelDVector(&model_f_distribution);
  DelDVector(&best_model_f_distribution);
  DelDVector(&modelsfitness);
  DelUIVector(&popvector);
  DelUIVector(&best_model);
  DelMatrix(&models2);
  DelMatrix(&models1);
}

void PLSSpearmannVariableSelection(matrix* mx, matrix* my, matrix* px, matrix* py,
                                   size_t xautoscaling, size_t yautoscaling, size_t nlv, int validation_type, size_t ngroup, size_t niter,
                                   double threshold,
                                   uivector** varselected, matrix** map, uivector **vardistribution, size_t nthreads, ssignal* s){
  size_t i, j, k, nstep, nvarON;
  double step, r2, q2, best_q2 = 0, rangemax;

  matrix *mxy, *spearmanncorrelmx;
  uivector *popvector, *bestmodel;
  dvector *rowmap;

  NewDVector(&rowmap, 4);
  NewUIVector(&popvector, mx->col);
  NewUIVector(&bestmodel, mx->col);
  NewMatrix(&mxy, mx->row, mx->col+my->col);
  UIVectorResize(vardistribution, mx->col);

  for(i = 0; i < mx->row; i++){
    k = 0;
    for(j = 0; j < mx->col; j++){
      setMatrixValue(mxy, i, k, getMatrixValue(mx, i, j));
      k++;
    }

    for(j = 0; j < my->col; j++){
      setMatrixValue(mxy, i, k, getMatrixValue(my, i, j));
      k++;
    }
  }

  initMatrix(&spearmanncorrelmx);

  SpearmanCorrelMatrix(mxy, spearmanncorrelmx);

  rangemax = fabs(getMatrixValue(spearmanncorrelmx, 0, mx->col)); /* the first y value in the spearmanncorrelmx */
  for(j = 0; j < my->col; j++){
    for(i = 1; i < mx->col; i++){
      if(rangemax < getMatrixValue(spearmanncorrelmx, i, mx->col+j)){
        rangemax = getMatrixValue(spearmanncorrelmx, i, mx->col+j);
      }
      else{
        continue;
      }
    }
  }

  if(threshold < 0 || threshold > 1){
    threshold = 0.1; /* 1 / 10 */
  }

  threshold *= rangemax;
  nstep = (size_t) ceil(rangemax/threshold);
  step = rangemax - threshold; /* From maximum of relaction to minimum */
  for(i = 0; i < nstep; i++){
    if(s != NULL && (*s) == SIGSCIENTIFICSTOP){
      break;
    }
    else{
      nvarON = 0;
      for(k = 0; k < mx->col; k++){
        double nyes = 0.f;
        for(j = 0; j < my->col; j++){
          if(fabs(getMatrixValue(spearmanncorrelmx, mx->col+j, k)) > step){
            nyes += 1;
          }
          else{
            continue;
          }
        }

        if((double)((my->col - nyes)/(double)my->col) > 0.5){ /* if the variable is no more good for the 50% of y reject*/
          setUIVectorValue(popvector, k, 0);
        }
        else{
          setUIVectorValue(popvector, k, 1);
          setUIVectorValue((*vardistribution), k, getUIVectorValue((*vardistribution), k)+1);
          nvarON++;
        }
      }

      if(nvarON > 0){
        PLSCostPopulation(mx, my, px, py, popvector, xautoscaling, yautoscaling, nlv, validation_type, ngroup, niter, &r2, &q2, NULL, nthreads, s);


        if(q2 > best_q2){
          for(j = 0; j < popvector->size; j++){
            setUIVectorValue(bestmodel, j, getUIVectorValue(popvector, j));
          }
          best_q2 = q2;
        }

        setDVectorValue(rowmap, 0, r2);
        setDVectorValue(rowmap, 1, q2);
        setDVectorValue(rowmap, 2, step);
        setDVectorValue(rowmap, 3, nvarON); /* number of variables */
        MatrixAppendRow(map, rowmap);
        DVectorSet(rowmap, 0);
        UIVectorSet(popvector, 0);
      }
      step -= threshold;
    }
  }

  for(i = 0; i < bestmodel->size; i++){
    UIVectorAppend(varselected, getUIVectorValue(bestmodel,  i));
  }

  DelUIVector(&bestmodel);
  DelDVector(&rowmap);
  DelUIVector(&popvector);
  DelMatrix(&spearmanncorrelmx);
  DelMatrix(&mxy);
}
