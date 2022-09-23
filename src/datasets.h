/* datasets.h
*
* Copyright (C) <2016>  Giuseppe Marco Randazzo
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef DATASET_H
#define DATASET_H
#include "matrix.h"
#include "scientificinfo.h"

/*
 * iris fisher dataset
 * - sepal length in cm
 * -sepal width in cm
 * -petal length in cm
 * -petal width in cm
 *
 * y variable 1
 * class names:
 * 0 Iris Setosa
 * 1 Iris Versicolour
 * 2 Iris Virginica
 */

void iris(matrix *x, matrix *y);

/*
 * Residential Building Dataset
 * X variables 13
 * X variable names
 * - crim
 * - zn
 * -indus
 * -chas
 * -nox
 * -rm
 * -age
 * -dis
 * -rad
 * -tax
 * -ptratio
 * -b
 * -lstat
 *
 * Y variable 1
 * Y variable name
 * -medv
 *
 * doi: 10.1016/0095-0696(78)90006-2.
 */

void boston_house_price(matrix *x, matrix *y);

/*
 * Residential Building Data set
 * X features:
 * 8 project physical and financial variables,
 * 19 economic variables and indices in 5 time lag numbers (5*19 = 95),
 * Y features
 * two output variables that are construction costs and sale prices
 * doi: doi.org/10.1061/(ASCE)CO.1943-7862.0001047
 *
 * X variable names:
 * - Project locality defined in terms of zip codes
 * - Total floor area of the building
 * - Lot area
 * - Total preliminary estimated construction cost based on the prices at the beginning of the project
 * - Preliminary estimated construction cost based on the prices at the beginning of the project
 * - Equivalent preliminary estimated construction cost based on the prices at the beginning of the project in a selected base year a
 * - Duration of construction
 * - The number of building permits issued
 * - Building services index (BSI) b for a preselected base year a
 * - Wholesale price index (WPI) c of building materials for the base year
 * - Total floor areas of building permits issued by the city/municipality
 * - Cumulative liquidity d
 * - Private sector investment in new buildings
 * - Land price index for the base year a
 * - The number of loans extended by banks in a time resolution e
 * - The amount of loans extended by banks in a time resolution e
 * - The interest rate for loan in a time resolution e
 * - The average construction cost of buildings by private sector at the time of completion of construction
 * - The average of construction cost of buildings by private sector at the beginning of the construction
 * - Official exchange rate with respect to dollars
 * - Nonofficial (street market) exchange rate with respect to dollars h
 * - Consumer price index (CPI) i in the base year a
 * - CPI of housing, water, fuel & power in the base year a
 * - Stock market index j
 * - Population of the city
 * - Gold price per ounce
 *
 *
 * Y variable names:
 * - Actual sales prices (output)
 * - Actual construction costs (output)
 *
 */
void residential_building(matrix *x, matrix *y);

#endif
