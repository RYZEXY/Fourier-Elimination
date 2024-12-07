#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <float.h>
#include <cstdlib>
#include <list>

/**
 * @brief solves the trivial system of inequalities with only 1 variable and appends it to the front of x
 *
 * @param coefficients Coefficient Matrix
 * @param upper_bound The inclusive upper bound for each equation
 * @param x Solution vector
 */
void trivial(std::vector<std::vector<double>> coefficients, std::vector<double> upper_bound, std::vector<double>& x);

/**
 * @brief Reduces a system of linear inequalities to a trivial inequality
 *
 * @param x Solution vector
 * @param system Data structure containing all created systems of linear inequalities
 * @param allBounds Data structure containing all upperbounds to it's respective system of linear inequalities
 */
void notTrivial(std::vector<double>& x, std::list<std::vector<std::vector<double>>>& system, std::list<std::vector<double>>& allBounds);

/**
 * @brief Reduces the leading coefficient to 1
 *
 * @param newCoefficients Newest coefficient matrix
 * @param newUpperBound Newest vector of upperbounds
 * @param system Data structure containing all created systems of linear inequalities
 * @param allBounds Data structure containing all upperbounds to it's respective system of linear inequalities
 */
void simplifyCoefficients(std::vector<std::vector<double>>& newCoefficients, std::vector<double>& newUpperBound, std::list<std::vector<std::vector<double>>>& system, std::list<std::vector<double>>& allBounds);

/**
 * @brief Creates a new system of linear inequalities by eliminating the last variable in the previous system
 *
 * @param coefficients Coefficient matrix
 * @param upper_bound Vector of upperbounds
 * @param system Data structure containing all created systems of linear inequalities
 * @param allBounds Data structure containing all upperbounds to it's respective system of linear inequalities
 */
void newSystem(std::vector<std::vector<double>>& coefficients, std::vector<double>& upper_bound, std::list<std::vector<std::vector<double>>>& system, std::list<std::vector<double>>& allBounds);

#endif