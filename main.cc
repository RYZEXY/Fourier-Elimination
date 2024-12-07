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
#include "main.h"
using namespace std;

void trivial(vector<vector<double>> coefficients, vector<double> upper_bound, vector<double>&x) { //not left to the reader as an exercise - the reader can sit back and relax!

   //find the smallest upperbound adn the largest lower bound. if lower>upper then infeasible. otherwise, pick x in between upper and lower. 

    // Initialize min_upper_bound to positive infinity
    double min_upper_bound = numeric_limits<double>::infinity();
    // Initialize max_lower_bound to negative infinity
    double max_lower_bound = -numeric_limits<double>::infinity();

    bool hasLowerBound = false;
    bool hasUpperBound = false;

    // Find the smallest upper bound and the largest lower bound starting from index 1
    for (size_t i = 0; i < coefficients[0].size(); i++) {
        if (coefficients[0][i] < 0) {
            double lower_bound = upper_bound[i] / (coefficients[0][i]);
            hasLowerBound = true;

            if (lower_bound > max_lower_bound) {
                max_lower_bound = lower_bound;
            }
        } else if (coefficients[0][i] == 0) {
            // Handle the case where the coefficient is zero (no dependence on x)
            if (upper_bound[i] < 0) {
                cout << "The system is infeasible." << endl;
                exit(1);
            }
        }
        else {
            double upper_bound_value = upper_bound[i];
            if (!isfinite(upper_bound_value)) {
                upper_bound_value = max_lower_bound;
            } else {
                hasUpperBound = true;
            }

            if (upper_bound_value < min_upper_bound) {
                min_upper_bound = upper_bound_value;
            }
        }
    }

    if (hasLowerBound && hasUpperBound && max_lower_bound > min_upper_bound) {
        cout << "The system is infeasible." << endl;
        exit(1);
    } 
    else {
        if (!hasLowerBound) {
            max_lower_bound = min_upper_bound;
        } else if (!hasUpperBound) {
            min_upper_bound = max_lower_bound;
        }
        
        cout << "The minimum upper bound is: " << min_upper_bound << " and the maximum lower bound is: " << max_lower_bound << endl;
        double x_value = (max_lower_bound + min_upper_bound) / 2;
        x.emplace_back(x_value);
    }
}


void notTrivial(vector<double>& x, list<vector<vector<double>>>& system, list<vector<double>>& allBounds) {
    /*
    Multiply the respective columns by the x value for each then subtract the upper_bound by x.
    The order of x is the same order as coefficients. i.e., the first element in x is x1, which will multiply with the first column in coefficients (matrix multiplication).
    Then, call trivial() with the last column of coefficient. trivial will push x_n to the end of x.
    */
   int size = system.size();
    for (int i = 0; i < size; i++) {
        vector<vector<double>>& currentCoeff = system.front();
        vector<double>& currentBounds = allBounds.front();

            cout << "New system is:" << endl;

        vector<double> oldBounds=currentBounds;


         for (size_t i = 0; i < currentCoeff.size(); i++) {
             for (size_t j = 0; j < currentCoeff[i].size(); j++) {
                 if(j!=currentCoeff[i].size()-1){
                     cout << currentCoeff[i][j] <<"("<< x[j]<<")"<<"+" << ' ';
                 }
                 else{
                     cout << currentCoeff[i][j] <<"x"<<j+1 << ' ';
                 }
             }
             cout << "<= " << oldBounds[i] << endl;
         }
         cout<<endl;
        // Multiply the columns of the coefficient matrix by the values in x then subtract the constant A(n)x(n) from b(n)
        for (size_t j = 0; j < currentCoeff.size(); j++) {
            for (size_t k = 0; k < currentCoeff[j].size() - 1; k++) {
                currentCoeff[j][k] *= x[k];
                 currentBounds[j] -= currentCoeff[j][k];
            }
           
        }


            cout << "New updated system is:" << endl;


            for (size_t i = 0; i < currentCoeff.size(); i++) {
                for (size_t j = 0; j < currentCoeff[i].size(); j++) {
                    if(j!=currentCoeff[i].size()-1){
                    }
                    else{
                        cout << currentCoeff[i][j] <<"x"<<j+1 << ' ';
                    }
                }
                cout << "<= " << currentBounds[i] << endl;
            }
            cout<<endl;

        // Call the trivial function with the last column of coefficients
        vector<vector<double>> lastColumn(1, vector<double>(currentCoeff.size()));
        for (size_t j = 0; j < currentCoeff.size(); j++) {
            lastColumn[0][j] = currentCoeff[j].back();
        }


        trivial(lastColumn, currentBounds, x);

        // Remove the processed system and bounds from the lists
        system.pop_front();
        allBounds.pop_front();
    }
}


void newSystem(vector<vector<double>>& coefficients, vector<double>& upper_bound, list<vector<vector<double>>>& system, list<vector<double>>& allBounds) {
    vector<vector<double>> newCoefficients;
    vector<double> newUpperBound;

    vector<vector<double>> positiveCoeffs; // Store positive coefficients (or 0)
    vector<double> positiveBounds;         // Store corresponding bounds

    vector<vector<double>> negativeCoeffs; // Store negative coefficients (or 0)
    vector<double> negativeBounds;         // Store corresponding bounds

    for (int i = 0; i < coefficients.size(); i++) {
        if (coefficients[i][coefficients[0].size() - 1] == 0) { // J0
            vector<double> newRow(coefficients[i].begin(), coefficients[i].end() - 1); // Remove the last column
            newCoefficients.emplace_back(newRow);
            newUpperBound.emplace_back(upper_bound[i]);
        } else if (coefficients[i][coefficients[0].size() - 1] < 0) { // J-
            negativeCoeffs.emplace_back(coefficients[i]);
            negativeBounds.emplace_back(upper_bound[i]);
        } else { // J+
            positiveCoeffs.emplace_back(coefficients[i]);
            positiveBounds.emplace_back(upper_bound[i]);
        }
    }

    // Check if either J+ or J- is empty
    if (negativeCoeffs.empty()) {
        // Add positive coefficients to cancel the last column
        vector<double> newRow(coefficients[0].size() - 1, 0.0);
        for (const auto& row : positiveCoeffs) {
            for (int k = 0; k < row.size() - 1; k++) {
                newRow[k] += row[k];
            }
        }
        newCoefficients.emplace_back(newRow);
        
        // Sum the upper bounds from J+
        double upperSum = 0.0;
        for (size_t i = 0; i < positiveCoeffs.size(); i++) {
            upperSum += positiveBounds[i];
        }
        newUpperBound.emplace_back(upperSum);
    } 
    else if (positiveCoeffs.empty()) {
        // Subtract negative coefficients to cancel the last column
        vector<double> newRow(coefficients[0].size() - 1, 0.0);
        for (const auto& row : negativeCoeffs) {
            for (int k = 0; k < row.size() - 1; k++) {
                newRow[k] -= row[k];
            }
        }
        newCoefficients.emplace_back(newRow);
        
        // Sum the upper bounds from J-
        double upperSum = 0.0;
        for (size_t i = 0; i < negativeCoeffs.size(); i++) {
            upperSum -= negativeBounds[i]; // Subtract the upper bounds
        }
        newUpperBound.emplace_back(upperSum);
    }
    else {
        // Both J+ and J- are non-empty, create new systems by adding each combination
        for (size_t i = 0; i < positiveCoeffs.size(); i++) {
            for (size_t j = 0; j < negativeCoeffs.size(); j++) {
                vector<double> newRow(coefficients[0].size() - 1, 0.0);
                for (int k = 0; k < newRow.size(); k++) {
                    newRow[k] = positiveCoeffs[i][k] + negativeCoeffs[j][k];
                }
                newCoefficients.emplace_back(newRow);

                // Calculate the upper bound as the sum of individual upper bounds
                double upperSum = positiveBounds[i] + negativeBounds[j];
                newUpperBound.emplace_back(upperSum);
            }
        }
    }
    // Print the coefficients and corresponding upper bounds before simplification
    cout << "New system is:" << endl;

    for (size_t i = 0; i < newCoefficients.size(); i++) {
        for (size_t j = 0; j < newCoefficients[i].size(); j++) {
            if(j!=newCoefficients[i].size()-1){
                cout << newCoefficients[i][j] <<"x"<<j+1 <<" +" << ' ';
            }
            else{
                cout << newCoefficients[i][j] <<"x"<<j+1 << ' ';
            }
        }
        cout << "<= " << newUpperBound[i] << endl;
    }
    cout<<endl;

    // Simplify the new coefficients and bounds
    simplifyCoefficients(newCoefficients, newUpperBound, system, allBounds);
}



void simplifyCoefficients(vector<vector<double>>& newCoefficients, vector<double>& newUpperBound, list<vector<vector<double>>>& system, list<vector<double>>& allBounds) {
    
    for (size_t i = 0; i < newCoefficients.size(); ++i) {
        if (newCoefficients[i].empty()) {
            continue;
        }

        double lastCoefficient = newCoefficients[i].back();
        double divisor = abs(lastCoefficient);

        if (divisor != 0) {
            for (size_t j = 0; j < newCoefficients[i].size(); ++j) {
                newCoefficients[i][j] /= divisor;
            }

            newUpperBound[i] /= divisor;
        } 
    }

    // Print the simplified coefficients
    cout<<"This is simplified to: "<<endl;

    for (size_t i = 0; i < newCoefficients.size(); i++) {
        for (size_t j = 0; j < newCoefficients[i].size(); j++) {
            if(j!=newCoefficients[i].size()-1){
                cout << newCoefficients[i][j] <<"x"<<j+1 <<" +" << ' ';
            }
            else{
                cout << newCoefficients[i][j] <<"x"<<j+1 << ' ';
            }
        }
        cout << "<= " << newUpperBound[i] << endl;
    }
    cout<<endl;


    // Push the simplified coefficients to the system
    system.push_front(newCoefficients);
    allBounds.push_front(newUpperBound);
    if (newCoefficients[0].size() == 1) {
        return;  // Exit the function when the trivial case is reached
    }
    else{
    newSystem(system.front(),allBounds.front(),system,allBounds);
    }
}





int main(int argc, char *argv[]) {
    // if (argc != 2) {
    //     cerr << "Usage: " << argv[0] << " <input_file>" << endl;
    //     return 1;
    // }

    ifstream inputFile;
    inputFile.open("test2.txt");

    if (!inputFile) {
        cerr << "Error opening the file." << endl;
        return 1;
    }

    int numRows, numCols;
    inputFile >> numRows >> numCols;

    vector<vector<double>> coefficients(numRows, vector<double>(numCols));
    list<vector<vector<double>>> system;
    vector<double> upper_bound(numRows);
    list<vector<double>> allBounds;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            if (!(inputFile >> coefficients[i][j])) {
                cerr << "Error reading double values from the file." << endl;
                return 1;
            }
        }
    }
    for(int i=0; i < numRows; i++){
        if (!(inputFile >> upper_bound[i])) {
            cerr << "Error reading double values from the file." << endl;
            return 1;
        }
    }
    inputFile.close();

    // system.push_back(coefficients);
    // allBounds.push_back(upper_bound);


    for (size_t i = 0; i < coefficients.size(); i++) {
        for (size_t j = 0; j < coefficients[i].size(); j++) {
            if(j!=coefficients[i].size()-1){
                cout << coefficients[i][j] <<"x"<<j+1 <<" +" << ' ';
            }
            else{
                cout << coefficients[i][j] <<"x"<<j+1 << ' ';
            }
        }
        cout << "<= " << upper_bound[i] << endl;
    }
    cout<<endl;


    vector<double> x;
    simplifyCoefficients(coefficients,upper_bound,system,allBounds);
    
     vector<vector<double>> lastColumn(1, vector<double>(system.front().size()));
     for (size_t j = 0; j < system.front().size(); j++) {
         lastColumn[0][j] = system.front()[j].back();
     }
    trivial(lastColumn,allBounds.front(),x);
    system.pop_front();
    allBounds.pop_front();
    if(!system.empty()){
        notTrivial(x, system, allBounds);
    }

    for(int i=0; i<x.size(); i++){
        cout<< "x"<<i+1<<" = " << x[i] << " ";
    }
    return 0;
}
