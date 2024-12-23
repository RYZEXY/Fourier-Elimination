This program is an implementation of the Fourier-Motzkin Elimination Algorithm, a technique used for solving systems of linear inequalities. The algorithm eliminates variables systematically, repeatedly
reducing a system of inequalities in n variables to a system in n−1 variables. Once we reach a base case of one variable, it is trivial to determine a solution range for that one variable. We can then recursively 
solve for the feasibility range for vectors of size n+1 until we reach the size of the original vector. Ultimately, this program provides a finite algorithm for solving systems of linear inequalities, which has
practical use for linear programming, and mathematical optimization problems.
