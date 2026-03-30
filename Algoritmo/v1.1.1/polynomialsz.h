#ifndef POLYNOMIALSZ_H_INCLUDED
#define POLYNOMIALSZ_H_INCLUDED

/**
 * @file polynomialsz.h
 * @brief Header file for polynomial factorization structures.
 *
 * This file provides the necessary data structures and declarations used by the
 * polynomial factorization algorithm, which supports:
 *  - Aberth-Ehrlich method;
 *  - Bhaskara's method (quadratic formula);
 *  - Briot-Ruffini (synthetic division);
 *  - Cyclotomic Factorization.
 *
 * The algorithm works with polynomials with integer coefficients.
 *
 * @author Isaque Passos
 * @version 1.1.1
 * @date 2026
 *
 * @note Version 1.0.0: Initial implementation.
 * @note Version 1.1.0: Added Aberth method and other improvements.
 * @note Version 1.1.1: Mild corrections and improvements.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> /**< Required for mathematical operations. */

#define _USE_MATH_DEFINES /**< Enables use of M_PI and other constants in some compilers. */

// polynomialsz.h

/**
 * @brief Symbolic variable used in polynomial expressions (default: 'x')
 */
extern char var;

/**
 * @brief Exponent of the factor x extracted by divideX()
 */
extern int degreeX;

/**
 * @brief GCD factor extracted by divideGCD()
 */
extern int divider;

/**
 * @brief Flag indicating if factorization was complete (1) or not (0)
 */
extern int sol;

/**
 * @struct term
 * @brief Represents a single term in a polynomial.
 *
 * A term is composed of a coefficient and an exponent.
 * For example, the term 3x² has coefficient = 3 and exponent = 2.
 */
typedef struct sterm {
    int coefficient; /**< The coefficient of the term. */
    int exponent;    /**< The exponent of the term. */
} term;

/**
 * @struct polynomial
 * @brief Represents a polynomial as an array of terms.
 *
 * The terms are not required to be sorted, but usually are ordered by descending exponents.
 */
typedef struct spolynomial {
    term *terms;        /**< Dynamic array of polynomial terms. */
    int numTerms;   /**< Number of terms in the polynomial. */
} polynomial;

void printCyclotomicRoots(int N, polynomial p);
int findN(polynomial p);
int dividesXPowerNMinusOne(polynomial p, int N);

//-----------------------------------------------------------------------------
/**
 * @brief Creates and sets the coefficient and exponent of a polynomial term.
 *
 * @param coef The coefficient of the term.
 * @param exp The expoent of the term.
 * @return T The complete term.
 */
term setTerms(int coef, int exp);
//-----------------------------------------------------------------------------
/**
 * @brief Returns the base of a power with an exponent equal to ind, which equals to rad.
 *
 * @param rad The radicand (power's result).
 * @param ind The index (power's exponent).
 * @return base The base (radical's result).
 */
double nrt(double rad, int ind);
//-----------------------------------------------------------------------------
/**
 * @brief Returns the exponent (index) of a power that has base rt and result rad.
 *
 * @param rad The radicand (power's result).
 * @param rt The root (power's base).
 * @return ind The index (power's exponent).
 */
int indOfRoot(int rad, double rt);
//-----------------------------------------------------------------------------
/**
 * @brief Creates and returns a polynomial with a specified number of terms.
 *
 * This function allocates memory for a polynomial structure and initializes
 * its internal array of terms based on t    // Formato Wolfram: remove parênteses e coeficiente 1he specified number of terms. It is
 * used as a starting point for defining or manipulating polynomials.
 *
 * @param numTerms The number of terms the polynomial should contain.
 * @return A polynomial structure with allocated space for the specified number of terms.
 */
polynomial pCreate(int numTerms);
//-----------------------------------------------------------------------------
/**
 * @brief Prints the polynomial.
 *
 * @param p The polynomial.
 */
void pPrint(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Calculates the greatest common divisor (GCD) of two integers using the Euclidean algorithm.
 *
 * @param a The first integer.
 * @param b The second integer.
 * @return The greatest common divisor of a and b.
 */
int gcd(int a, int b);
//-----------------------------------------------------------------------------
/**
 * @brief Generates a string representation of the factored form for a quadratic
 *        equation with integer (rational) roots and positive discriminant.
 *
 * @param aexp The exponent of the first term (2 or 4, indicating linear or quadratic factors)
 * @param den The denominator value (2a) from Bhaskara's formula
 * @param b The coefficient b from the quadratic equation
 * @param delta The discriminant value (b² - 4ac)
 * @param absDelta The absolute value of delta
 * @param bSimplify Pre-allocated string buffer to store the result
 * @param powerRoot String indicating root exponent format ("^2" for quartic, empty for quadratic)
 * @return Pointer to the formatted string buffer
 */
char* intSrPositiveDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Generates a string representation of the factored form for a quadratic
 *        equation with integer (rational) roots and negative discriminant.
 *
 * @param aexp The exponent of the first term (2 or 4, indicating linear or quadratic factors)
 * @param den The denominator value (2a) from Bhaskara's formula
 * @param b The coefficient b from the quadratic equation
 * @param delta The discriminant value (b² - 4ac)
 * @param absDelta The absolute value of delta
 * @param bSimplify Pre-allocated string buffer to store the result
 * @param powerRoot String indicating root exponent format ("^2" for quartic, empty for quadratic)
 * @return Pointer to the formatted string buffer
 */
char* intSrNegativeDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Generates a string representation of the factored form for a quadratic
 *        equation with irrational roots and positive discriminant.
 *
 * @param aexp The exponent of the first term (2 or 4, indicating linear or quadratic factors)
 * @param den The denominator value (2a) from Bhaskara's formula
 * @param b The coefficient b from the quadratic equation
 * @param delta The discriminant value (b² - 4ac)
 * @param absDelta The absolute value of delta
 * @param bSimplify Pre-allocated string buffer to store the result
 * @param powerRoot String indicating root exponent format ("^2" for quartic, empty for quadratic)
 * @return Pointer to the formatted string buffer
 */
char* floatSrPositiveDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Generates a string representation of the factored form for a quadratic
 *        equation with complex irrational roots and negative discriminant.
 *
 * @param aexp The exponent of the first term (2 or 4, indicating linear or quadratic factors)
 * @param den The denominator value (2a) from Bhaskara's formula
 * @param b The coefficient b from the quadratic equation
 * @param delta The discriminant value (b² - 4ac)
 * @param absDelta The absolute value of delta
 * @param bSimplify Pre-allocated string buffer to store the result
 * @param powerRoot String indicating root exponent format ("^2" for quartic, empty for quadratic)
 * @return Pointer to the formatted string buffer
 */
char* floatSrNegativeDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Determines the appropriate formatting function for rational roots
 *        based on discriminant sign.
 *
 * @param aexp The exponent of the first term (2 or 4)
 * @param den The denominator value (2a) from Bhaskara's formula
 * @param b The coefficient b from the quadratic equation
 * @param delta The discriminant value (b² - 4ac)
 * @param absDelta The absolute value of delta
 * @param bSimplify Pre-allocated string buffer to store the result
 * @param powerRoot String indicating root exponent format
 * @return Pointer to the formatted string buffer
 */
char* rationalSqRoots(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Determines the appropriate formatting function for irrational roots
 *        based on discriminant sign.
 *
 * @param aexp The exponent of the first term (2 or 4)
 * @param den The denominator value (2a) from Bhaskara's formula
 * @param b The coefficient b from the quadratic equation
 * @param delta The discriminant value (b² - 4ac)
 * @param absDelta The absolute value of delta
 * @param bSimplify Pre-allocated string buffer to store the result
 * @param powerRoot String indicating root exponent format
 * @return Pointer to the formatted string buffer
 */
char* irrationalSqRoots(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Main simplification function that coordinates rational/irrational
 *        root formatting based on discriminant properties.
 *
 * This function checks if the discriminant is a perfect square to determine
 * whether roots are rational or irrational, then calls the appropriate
 * formatting functions.
 *
 * @param aexp The exponent of the first term (2 or 4)
 * @param den The denominator value (2a) from Bhaskara's formula
 * @param b The coefficient b from the quadratic equation
 * @param delta The discriminant value (b² - 4ac)
 * @param absDelta The absolute value of delta
 * @param rootsPair Array containing the two calculated roots
 * @return Pointer to the formatted string buffer (must be freed by caller)
 */
char* bhaskaraSimplify(int aexp, int den, int b, double delta, int absDelta, double* rootsPair);
//-----------------------------------------------------------------------------
/**
 * @brief Calculates and prints the factorization of a second-degree polynomial
 *        using the Bhaskara (quadratic) formula.
 *
 * This function simplifies the result depending on the values ​​of "b" and delta.
 *
 * @param p A polynomial of degree 2 or 4 with three terms in descending order
 *          of exponent.
 */
void bhaskara(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Factors a polynomial using the Rational Root Theorem and synthetic division.
 *
 * Applies the Briot-Ruffini method to find integer roots in the range
 * [-√maxNum, √maxNum]. When a quadratic factor is found, calls bhaskara()
 * to complete factorization (including complex roots).
 *
 * @param p Polynomial to be factored (modified during the process)
 * @param maxNum Largest absolute coefficient value of the original polynomial
 *
 * @note The input polynomial is reduced as roots are found.
 * @note Global variable sol indicates if factorization was complete.
 */
void briotRuffini(polynomial p, int maxNum);
//-----------------------------------------------------------------------------
/**
 * @brief Main factorization function that coordinates all strategies.
 *
 * Orchestrates the complete factorization process in the following order:
 * 1. Extracts x factor (if any) via divideX() [recursive]
 * 2. Extracts coefficient GCD via divideGCD()
 * 3. If quadratic/biquadratic: uses bhaskara()
 * 4. If coefficients are sequential: uses briotRuffini()
 * 5. If coefficients are ±1: uses cyclotomicFac()
 * 6. Otherwise: uses aberth() for numerical approximation
 *
 * @param p Polynomial to be factored
 *
 * @note Global variable sol indicates if factorization was complete
 * @note Global variables divider and degreeX store extracted factors
 */
void fac(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Prints the polynomial followed by its factorization.
 *
 * @param p Polynomial to be printed and factored
 */
void printFac(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Unlike `printFac()`, this alternative also deallocates the polynomial's terms
 * memory using `free(p.terms)`.
 *
 * @param p The polynomial to be factored.
 */
 void pFree(polynomial p);
 //-----------------------------------------------------------------------------
 /**
 * @brief Extracts the factor x^(minExp) from the polynomial and calls fac() recursively.
 *
 * Finds the smallest exponent among all terms, divides all terms by x^(minExp),
 * and calls fac() on the reduced polynomial. The extracted factor is stored
 * in the global variable degreeX for later printing.
 *
 * @param p Original polynomial (not modified directly)
 *
 * @note Global variable degreeX stores the exponent of the extracted x factor
 */
void divideX(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Divides all polynomial coefficients by their GCD.
 *
 * Calculates the greatest common divisor (GCD) of all polynomial coefficients
 * and divides each coefficient by this value. The extracted factor is stored
 * in the global variable divider for later printing.
 *
 * @param p Polynomial to be simplified (modified in-place)
 *
 * @note Global variable divider stores the extracted GCD
 */
void divideGCD(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Removes zero-coefficient terms from the polynomial and reallocates memory.
 *
 * Scans the term array, keeping only those with coefficient != 0.
 * If all terms are zero, completely frees the polynomial.
 * Reallocates the array to the exact needed size.
 *
 * @param p Pointer to the polynomial to be cleaned (modified in-place)
 */
void removeZeros(polynomial *p);
//-----------------------------------------------------------------------------
/**
 * @brief Factors cyclotomic polynomials (coefficients ±1).
 *
 * For polynomials of the form x^n + x^(n-1) + ... + x + 1 (or variations with ±1),
 * finds complex roots on the unit circle: exp(2πik/n).
 * Prints linear complex factors using exponential notation.
 *
 * @param p Cyclotomic polynomial to be factored
 *
 * @note Assumes polynomial has coefficients only ±1
 * @note Prints factors in format (x - Exp[i*π/m])
 */
void cyclotomicFac(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Approximates all roots of a polynomial using Aberth's method.
 *
 * Implements the Aberth (or Aberth–Ehrlich) method for simultaneous approximation
 * of all roots (real and complex) of a polynomial.
 *
 * Algorithm:
 * 1. Converts coefficients to double
 * 2. Estimates maximum root radius: R = 1 + max|coef[i]/coef[0]|
 * 3. Initializes roots equally spaced on circle of radius R
 * 4. Iterates Aberth correction until convergence:
 *    z_k^(new) = z_k - [P(z_k)/P'(z_k)] / [1 - (P(z_k)/P'(z_k)) * Σ_{j≠k} 1/(z_k - z_j)]
 *
 * @param p Polynomial whose roots will be approximated
 *
 * @note Uses complex numbers from <complex.h>
 * @note Convergence precision: 1e-12
 * @note Maximum of 30 iterations
 */
void aberth(polynomial p);

#endif // POLYNOMIALSZ_H_INCLUDED
