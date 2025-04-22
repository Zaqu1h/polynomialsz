#ifndef POLYNOMIALSZ_H_INCLUDED
#define POLYNOMIALSZ_H_INCLUDED

/**
 * @file polynomialsz.h
 * @brief Header file for polynomial factorization structures.
 *
 * This file provides the necessary data structures and declarations used by the
 * polynomial factorization algorithm, which supports:
 *  - Eisenstein's Criterion
 *  - Briot-Ruffini (synthetic division)
 *  - Bhaskara's method (quadratic formula)
 *  - Power difference factoring
 *  - Complex root analysis
 *
 * The algorithm works with polynomials with integer coefficients.
 *
 * @author Isaque Passos
 * @version 1.0.0
 * Initial public release. Internal versions and tests existed prior to this.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> /**< Required for mathematical operations. */

#define _USE_MATH_DEFINES /**< Enables use of M_PI and other constants in some compilers. */

/**
 * @brief Global symbolic variable used in polynomial expressions.
 *
 * This variable is used for displaying polynomial roots symbolically.
 * For example, if var = 'x', the roots of x² - 2 are printed as x + √2 and x - √2.
 */
extern char var;

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
    int num_of_terms;   /**< Number of terms in the polynomial. */
} polynomial;

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
 * its internal array of terms based on the specified number of terms. It is
 * used as a starting point for defining or manipulating polynomials.
 *
 * @param num_of_terms The number of terms the polynomial should contain.
 * @return A polynomial structure with allocated space for the specified number of terms.
 */
polynomial pCreate(int num_of_terms);
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
 * @brief Returns the n-th isPrime number, where n is less than or equal to 41.
 *
 * This function only supports values of n up to 41 because the 41st isPrime number is 179,
 * which is the largest isPrime whose square (179^2 = 32,041) still fits within a 32-bit signed integer.
 *
 * @param num The position of the isPrime number to return (must be ≤ 41).
 * @return The n-th isPrime number, or 0 if num > 41.
 */
int isPrime(int num);
//-----------------------------------------------------------------------------
/**
 * @brief Applies Eisenstein's criterion to determine whether a polynomial is irreducible.
 *
 * Checks if the polynomial satisfies the conditions of Eisenstein's Criterion for irreducibility
 * using the isPrime number p. According to the criterion, all coefficients except the leading one
 * must be divisible by p, and the constant term must not be divisible by p^2.
 *
 * Although Eisenstein's Criterion is a powerful tool, it has limitations — for example,
 * the polynomial (x^4 + 1) is irreducible over the rationals, despite not satisfying the criterion.
 *
 * @param x The polynomial to be tested.
 * @param p A isPrime number used for the Eisenstein's test.
 * @return 1 if the polynomial is irreducible by Eisenstein's criterion,
 *        -1 otherwise.
 */
int eisenstein(polynomial x, double p);
//-----------------------------------------------------------------------------
/**
 * @brief Prints all non-real complex roots of a polynomial, assuming unit modulus.
 *
 * This function assumes that the polynomial has only complex roots located on the unit circle
 * and calculates the complex roots in exponential form e^{iθ} and e^{-iθ}, where θ is an angle
 * in radians represented as a multiple of π. It handles both cases when the polynomial is
 * monic or not and adjusts the number of printed roots accordingly.
 *
 * The roots are displayed using Euler's formula and shown as factorized expressions.
 *
 * @param p The polynomial from which the complex roots will be extracted and printed.
 */
void complexp(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Generates a symbolic factorized string representation of a second-degree polynomial.
 *
 * This function receives the coefficients and discriminant of a quadratic equation and returns a
 * heap-allocated string with a symbolic factorization of its roots. It handles various cases, including:
 *  - Real and rational roots
 *  - Irrational roots (containing square roots)
 *  - Complex conjugate roots
 *  - Pure imaginary roots
 *
 * It also supports optional squaring of the variable (e.g., x^2) if the polynomial is a square
 * of a quadratic term (useful for degree 4 cases).
 *
 * @warning The returned string is dynamically allocated using `malloc()` and must be freed
 * by the caller to avoid memory leaks.
 *
 * @param aexp The exponent of the leading term (typically 2 or 4).
 * @param den The denominator used to normalize the roots (usually 2 * a).
 * @param b The coefficient of the linear term in the polynomial.
 * @param delta The discriminant value (b² - 4ac) of the polynomial.
 * @param rootsPair A pointer to an array of two doubles representing the numeric roots (used for formatting).
 *
 * @return A pointer to a dynamically allocated string containing the symbolic factorization.
 *         The string must be freed manually by the caller.
 */
char* bhaskaraSimplify(int aexp, int den, int b, double delta, double* rootsPair);
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
 * @brief Applies the sum or difference of powers identity (a^n ± b^n).
 *        For odd n, both sum and difference can be factorized.
 *        For even n, only the difference can be factorized.
 *
 * Assumes the polynomial is in the form x^n ± a and prints the factorization steps based on the formula.
 * Handles both even and odd exponents.
 *
 * @param p The polynomial to be factored.
 */
void powerDiffFactoring(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Attempts to factor a polynomial using the Rational Root Theorem and synthetic division.
 *
 * Applies Briot–Ruffini's method (synthetic division) to find integer and rational roots of a polynomial,
 * as well as approximated square root roots when possible. Displays factored linear and radical terms.
 *
 * @param p The polynomial to be factored.
 */
void briot_ruffini(polynomial p);

//-----------------------------------------------------------------------------

/**
 * @brief Prints the factorization of a polynomial or displays the appropriate message
 *        based on the type of factorization that can be performed on it.
 *
 * This function evaluates a polynomial to determine the appropriate factorization technique:
 * - If the polynomial is irreducible based on Eisenstein's criterion, it outputs a message
 *   indicating that there is no non-imaginary factorization or roots representable by radicals.
 * - If the polynomial is a quadratic, it applies potential difference factoring, if applicable.
 * - If the polynomial has an even degree and three or fewer terms, it applies the Bhaskara method.
 * - If the polynomial does not meet any of the above criteria, it attempts Briot-Ruffini's method
 *   or performs complex factorization if needed.
 *
 * The function handles polynomials with a varying number of terms and coefficients,
 * including checking for special cases such as unique coefficients and evaluating Eisenstein's criterion.
 *
 * @param p The polynomial to be factored.
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

#endif // POLYNOMIALSZ_H_INCLUDED
