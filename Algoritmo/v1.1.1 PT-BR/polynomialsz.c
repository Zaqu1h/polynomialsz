#include "polynomialsz.h"

/**
 * @file polynomialsz.c
 * @brief Implementação das funcionalidades de fatoração de polinômios.
 *
 * Este arquivo implementa a lógica principal para fatoração de polinômios, incluindo:
 *  - Método de Aberth-Ehrlich (aproximação numérica de raízes);
 *  - Fórmula de Bhaskara (equações quadráticas);
 *  - Briot-Ruffini (divisão sintética e teorema das raízes racionais);
 *  - Fatoração de polinômios ciclotômicos.
 *
 * O algoritmo trabalha com polinômios de coeficientes inteiros.
 *
 * @author Isaque Passos
 * @version 1.1.1
 * @date 2026
 *
 * @note Versão 1.0.0: Implementação inicial.
 * @note Versão 1.1.0: Adicionado método de Aberth e outras melhorias.
 * @note Versão 1.1.1: Correções e melhorias leves.
 */

char var = 'x';
int degreeX = 0;
int sol = 0;
int divider = 0;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int gcd(int a, int b);
int indOfRoot(int rad, double rt);

double nrt(double rad, int ind);

term setTerms(int coef, int exp);
polynomial pCreate(int numTerms);

void pPrint(polynomial p);
void printFac(polynomial p);
void pFree(polynomial p);
void divideX(polynomial p);
void divideGCD(polynomial p);
void removeZeros(polynomial* p);
void fac(polynomial p);

void bhaskara(polynomial p);
void briotRuffini(polynomial p, int maxNum);
void cyclotomicFac(polynomial p);
void aberth(polynomial p);

char* bhaskaraSimplify(int aexp, int den, int b, double delta, int absDelta, double* rootsPair);
char* rationalSqRoots(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
char* irrationalSqRoots(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
char* intSrPositiveDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
char* intSrNegativeDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
char* floatSrPositiveDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
char* floatSrNegativeDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);

//-----------------------------------------------------------------------------

term setTerms(int coef, int exp) {

	term T;
	T.coefficient = coef;
	T.exponent = exp;

	return T;
}

//-----------------------------------------------------------------------------

double nrt(double rad, int ind){

    double base = (pow((fabs)(rad), 1.0/(double)(ind)));

    return base;
}

//-----------------------------------------------------------------------------

int indOfRoot(int rad, double rt){

    int ind = 1;

    while((int)pow(rt, ind) != (fabs)(rad) && ind != rad){
        ind++;
    }

    return ind;
}

//-----------------------------------------------------------------------------

polynomial pCreate(int numTerms) {

    int j = 0;

	polynomial p;
	p.numTerms = numTerms;
	p.terms = (term*)malloc(numTerms * sizeof(term));

	if(p.terms == NULL){

        printf("\nError allocating memory.");
        exit(1);
	}

	for(int i = numTerms-1; i >= 0; i--){

        p.terms[j].coefficient = 0;
        p.terms[j].exponent = i;
        j++;
	}

	return p;
}

//-----------------------------------------------------------------------------

void pPrint(polynomial p){

    printf("%i%c^%i", p.terms[0].coefficient, var, p.terms[0].exponent);

    for(int i = 1; i < p.numTerms; i++){
        printf(" %c %i%c^%i",
            (p.terms[i].coefficient >= 0) ? '+' : '-',
            (abs)(p.terms[i].coefficient), var, p.terms[i].exponent);
    }
}

//-----------------------------------------------------------------------------

int gcd(int a, int b) {

    if(b == 0)
        return a;

    return gcd(b, a % b);
}

//-----------------------------------------------------------------------------

void cyclotomicFac(polynomial p) {

    int grau = p.terms[0].exponent;

    int countRaizes = grau + 1;

    if ((countRaizes % 2 == 0)) {

        printf("(%c + 1)", var);
    }

    for (int k = 1; k <= grau; k++) {

        if (k % countRaizes == 0) continue;

        int conjugadaK = countRaizes - k;

        if (k >= conjugadaK) continue;

        int num = 2 * k;
        int den = countRaizes;
        int divisor = gcd(num, den);

        num /= divisor;
        den /= divisor;

        if (num == 1)
            printf("(%c - Exp[i*Pi/%d])", var, den);
        else
            printf("(%c - Exp[%d*i*Pi/%d])", var, num, den);

        if (num == 1)
            printf("(%c - Exp[-i*Pi/%d])", var, den);
        else
            printf("(%c - Exp[-%d*i*Pi/%d])", var, num, den);
    }
}

//-----------------------------------------------------------------------------

char* intSrPositiveDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot){

    int numeratorSum = (-b + (int)sqrt(absDelta));
    int numeratorSub = (-b - (int)sqrt(absDelta));
    unsigned char divisible = (numeratorSum % den == 0 && numeratorSub % den == 0) ? 1 : 0;

    if(divisible){

            snprintf(bSimplify, 100, "(%c%s %c %i)(%c%s %c %i)", var, powerRoot,
            ((-b + sqrt(absDelta))/den >= 0) ? '-' : '+', (abs)(numeratorSum)/den, var, powerRoot,
            ((-b - sqrt(absDelta))/den >= 0) ? '-' : '+', (abs)(numeratorSub)/den);
        }

        else if(!divisible && b != 0){

            snprintf(bSimplify, 100, "(%c%s - (%i + %i)/%i)(%c%s - (%i - %i)/%i)", var, powerRoot,
             -b, (int)sqrt(absDelta), den, var, powerRoot,
             -b, (int)sqrt(absDelta), den);
        }

        else{

            snprintf(bSimplify, 100, "(%c%s %c (%i/%i))(%c%s %c (%i/%i))", var, powerRoot,
            (-b + sqrt(absDelta) >= 0) ? '-' : '+', (int)sqrt(absDelta), den, var, powerRoot,
            (-b + sqrt(absDelta) >= 0) ? '+' : '-', (int)sqrt(absDelta), den);
        }


    return bSimplify;
}

//-----------------------------------------------------------------------------

char* intSrNegativeDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot){

    if(-b % den == 0 && (int)sqrt(absDelta) % den == 0 && (int)sqrt(absDelta) / den != 1 && b != 0){
        snprintf(bSimplify, 100, "(%c%s - (%i + %ii))(%c%s - (%i - %ii))", var, powerRoot,
                -b/den, (int)sqrt(absDelta)/den, var, powerRoot, -b/den, (int)sqrt(absDelta)/den);
    }
    else if(-b % den == 0 && (int)sqrt(absDelta) % den == 0 && b != 0){
        snprintf(bSimplify, 100, "(%c%s - (%i + i))(%c%s - (%i - i))", var, powerRoot,
                -b/den, var, powerRoot, -b/den);
    }
    else if((-b % den != 0 || (int)sqrt(absDelta) % den != 0) && b != 0){
        snprintf(bSimplify, 100, "(%c%s - (%i + %ii)/%i)(%c%s - (%i - %ii)/%i)", var, powerRoot,
                -b, (int)sqrt(absDelta), den, var, powerRoot, -b, (int)sqrt(absDelta), den);
    }
    else if((int)sqrt(absDelta) % den == 0 && (int)sqrt(absDelta) / den != 1 && b == 0){
        snprintf(bSimplify, 100, "(%c%s - %ii)(%c%s + %ii)", var, powerRoot,
                (int)sqrt(absDelta)/den, var, powerRoot, (int)sqrt(absDelta)/den);
    }
    else if((int)sqrt(absDelta) % den == 0 && b == 0){
        snprintf(bSimplify, 100, "(%c%s - i)(%c%s + i)", var, powerRoot, var, powerRoot);
    }
    else{
        snprintf(bSimplify, 100, "(%c%s - %ii/%i)(%c%s + %ii/%i)", var, powerRoot,
                (int)sqrt(absDelta), den, var, powerRoot, (int)sqrt(absDelta), den);
    }

    return bSimplify;
}

//-----------------------------------------------------------------------------

char* floatSrPositiveDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot){

    if(b != 0){
        snprintf(bSimplify, 100, "(%c%s - ((%i + Sqrt[%i])/%i))(%c%s - ((%i - Sqrt[%i])/%i))", var, powerRoot,
                -b, absDelta, den, var, powerRoot,
                -b, absDelta, den);
    }
    else{
        snprintf(bSimplify, 100, "(%c%s %c (Sqrt[%i]/%i))(%c%s %c (Sqrt[%i]/%i))", var, powerRoot,
               (-b + sqrt(absDelta) >= 0) ? '-' : '+', absDelta, den, var, powerRoot,
               (-b - sqrt(absDelta) >= 0) ? '+' : '-', absDelta, den);
    }

    return bSimplify;
}

//-----------------------------------------------------------------------------

char* floatSrNegativeDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot){

    if(b != 0){
        snprintf(bSimplify, 100, "(%c%s - ((%i + iSqrt[%i])/%i))(%c%s - ((%i - iSqrt[%i])/%i))", var, powerRoot,
                -b, absDelta, den, var, powerRoot,
                -b, absDelta, den);
    }
    else{
        snprintf(bSimplify, 100, "(%c%s %c (iSqrt[%i]/%i))(%c%s %c (iSqrt[%i]/%i))", var, powerRoot,
               (-b + sqrt(absDelta) >= 0) ? '-' : '+', absDelta, den, var, powerRoot,
               (-b - sqrt(absDelta) >= 0) ? '+' : '-', absDelta, den);
    }

    return bSimplify;
}

//-----------------------------------------------------------------------------

char* rationalSqRoots(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot){

    if(delta >= 0)
        return intSrPositiveDelta(aexp, den, b, delta, absDelta, bSimplify, powerRoot);

    return intSrNegativeDelta(aexp, den, b, delta, absDelta, bSimplify, powerRoot);
}

//-----------------------------------------------------------------------------

char* irrationalSqRoots(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot){

    if(delta >= 0){

        return floatSrPositiveDelta(aexp, den, b, delta, absDelta, bSimplify, powerRoot);
    }

    return floatSrNegativeDelta(aexp, den, b, delta, absDelta, bSimplify, powerRoot);
}

//-----------------------------------------------------------------------------

char* bhaskaraSimplify(int aexp, int den, int b, double delta, int absDelta, double* rootsPair){

    char* bSimplify = malloc(sizeof(char) * 100);
    char powerRoot[3] = "";
    unsigned char perfectSquare = (round((int)sqrt(absDelta)) * round((int)sqrt(absDelta)) == absDelta) ? 1 : 0;

    if(aexp == 4){

        snprintf(powerRoot, 3, "^2");
    }

    if(perfectSquare){

        return rationalSqRoots(aexp, den, b, delta, absDelta, bSimplify, powerRoot);
    }
    else{

        return irrationalSqRoots(aexp, den, b, delta, absDelta, bSimplify, powerRoot);
    }
}

//-----------------------------------------------------------------------------

void bhaskara(polynomial p) {

    int aexp = p.terms[0].exponent;

    int a = 0, b = 0, c = 0;

    for (int i = 0; i < p.numTerms; i++) {

        int expo = p.terms[i].exponent;
        int coef = p.terms[i].coefficient;

        if (i == 0)      a = coef;
        else if (i == 1) b = coef;
        else if (i == 2) c = coef;
    }

    int den = 2 * p.terms[0].coefficient;

    double delta = pow(b, 2) - (4 * a * c);
    double rootsPair[2] = {0, 0};

    int absDelta = (int)(fabs)(delta);

    rootsPair[0] = (-b + sqrt(absDelta)) / den;
    rootsPair[1] = (-b - sqrt(absDelta)) / den;

    int asbRoot1 = (int)(fabs)(rootsPair[0]);
    int asbRoot2 = (int)(fabs)(rootsPair[1]);

    char* bSimplify = bhaskaraSimplify(aexp, den, b, delta, absDelta, rootsPair);

    int negRoot = (rootsPair[0] <= 0 || rootsPair[1] <= 0 || (aexp == 4 && delta >= 0)) ? 1 : -1;


    if(negRoot == 1){

        printf("%s", bSimplify);
    }
    else if(aexp == 4){

        printf("(%c %c %i)(%c %c %i)(%c %c %i)(%c %c %i)", var, (-b + sqrt(absDelta) >= 0) ? '-' : '+', asbRoot1,
                                                   var, (-b + sqrt(absDelta) >= 0) ? '+' : '-', asbRoot1,
                                                   var, (-b - sqrt(absDelta) >= 0) ? '-' : '+', asbRoot2,
                                                   var, (-b - sqrt(absDelta) >= 0) ? '+' : '-', asbRoot2);
    }
    else{

        printf("(%c %c %i)(%c %c %i)", var, (-b + sqrt(absDelta) >= 0) ? '-' : '+', asbRoot1,
                                   var, (-b - sqrt(absDelta) >= 0) ? '-' : '+', asbRoot2);
    }

    free(bSimplify);
}

//-----------------------------------------------------------------------------

void briotRuffini(polynomial p, int maxNum){

    int step = 1, i = 0, r = 0, rnum = 1, cont = 1;
    int aexp = p.terms[0].exponent;
    int numTermsATM = p.numTerms;
    double ind = 1, stepf = 1;

    int *root = (int*)calloc(p.numTerms, sizeof(int));

	while (i <= sqrt(maxNum) && aexp != r) {

		step = p.terms[0].coefficient;
        polynomial aux = pCreate(numTermsATM);
        aux.terms[0] = setTerms(step, p.terms[0].exponent - 1);

		for (int k = 1; k < numTermsATM; k++) {

			step = i * step + p.terms[k].coefficient;
			aux.terms[k] = setTerms(step, p.terms[k].exponent - 1);
        }


        if (step == 0) {

            root[r] = i;
            r++;

            for(int j = 0; j < numTermsATM - 1; j++){

                p.terms[j] = aux.terms[j];
            }

            numTermsATM--;
        }
        else i++;

        free(aux.terms);
	}

	if (r < aexp) {

		i = -(int)(sqrt(maxNum));

		while (i != 0 && aexp != r) {

            step = p.terms[0].coefficient;
            polynomial aux = pCreate(numTermsATM);
            aux.terms[0] = setTerms(step, p.terms[0].exponent - 1);

            for (int k = 1; k < numTermsATM; k++) {

                step = i * step + p.terms[k].coefficient;
                aux.terms[k] = setTerms(step, p.terms[k].exponent - 1);
            }

            if (step == 0) {

                root[r] = i;
                r++;

                for(int j = 0; j < numTermsATM - 1; j++){

                    p.terms[j] = aux.terms[j];
                }

                numTermsATM--;
            }
            else i++;

            free(aux.terms);
        }
    }

    p.numTerms = numTermsATM;

    float zeroAux = 0;

	if (r > 0) {

		for (int l = 0; ((root[l] == 0 && l == 0) || root[l] != 0) && l < aexp; l++) {

            if(root[l] == 0){

                zeroAux = 1.5;

                printf("%c", var);
            }
            else{
                printf("(%c %c %i)", var, ((root[l] > 0) ? '-' : '+'), abs(root[l]));
            }
        }

            if(r < aexp){

                removeZeros(&p);

                if (p.terms[0].exponent == 2 && p.numTerms == 3) {

                    if (degreeX > 0) printf("%c^%d", var, degreeX);

                    if (zeroAux != 1.5 && zeroAux != 0) {

                        printf("(%c %c %d)", var, (zeroAux > 0) ? '-' : '+', (int)fabs(zeroAux));
                    }

                    bhaskara(p);
                    sol = 1;
                    return;
                }

                if (p.terms[0].exponent == 4 && p.numTerms == 3 && p.terms[1].exponent == 2) {

                    polynomial quad = pCreate(3);

                    quad.terms[0] = setTerms(p.terms[0].coefficient, 4);
                    quad.terms[1] = setTerms(p.terms[1].coefficient, 2);
                    quad.terms[2] = setTerms(p.terms[2].coefficient, 0);

                    if (degreeX > 0) printf("%c^%d", var, degreeX);

                    if (zeroAux != 1.5 && zeroAux != 0) {

                        printf("(%c %c %d)", var, (zeroAux > 0) ? '-' : '+', (int)fabs(zeroAux));
                    }

                    bhaskara(quad);
                    sol = 1;
                    return;
            }
        }
	}

	if(r == aexp) sol = 1;
	else{
        aberth(p);
        sol = 1;
	}

    free(root);
    return;
}

//-----------------------------------------------------------------------------

void removeZeros(polynomial *p) {

    if (p == NULL || p->terms == NULL || p->numTerms == 0) return;

    int indexWrite = 0;

    for (int indexRead = 0; indexRead < p->numTerms; indexRead++) {
        if (p->terms[indexRead].coefficient != 0) {
            if (indexWrite != indexRead) {

                p->terms[indexWrite] = p->terms[indexRead];
            }

            indexWrite++;
        }
    }

    if (indexWrite == p->numTerms) return;

    if (indexWrite == 0) {

        free(p->terms);
        p->terms = NULL;
        p->numTerms = 0;
        return;
    }

    term *temp = (term*)realloc(p->terms, indexWrite * sizeof(term));

    if (temp != NULL) {

        p->terms = temp;
        p->numTerms = indexWrite;
    } else {

        p->numTerms = indexWrite;
    }
}

//-----------------------------------------------------------------------------

void fac(polynomial p) {

    int qtZeros = 0;
    int maxNum = 0;
    int maxExp = 0;
    int binary = 1;
    int sequence = 1;
    int maxCoefOne = 1;

    sol = 0;

    for (int i = 0; i < p.numTerms; i++) {

        binary *= p.terms[i].exponent;
    }

    if (binary != 0){

        divideX(p);
        return;
    }

    divideGCD(p);

    for (int j = 0; j < p.numTerms; j++) {

        if (j > 0 && p.terms[j-1].exponent != p.terms[j].exponent + 1) {

            maxCoefOne = 0;
            sequence = 0;
        }

        maxCoefOne = ((abs)(p.terms[j].coefficient) * maxCoefOne);

		maxNum = (abs(p.terms[j].coefficient) > abs(maxNum)) ? abs(p.terms[j].coefficient) : abs(maxNum);
		maxExp = (p.terms[j].exponent > maxExp) ? p.terms[j].exponent : maxExp;
    }

    maxCoefOne = (abs)(maxCoefOne);

    for (int i = 0; i < p.numTerms; i++){

        if(p.terms[i].coefficient == 0) qtZeros++;
    }

    if(divider > 1 || divider < -1) printf("%i", divider);

    if(degreeX != 0){

        printf("%c^%d", var, degreeX);

        if(maxExp == 0){

            printf("\n\n");
            return;
        }
    }

    if((p.numTerms == 3) && (p.terms[0].exponent < 5) && (p.terms[0].exponent % 2 == 0) && (p.terms[1].exponent == 1 && p.terms[0].exponent == 2 || p.terms[1].exponent == 2)) {

        bhaskara(p);
        sol = 1;
    }
    if(sequence == 1 && maxCoefOne != 1 && sol == 0){

        briotRuffini(p, maxNum);
    }
    else if(maxCoefOne == 1 && sol == 0){

        cyclotomicFac(p);
        sol = 1;
        degreeX = 0;
        printf("\n\n");
        return;
    }
    if(sol == 0){

        aberth(p);
    }

    degreeX = 0;
    printf("\n\n");

    return;
}

//-----------------------------------------------------------------------------

void divideX(polynomial p) {

    int minExp = p.terms[0].exponent;

    for (int i = 1; i < p.numTerms; i++) {

        if (p.terms[i].exponent < minExp) minExp = p.terms[i].exponent;
    }

    polynomial remainder = pCreate(p.numTerms);

    for (int i = 0; i < p.numTerms; i++) {

        remainder.terms[i] = setTerms(p.terms[i].coefficient, p.terms[i].exponent - minExp);
    }

    removeZeros(&remainder);

    degreeX = minExp;

    fac(remainder);

    free(remainder.terms);
}

//-----------------------------------------------------------------------------

void divideGCD(polynomial p){

    divider = p.terms[0].coefficient;

    for (int i = 1; i < p.numTerms; i++){

        divider = gcd(divider, p.terms[i].coefficient);
    }

    for (int j = 0; j < p.numTerms; j++){

        p.terms[j].coefficient /= divider;
    }
}

//-----------------------------------------------------------------------------

void aberth(polynomial p) {

    int aexp = p.terms[0].exponent;
    int converge = 0;
    double R = 1.0, real = 0.0, imag = 0.0, val = 0.0, angle = 0.0;
    double complex pVal = 0.0, pDer = 0.0, sum = 0.0, q = 0.0, adjustment = 0.0, newRoot = 0.0;
    double *coef = NULL;
    double complex *roots = NULL;

    coef = (double*)calloc(aexp + 1, sizeof(double));

    for (int i = 0; i < p.numTerms; i++) {

        int exp = p.terms[i].exponent;

        if (exp >= 0 && exp <= aexp) coef[aexp - exp] = p.terms[i].coefficient;
    }

    for (int i = 1; i <= aexp; i++) {

        val = fabs(coef[i] / coef[0]);
        if (val > R) R = val;
    }

    R += 1.0;

    roots = malloc(aexp * sizeof(double complex));

    for (int i = 0; i < aexp; i++) {

        angle = 2.0 * M_PI * i / aexp;
        roots[i] = R * (cos(angle) + I * sin(angle));
    }

    for (int iter = 0; iter < 30; iter++) {

        converge = 1;

        for (int i = 0; i < aexp; i++) {

            pVal = coef[0];
            pDer = coef[0] * aexp;

            for (int j = 1; j < aexp; j++) {

                pVal = pVal * roots[i] + coef[j];
                pDer = pDer * roots[i] + coef[j] * (aexp - j);
            }

            pVal = pVal * roots[i] + coef[aexp];

            sum = 0.0;

            for (int j = 0; j < aexp; j++) {

                if (j != i) sum += 1.0 / (roots[i] - roots[j]);
            }

            q = pVal / pDer;

            adjustment = q / (1.0 - q * sum);

            newRoot = roots[i] - adjustment;

            if (cabs(newRoot - roots[i]) > 1e-12) converge = 0;

            roots[i] = newRoot;
        }

        if (converge) break;
    }

    printf("(");

    for (int i = 0; i < aexp; i++) {

        real = creal(roots[i]);
        imag = cimag(roots[i]);

        if (fabs(real) < 1e-6) real = 0.0;
        if (fabs(imag) < 1e-6) imag = 0.0;

        if (imag == 0.0) printf("(%c %c %.6f)", var, (real >= 0) ? '-' : '+', fabs(real));
        else if (real == 0.0) printf("(%c %c %.6fi)", var, (imag >= 0) ? '-' : '+', fabs(imag));
        else printf("(%c - (%.6f %c %.6fi))", var, real, (imag >= 0) ? '+' : '-', fabs(imag));
    }

    printf(")");

    free(coef);
    free(roots);
}

//-----------------------------------------------------------------------------

void printFac(polynomial p){

    pPrint(p);
    printf(" = ");

    if(p.numTerms == 1 && p.terms[0].exponent == 0){

        printf("%i\n\n", p.terms[0].coefficient);
        return;
    }

    fac(p);
    degreeX = 0;
}

//-----------------------------------------------------------------------------

void pFree(polynomial p) {

    printFac(p);
    free(p.terms);
    p.terms = NULL;
}
