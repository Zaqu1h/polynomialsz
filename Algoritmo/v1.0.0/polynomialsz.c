#include "polynomialsz.h"

/**
 * @file polynomialsz.c
 * @brief Source file for polynomial factorization functionalities.
 *
 * This file implements the core logic for factoring polynomials, including
 * algorithms such as Briot-Ruffini, Eisenstein's Criterion, and power difference
 * factoring. It uses the data structures and declarations defined in the
 * corresponding header file (polynomialsz.h).
 */

/**
 * @brief The symbolic variable used to represent roots and polynomial expressions.
 *
 * This variable is initialized with 'x' by default, but can be changed as needed.
 */
char var = 'x';

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

    while((int)pow(rt, ind) != (fabs)(rad)){
        ind++;
    }

    return ind;
}

//-----------------------------------------------------------------------------

polynomial pCreate(int num_of_terms) {

    int j = 0;
	polynomial p;
	p.num_of_terms = num_of_terms;
	p.terms = (term*)malloc(num_of_terms * sizeof(term));

	if(p.terms == NULL){

        printf("\nError allocating memory.");
        return p;
	}

	for(int i = num_of_terms-1; i >= 0; i--){

        p.terms[j].coefficient = 0;
        p.terms[j].exponent = i;
        j++;
	}

	return p;
}

//-----------------------------------------------------------------------------

void pPrint(polynomial p){

    printf("(%i%c^%i", p.terms[0].coefficient, var, p.terms[0].exponent);

    for(int i = 1; i < p.num_of_terms; i++){

        printf("%c%i%c^%i", (p.terms[i].coefficient >= 0) ? '+' : '-', (abs)(p.terms[i].coefficient), var, p.terms[i].exponent);
    }

    printf(")\n\n");
}

//-----------------------------------------------------------------------------

int gcd(int a, int b) {

    if(b == 0)
        return a;

    return gcd(b, a % b);
}
//-----------------------------------------------------------------------------

int isPrime(int num){

    if(num > 41)
        return 0;

    int i = 2;
    int pcont = 0;
    int dcont = 0;

    for(i = 2; i < 180 && pcont != num; i++){

        for(int k = 1; k <= i; k++){

            dcont = (i % k == 0) ? dcont + 1 : dcont;
        }
        pcont = (dcont == 2) ? pcont + 1 : pcont;
        dcont = 0;
    }

    return i - 1;
}

//-----------------------------------------------------------------------------

int eisenstein(polynomial x, double p){

    int irr = 1;

    for(int i = 1; i < x.num_of_terms; i++){

        irr = (x.terms[i].coefficient % (int)p == 0) ? irr * 1 : 0;
    }

    if(x.terms[x.num_of_terms - 1].coefficient % (int)pow(p, 2) != 0)
        irr *= 1;
    else
        irr = 0;


    if(irr == 1)
        return 1;
    else
        return -1;
}

//-----------------------------------------------------------------------------

void complexp(polynomial p) {

    double angle = 0;
    int numerator = 0;
    int den = 1;

    for(int l = 0; l < p.num_of_terms; l++){

        den *= p.terms[l].coefficient;
    }

    int degree = (den == 1) ? p.num_of_terms : p.terms[0].exponent;

    int total_roots = (degree % 2 == 0) ? degree - 1 : degree;

    if(den != 1) total_roots++;

    int common_divisor = 0;

    if (degree % 2 == 0) {
        printf("(%c + 1)", var);
    }

    if (den == 1){

        for (int k = 1; k <= total_roots / 2; k++) {

            den = degree;
            angle = 2 * M_PI * k / degree;
            numerator = 2 * k;
            common_divisor = gcd(numerator, den);

            numerator /= common_divisor;
            den /= common_divisor;

            char angle_str[50];

            if (den == 1){
                snprintf(angle_str, 50, "%i*pi", numerator);
            }
            else if (numerator == 1){
                snprintf(angle_str, 50, "*pi/%i", den);
            }
            else{
                snprintf(angle_str, 50, "%i*pi/%i", numerator, den);
            }

            if(den == 2 && numerator == 1){

                printf("(%c - i)", var);
                printf("(%c + i)", var);
                return;
            }

            printf("(%c - e^{i%s})", var, angle_str);
            printf("(%c - e^{-i%s})", var, angle_str);
        }
    }
    else{

        for (int k = 0; k < p.terms[0].exponent; k++) {

            den = degree;
            numerator = (2 * k) + 1;
            common_divisor = gcd(numerator, den);

            numerator /= common_divisor;
            den /= common_divisor;

            char angle_str[50];

            if (den == 1){
                snprintf(angle_str, 50, "%i*pi", numerator);
            }
            else if (numerator == 1){
                snprintf(angle_str, 50, "*pi/%i", den);
            }
            else{
                snprintf(angle_str, 50, "%i*pi/%i", numerator, den);
            }

            if(numerator != den) printf("(%c - e^{i%s})", var, angle_str);
            else printf("(%c + 1)", var);
        }
    }
}

//-----------------------------------------------------------------------------

char* bhaskaraSimplify(int aexp, int den, int b, double delta, double* rootsPair){

    int absDelta = (int)(fabs)(delta);
    char* bSimplify = malloc(sizeof(char) * 60);
    char powerRoot[3] = "";

    if(aexp == 4){

        snprintf(powerRoot, 3, "^2");
    }

    if((int)sqrt(absDelta) * (int)sqrt(absDelta) == absDelta){

        if(delta >= 0){

            if((-b + (int)sqrt(absDelta)) % den == 0){

                snprintf(bSimplify, 60, "(%c%s %c %i)(%c%s %c %i)", var, powerRoot,
               (-b + absDelta >= 0) ? '-' : '+', (-b + (int)sqrt(absDelta))/den, var, powerRoot,
               (-b + absDelta >= 0) ? '+' : '-', (-b + (int)sqrt(absDelta))/den);
            }
            else if((-b + (int)sqrt(absDelta)) % den != 0 && b != 0){

                snprintf(bSimplify, 60, "(%c%s %c (%i + %i)/%i)(%c%s %c (%i - %i)/%i)", var, powerRoot,
               (-b + absDelta >= 0) ? '-' : '+', -b, (int)sqrt(absDelta), den, var, powerRoot,
               (-b + absDelta >= 0) ? '+' : '-', -b, (int)sqrt(absDelta), den);
            }
            else{

                snprintf(bSimplify, 60, "(%c%s %c (%i/%i))(%c%s %c (%i/i))", var, powerRoot,
               (-b + absDelta >= 0) ? '-' : '+', (int)sqrt(absDelta), den, var, powerRoot,
               (-b + absDelta >= 0) ? '+' : '-', (int)sqrt(absDelta), den);
            }
        }
        else{
            if(-b % den == 0 && (int)sqrt(absDelta) % den == 0 && (int)sqrt(absDelta) / den != 1 && b != 0){

                snprintf(bSimplify, 60, "(%c%s %c (%i + i%i))(%c%s %c (%i - i%i))", var, powerRoot,
               (-b + absDelta >= 0) ? '-' : '+', -b/den, (int)sqrt(absDelta)/den, var, powerRoot,
               (-b + absDelta >= 0) ? '+' : '-', -b/den, (int)sqrt(absDelta)/den);
            }
            else if(-b % den == 0 && (int)sqrt(absDelta) % den == 0 && b != 0){

                snprintf(bSimplify, 60, "(%c%s %c (%i + i))(%c%s %c (%i - i))", var, powerRoot,
               (-b + absDelta >= 0) ? '-' : '+', -b/den, var, powerRoot,
               (-b + absDelta >= 0) ? '+' : '-', -b/den);
            }

            else if(-b % den != 0 || (int)sqrt(absDelta) % den != 0 && b != 0){

                snprintf(bSimplify, 60, "(%c%s %c (%i + i%i)/%i)(%c%s %c (%i - i%i)/%i)", var, powerRoot,
               (-b + absDelta >= 0) ? '-' : '+', -b, (int)sqrt(absDelta), den, var, powerRoot,
               (-b + absDelta >= 0) ? '+' : '-', -b, (int)sqrt(absDelta), den);
            }
            else if((int)sqrt(absDelta) % den == 0 && (int)sqrt(absDelta) / den != 1 && b == 0){

                snprintf(bSimplify, 60, "(%c%s %c i%i)(%c%s %c i%i)", var, powerRoot,
               (-b + absDelta >= 0) ? '-' : '+', (int)sqrt(absDelta)/den, var, powerRoot,
               (-b + absDelta >= 0) ? '+' : '-', (int)sqrt(absDelta)/den);
            }
            else if((int)sqrt(absDelta) % den == 0 && b == 0){

                snprintf(bSimplify, 60, "(%c%s %c i)(%c%s %c i)", var, powerRoot,
               (-b + absDelta >= 0) ? '-' : '+', var, powerRoot,
               (-b + absDelta >= 0) ? '+' : '-');
            }
            else{

                snprintf(bSimplify, 60, "(%c%s %c i%i/%i)(%c%s %c i%i/%i)", var, powerRoot,
               (-b + absDelta >= 0) ? '-' : '+', (int)sqrt(absDelta), den, var, powerRoot,
               (-b + absDelta >= 0) ? '+' : '-', (int)sqrt(absDelta), den);
            }
        }
    }
    else{
        if(delta >= 0){

            if(b != 0){

                snprintf(bSimplify, 60, "(%c%s %c ((%i + %i^{1/2})/%i)) (%c%s %c ((%i - %i^{1/2})/%i))", var, powerRoot,
               (rootsPair[0] >= 0) ? '-' : '+', -b, absDelta, den, var, powerRoot,
               (rootsPair[1] >= 0) ? '+' : '-', -b, absDelta, den);
            }
            else{

                snprintf(bSimplify, 60, "(%c%s %c (%i^{1/2})/%i) (%c%s %c (%i^{1/2})/%i)", var, powerRoot,
               (rootsPair[0] >= 0) ? '-' : '+', absDelta, den, var, powerRoot,
               (rootsPair[1] >= 0) ? '+' : '-', absDelta, den);
            }
        }
        else{

            if(b != 0){

                snprintf(bSimplify, 60, "(%c%s %c ((%i + i%i^{1/2})/%i)) (%c%s %c ((%i - i%i^{1/2})/%i))", var, powerRoot,
               (rootsPair[0] >= 0) ? '-' : '+', -b, absDelta, den, var, powerRoot,
               (rootsPair[1] >= 0) ? '+' : '-', -b, absDelta, den);
            }
            else{

                snprintf(bSimplify, 60, "(%c%s %c (i%i^{1/2})/%i) (%c%s %c (i%i^{1/2})/%i)", var, powerRoot,
               (rootsPair[0] >= 0) ? '-' : '+', absDelta, den, var, powerRoot,
               (rootsPair[1] >= 0) ? '+' : '-', absDelta, den);
            }
        }
    }

    return bSimplify;
}

//-----------------------------------------------------------------------------

void bhaskara(polynomial p) {

    int aexp = p.terms[0].exponent;
    int a = p.terms[0].coefficient;
    int b = p.terms[1].coefficient;
    int c = p.terms[2].coefficient;

    int den = 2 * a;
    double delta = pow(b, 2) - (4 * a * c);
    double rootsPair[2] = {0, 0};

    rootsPair[0] = (-b + (sqrt((fabs)(delta)))) / den;
    rootsPair[1] = (-b - (sqrt((fabs)(delta)))) / den;

    int asbRoot1 = (int)(fabs)(rootsPair[0]);
    int asbRoot2 = (int)(fabs)(rootsPair[1]);

    char* bSimplify = bhaskaraSimplify(aexp, den, b, delta, rootsPair);

    int negRoot = (rootsPair[0] <= 0 || rootsPair[1] <= 0 || (aexp == 4 && delta >= 0)) ? 1 : -1;

    if(negRoot == 1){

        printf("%s", bSimplify);
    }
    else if(aexp == 4){

        printf("(%c%c%i)(%c%c%i)(%c%c%i)(%c%c%i)", var, (rootsPair[0] >= 0) ? '-' : '+', asbRoot1,
                                                   var, (rootsPair[0] >= 0) ? '+' : '-', asbRoot1,
                                                   var, (rootsPair[1] >= 0) ? '-' : '+', asbRoot2,
                                                   var, (rootsPair[1] >= 0) ? '+' : '-', asbRoot2);
    }
    else{

        printf("(%c%c%i)(%c%c%i)", var, (rootsPair[0] >= 0) ? '-' : '+', asbRoot1,
                                   var, (rootsPair[1] >= 0) ? '-' : '+', asbRoot2);
    }

    free(bSimplify);
}

//-----------------------------------------------------------------------------

void powerDiffFactoring(polynomial p){

    int i = p.terms[0].exponent - 1, even_num = 1;
    int count = 1;
    double a = (nrt(p.terms[1].coefficient, p.terms[0].exponent));

    printf("(%c%c%i)", var, (p.terms[1].coefficient < 0) ? '-' : '+', (int)a);

    if(p.terms[0].exponent % 2 == 0){
        printf("(%c%c%i)", var, (p.terms[1].coefficient < 0) ? '+' : '-', (int)a);
        even_num++;
    }

    if(p.terms[0].exponent < 2){
        return;
    }

    printf("(%c^{%i}", var, (even_num == 2) ? --i : i);

    while(i - even_num > 0){
        printf(" %c %i%c^{%i}",
               (p.terms[1].coefficient > 0 && i % 2 == 0) ? '-' : '+',
               (even_num == 2) ? (int)(pow(a, ++count)) : (int)(pow(a, count++)),
               var, i - even_num);

        if(even_num == 2){
            i = i - even_num;
            count++;
        }
        else i--;
    }

    printf(" + %i)", (int)(pow(a, p.terms[0].exponent - even_num)));

    return;
}

//-----------------------------------------------------------------------------

void briot_ruffini(polynomial p){

    int max_num = 0, start = 1, i = 1, r = 1, rq = 0, end_int = 0, rnum = 1, rI = 1, cont = 1;
    double ind = 1, startf = 1;

    int *root = (int*)calloc(p.num_of_terms, sizeof(int));

	for (int j = 1; j < p.num_of_terms; j++) {
		max_num = (p.terms[j].coefficient > max_num) ? p.terms[j].coefficient : max_num;
	}

	while (i != max_num) {

		start = i * 1 + p.terms[1].coefficient;

		for (int k = 2; k < p.num_of_terms; k++) {
			start = i * start + p.terms[k].coefficient;
		}

		if (start == 0) {
			root[r] = i;
			r++;
		}

		i++;
	}

	if (r < p.num_of_terms) {
		i = -max_num;

		while (i != 0) {
			start = i * 1 + p.terms[1].coefficient;

			for (int k = 2; k < p.num_of_terms; k++) {
				start = i * start + p.terms[k].coefficient;
			}

			if (start == 0) {
				root[r] = i;
				r++;
			}

			i++;
		}
	}

	end_int = r;

	if (r < p.num_of_terms) {
		rq = 2;

		while (rq < max_num && r < p.num_of_terms) {

			ind = sqrt(rq);
			startf = ind + p.terms[1].coefficient;

			for (int k = 2; k <= p.num_of_terms; k++) {
				startf = ind * startf + p.terms[k].coefficient;

				if ((int)(startf * 1000) == 0) {
					root[r] = rq;
					root[r + 1] = rq;
					r += 2;
					break;
				}
			}

			rq++;

			for (int rr = 0; rr < end_int; rr++) {
				if(root[rr] == sqrt(rq)) {
					rq++;
				}
			}
		}
	}

	int pr = 0;

	if (r != p.num_of_terms-1) {
		for (int l = 1; root[l] != 0 && l < p.num_of_terms; l++) {

			if (l < end_int) {
				printf("(%c%c%d)", var, ((root[l] > 0) ? '-' : '+'), abs(root[l]));
				pr++;
			}
			else if (l >= end_int && p.terms[0].exponent != 3) {
				printf("(%c + %d^{1/2})(%c - %d^{1/2})", var, root[l], var, root[l]);
				l++;
				pr++;
			}
			if (rq > 0 && p.terms[0].exponent % 2 != 0 && p.terms[0].exponent < 6 && (p.terms[p.num_of_terms - 1].coefficient >= 0)){

                int px = (p.terms[0].exponent == 3) ? p.terms[1].coefficient + root[1] : p.terms[3].coefficient / -root[1];
                int qx = (p.terms[0].exponent == 3) ? p.terms[2].coefficient + px*root[1] : p.terms[5].coefficient / -root[1];

                polynomial aux = pCreate(3);
                aux.terms[0] = setTerms(1, p.terms[0].exponent - 1);
                aux.terms[1] = setTerms(px, p.terms[0].exponent / 2);
                aux.terms[2] = setTerms(qx, 0);

                bhaskara(aux);
                printf("\n\n");
                return;
			}
		}
	}

    if(pr == 0){

        printf("There is no non-imaginary factorization and/or no root representable by radicals or integers for this polynomial.\n\n");
        free(root);
        return;
    }

    printf(";\n\n");
    free(root);
    return;
}

//-----------------------------------------------------------------------------

void printFac(polynomial p) {

    int check_coefficient = (p.num_of_terms > 2 && p.terms[p.num_of_terms - 1].coefficient == 1) ? 1 : 0;

    for (int j = 0; j < p.num_of_terms; j++) {
        check_coefficient = ((p.terms[j].coefficient) == -1) ? (check_coefficient + 2) : check_coefficient;
        check_coefficient = ((p.terms[j].coefficient) == 1) ? (check_coefficient * check_coefficient) : check_coefficient;
    }

    int irrdcble = 0;

    for (int i = 1; i < 42 && irrdcble != 1; i++) {
        irrdcble = eisenstein(p, isPrime(i));
    }

    /*
     * This condition checks if the polynomial meets one of two criteria:
     *
     * 1. The polynomial is irreducible, has 5 or more terms, and satisfies one of the following:
     *    - The highest degree is even, and the constant term (last coefficient) is positive.
     *    - The highest degree is odd, and the constant term is negative.
     *
     * 2. The polynomial is not irreducible (irrdcble == -1), has 5 or more terms, and its coefficient check value (check_coefficient) is greater than 1.
     *    This implies the polynomial may not follow a simple factorization or progression pattern.
     *
     * If either of these conditions is met, the function prints a message indicating that
     * no non-imaginary factorization and/or no root representable by radicals or integers
     * is possible for the polynomial, and the function returns.
     */
    if (irrdcble == 1 && p.num_of_terms >= 5 && ((p.num_of_terms % 2 == 0 && p.terms[p.num_of_terms - 1].coefficient > 0)
        || p.num_of_terms % 2 != 0 && p.terms[p.num_of_terms - 1].coefficient < 0) || irrdcble == -1 && p.num_of_terms >= 5 && check_coefficient > 1){

        printf("There is no non-imaginary factorization and/or no root representable by radicals or integers for this polynomial.\n\n");

        return;
    } else {
        // Handle case for a quadratic polynomial
        if (p.num_of_terms == 2 && p.terms[0].exponent == indOfRoot(p.terms[1].coefficient, (nrt(p.terms[1].coefficient, p.terms[0].exponent)))) {

            if ((p.terms[1].coefficient < 0 && p.terms[0].exponent % 2 == 0) || p.terms[0].exponent % 2 != 0) {

                // Apply power difference factoring
                powerDiffFactoring(p);

                // Apply complex factorization if needed
                if (check_coefficient == 1) complexp(p);

                printf(";\n\n");
                return;
            }
        }

        // Apply Bhaskara's method if the polynomial has an even degree and 3 or fewer terms
        if (p.terms[0].exponent % 2 == 0 && p.num_of_terms < 5 && p.terms[1].exponent != 3) {

            bhaskara(p);

            printf(";\n\n");
            return;
        }

        // Apply Briot-Ruffini's method for non-unique coefficients
        if (check_coefficient != 1) briot_ruffini(p);

        // Apply complex factorization if needed
        if (check_coefficient == 1){

            complexp(p);
            printf(";\n\n");
        }

        return;
    }
}

//-----------------------------------------------------------------------------

void pFree(polynomial p) {

    int check_coefficient = (p.num_of_terms > 2 && p.terms[p.num_of_terms - 1].coefficient == 1) ? 1 : 0;

    for (int j = 0; j < p.num_of_terms; j++) {
        check_coefficient = ((p.terms[j].coefficient) == -1) ? (check_coefficient + 2) : check_coefficient;
        check_coefficient = ((p.terms[j].coefficient) == 1) ? (check_coefficient * check_coefficient) : check_coefficient;
    }

    int irrdcble = 0;

    for (int i = 1; i < 42 && irrdcble != 1; i++) {
        irrdcble = eisenstein(p, isPrime(i));
    }

    if (irrdcble == 1 && p.num_of_terms >= 5 && ((p.num_of_terms % 2 == 0 && p.terms[p.num_of_terms - 1].coefficient > 0)
        || p.num_of_terms % 2 != 0 && p.terms[p.num_of_terms - 1].coefficient < 0) || irrdcble == -1 && p.num_of_terms >= 5 && check_coefficient > 1){

        printf("There is no non-imaginary factorization and/or no root representable by radicals or integers for this polynomial.\n\n");
        free(p.terms);

        return;
    } else {

        if (p.num_of_terms == 2 && p.terms[0].exponent == indOfRoot(p.terms[1].coefficient, (nrt(p.terms[1].coefficient, p.terms[0].exponent)))) {

            if ((p.terms[1].coefficient < 0 && p.terms[0].exponent % 2 == 0) || p.terms[0].exponent % 2 != 0) {

                powerDiffFactoring(p);

                if (check_coefficient == 1) complexp(p);

                free(p.terms);

                printf(";\n\n");
                return;
            }
        }

        if (p.terms[0].exponent % 2 == 0 && p.num_of_terms < 5 && p.terms[1].exponent != 3) {

            bhaskara(p);

            free(p.terms);

            printf(";\n\n");
            return;
        }

        if (check_coefficient != 1) briot_ruffini(p);

        if (check_coefficient == 1){

            complexp(p);
            printf(";\n\n");
        }

        free(p.terms);

        return;
    }

    free(p.terms);
}
