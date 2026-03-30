#include <stdio.h>
#include <stdlib.h>
#include "polynomialsz.h"

#define MAX_DEGREE 5

int main() {

    char filename[50];

    sprintf(filename, "primesg%d.txt", MAX_DEGREE);

    FILE* file = fopen(filename, "w");

    if (!file) return 1;

    int coefs[] = {0, 2, 3, 5, 7};

    int base = sizeof(coefs) / sizeof(coefs[0]);

    int total = 1;

    for (int i = 0; i <= MAX_DEGREE; i++) {

        total *= base;
    }

    int arrayCoefs[MAX_DEGREE + 1];

    for (int n = 0; n < total; n++) {

        int temp = n;

        for (int i = 0; i <= MAX_DEGREE; i++) {

            arrayCoefs[i] = coefs[temp % base];
            temp /= base;
        }

        int numTerms = 0;

        for (int i = 0; i <= MAX_DEGREE; i++) {

            if (arrayCoefs[i] != 0) numTerms++;
        }

        if (numTerms == 0) continue;

        polynomial p = pCreate(numTerms);

        int index = 0;

        for (int exp = MAX_DEGREE; exp >= 0; exp--) {

            if (arrayCoefs[exp] != 0) {

                p.terms[index] = setTerms(arrayCoefs[exp], exp);
                index++;
            }
        }

        p.numTerms = numTerms;

        FILE* original_stdout = stdout;
        stdout = file;
        pFree(p);
        stdout = original_stdout;
    }

    fclose(file);
    return 0;
}
