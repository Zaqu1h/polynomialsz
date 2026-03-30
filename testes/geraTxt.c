#include <stdio.h>
#include <stdlib.h>
#include "polynomialsz.h"

#define MAX_COEF 3
#define MAX_DEGREE 5

int main() {

    char filename[50];

    sprintf(filename, "teste%dg%d.txt", MAX_COEF, MAX_DEGREE);

    FILE* file = fopen(filename, "w");
    if (!file) return 1;

    int base = MAX_COEF + 1;
    int total = 1;

    for (int i = 0; i <= MAX_DEGREE; i++) {

        total *= base;
    }

    int arrayCoefs[MAX_DEGREE + 1];

    for (int n = 0; n < total; n++) {

        int temp = n;

        for (int i = 0; i <= MAX_DEGREE; i++) {
            arrayCoefs[i] = temp % base;
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
