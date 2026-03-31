#include "polynomialsz.h"

int main() {

	polynomial U = pCreate(4);
	U.terms[0] = setTerms(1, 3);
	U.terms[1] = setTerms(-6, 2);
	U.terms[2] = setTerms(11, 1);
	U.terms[3] = setTerms(-6, 0);

    polynomial E = pCreate(5);
    E.terms[0] = setTerms(1, 4);
    E.terms[1] = setTerms(2, 3);
    E.terms[2] = setTerms(-7, 2);
    E.terms[3] = setTerms(-8, 1);
    E.terms[4] = setTerms(12, 0);

    polynomial S = pCreate(6);
    S.terms[0] = setTerms(1, 5);
    S.terms[1] = setTerms(-7, 4);
    S.terms[2] = setTerms(-8, 3);
    S.terms[3] = setTerms(56, 2);
    S.terms[4] = setTerms(15, 1);
    S.terms[5] = setTerms(-105, 0);

	polynomial C = pCreate(7);
	C.terms[0] = setTerms(1, 6);
	C.terms[1] = setTerms(-5, 5);
	C.terms[2] = setTerms(-6, 4);
	C.terms[3] = setTerms(60, 3);
	C.terms[4] = setTerms(-37, 2);
	C.terms[5] = setTerms(-175, 1);
	C.terms[6] = setTerms(210, 0);

	polynomial X = pCreate(3);
	X.terms[0] = setTerms(4, 4);
	X.terms[1] = setTerms(1, 2);
	X.terms[2] = setTerms(1, 0);

	polynomial P1 = pCreate(4);
	P1.terms[0] = setTerms(1, 3);
	P1.terms[1] = setTerms(0, 2);
	P1.terms[2] = setTerms(-2, 1);
	P1.terms[3] = setTerms(4, 0);

    polynomial P2 = pCreate(5);
	P2.terms[0] = setTerms(1, 4);
	P2.terms[1] = setTerms(1, 3);
	P2.terms[2] = setTerms(1, 2);
	P2.terms[3] = setTerms(1, 1);
	P2.terms[4] = setTerms(1, 0);

    polynomial P3 = pCreate(6);
	P3.terms[0] = setTerms(1, 5);
	P3.terms[1] = setTerms(1, 4);
	P3.terms[2] = setTerms(1, 3);
	P3.terms[3] = setTerms(1, 2);
	P3.terms[4] = setTerms(1, 1);
	P3.terms[5] = setTerms(1, 0);

    polynomial P4 = pCreate(7);
	P4.terms[0] = setTerms(1, 6);
	P4.terms[1] = setTerms(1, 5);
	P4.terms[2] = setTerms(1, 4);
	P4.terms[3] = setTerms(1, 3);
	P4.terms[4] = setTerms(1, 2);
    P4.terms[5] = setTerms(1, 1);
	P4.terms[6] = setTerms(1, 0);

    polynomial Y = pCreate(2);
	Y.terms[0] = setTerms(1, 5);
	Y.terms[1] = setTerms(-32, 0);

	polynomial Fib = pCreate(3);
	Fib.terms[0] = setTerms(1, 2);
	Fib.terms[1] = setTerms(1, 1);
	Fib.terms[2] = setTerms(-1, 0);

	polynomial D = pCreate(17);
    D.terms[0] = setTerms(1, 16);
    D.terms[1] = setTerms(1, 15);
    D.terms[2] = setTerms(1, 14);
    D.terms[3] = setTerms(1, 13);
    D.terms[4] = setTerms(1, 12);
    D.terms[5] = setTerms(1, 11);
    D.terms[6] = setTerms(1, 10);
    D.terms[7] = setTerms(1, 9);
    D.terms[8] = setTerms(1, 8);
    D.terms[9] = setTerms(1, 7);
    D.terms[10] = setTerms(1, 6);
    D.terms[11] = setTerms(1, 5);
    D.terms[12] = setTerms(1, 4);
    D.terms[13] = setTerms(1, 3);
    D.terms[14] = setTerms(1, 2);
    D.terms[15] = setTerms(1, 1);
    D.terms[16] = setTerms(1, 0);

    polynomial T = pCreate(3);
	T.terms[0] = setTerms(1, 5);
	T.terms[1] = setTerms(-1, 1);
	T.terms[2] = setTerms(1, 0);

    polynomial F = pCreate(2);
	F.terms[0] = setTerms(1, 5);
	F.terms[1] = setTerms(1, 0);

	pFree(U);
	pFree(E);
	pFree(S);
	pFree(C);
	pFree(X);
	pFree(P1);
	pFree(P2);
	pFree(P3);
	pFree(P4);
    pFree(Y);
	pFree(Fib);
	pFree(D);
    pFree(T);
	pFree(F);

	return 0;
}
