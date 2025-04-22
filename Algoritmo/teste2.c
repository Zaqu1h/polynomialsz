#include "polynomialsz.h"

int main() {

	polynomial U = pCreate(3);
	U.terms[0] = setTerms(1, 2);
    U.terms[1] = setTerms(0, 1);
	U.terms[2] = setTerms(4, 0);

	//r = raiz = -2

	//a = 0
	//b = -2
	//c = 4

	//p = a + raiz = 0 - 2 = -2
	//q = b + pr = -2 + (-2 * -2) = 2

	//x^3 + ax^2 + bx + c = x^3 + 0x^2 + -2x + 4
	//x^2 + px + q = x^2 - 2x + 2

    printFac(U);

    polynomial E = pCreate(6);
	E.terms[0] = setTerms(1, 5);
	E.terms[1] = setTerms(2, 4);
	E.terms[2] = setTerms(4, 3);
	E.terms[3] = setTerms(8, 2);
	E.terms[4] = setTerms(16, 1);
	E.terms[5] = setTerms(32, 0);


    printFac(E);

	return 0;
}
