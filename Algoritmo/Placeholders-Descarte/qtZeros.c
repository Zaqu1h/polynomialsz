int qtZeros = 0;

    for(int i = 0; i < p.num_of_terms; i++){

        qtZeros = (p.terms[i].coefficient == 0) ? qtZeros + 1 : qtZeros;
    }

    if(qtZeros != 0 && p.terms[0].exponent % 2 == 0){

        polynomial aux = pCreate(p.num_of_terms - qtZeros);
        int coefAux = p.num_of_terms - 1;
        int expAux = 0;

        for(int j = p.num_of_terms - qtZeros - 1; j >= 0; j--){

            expAux = (p.terms[coefAux].coefficient == 0) ? expAux + 1 : expAux;
            coefAux = (p.terms[coefAux].coefficient == 0) ? coefAux - 1 : coefAux;

            aux.terms[j] = setTerms(p.terms[coefAux--].coefficient, expAux++);
        }

        printFac(aux);
        return;
    }
