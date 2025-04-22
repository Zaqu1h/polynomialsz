#ifndef POLYNOMIALSZ_H_INCLUDED
#define POLYNOMIALSZ_H_INCLUDED

/**
 * @file polynomialsz.h
 * @brief Arquivo de cabeçalho para estruturas de fatoração polinomial.
 *
 * Este arquivo fornece as estruturas de dados necessárias e declarações usadas pelo
 * algoritmo de fatoração de polinômios, que suporta:
 *  - Critério de Eisenstein
 *  - Briot-Ruffini (divisão sintética)
 *  - Método de Bhaskara (fórmula quadrática)
 *  - Fatoração por diferença de potências
 *  - Análise de raízes complexas
 *
 * O algoritmo trabalha com polinômios com coeficientes inteiros.
 *
 * @author Isaque Passos
 * @version 1.0.0
 * Lançamento público inicial. Versões internas e testes existiram antes disso.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> /**< Necessário para operações matemáticas. */
#include <locale.h> /**< Para imprimir a mensagem de erro em Pt-BR adequadamente. */

#define _USE_MATH_DEFINES /**< Habilita o uso de M_PI e outras constantes em alguns compiladores. */

/**
 * @brief Variável simbólica global usada nas expressões polinomiais.
 *
 * Essa variável é usada para exibir as raízes polinomiais simbolicamente.
 * Por exemplo, se var = 'x', as raízes de x² - 2 são impressas como x + √2 e x - √2.
 */
extern char var;

/**
 * @struct term
 * @brief Representa um único termo de um polinômio.
 *
 * Um termo é composto por um coeficiente e um expoente.
 * Por exemplo, o termo 3x² tem coeficiente = 3 e expoente = 2.
 */
typedef struct sterm {
    int coefficient; /**< O coeficiente do termo. */
    int exponent;    /**< O expoente do termo. */
} term;

/**
 * @struct polynomial
 * @brief Representa um polinômio como um array de termos.
 *
 * Os termos não precisam ser ordenados, mas geralmente são ordenados por expoentes decrescentes.
 */
typedef struct spolynomial {
    term *terms;        /**< Array dinâmico de termos do polinômio. */
    int num_of_terms;   /**< Número de termos no polinômio. */
} polynomial;
//-----------------------------------------------------------------------------
/**
 * @brief Cria e define o coeficiente e o expoente de um termo polinomial.
 *
 * @param coef O coeficiente do termo.
 * @param exp O expoente do termo.
 * @return T O termo completo.
 */
term setTerms(int coef, int exp);
//-----------------------------------------------------------------------------
/**
 * @brief Retorna a base de uma potência com expoente igual a `ind`, que é igual a `rad`.
 *
 * @param rad O radicando (resultado da potência).
 * @param ind O índice (expoente da potência).
 * @return base A base (resultado do radical).
 */
double nrt(double rad, int ind);
//-----------------------------------------------------------------------------
/**
 * @brief Retorna o expoente (índice) de uma potência que tem base `rt` e resultado `rad`.
 *
 * @param rad O radicando (resultado da potência).
 * @param rt A raiz (base da potência).
 * @return ind O índice (expoente da potência).
 */
int indOfRoot(int rad, double rt);
//-----------------------------------------------------------------------------
/**
 * @brief Cria e retorna um polinômio com um número especificado de termos.
 *
 * Esta função aloca memória para uma estrutura de polinômio e inicializa
 * seu array interno de termos com base no número especificado de termos. Ela é
 * usada como ponto de partida para definir ou manipular polinômios.
 *
 * @param num_of_terms O número de termos que o polinômio deve conter.
 * @return Uma estrutura polinomial com espaço alocado para o número especificado de termos.
 */
polynomial pCreate(int num_of_terms);
//-----------------------------------------------------------------------------

/**
 * @brief Calcula o maior divisor comum (MDC) de dois inteiros usando o algoritmo de Euclides.
 *
 * @param a O primeiro inteiro.
 * @param b O segundo inteiro.
 * @return O maior divisor comum entre a e b.
 */
int gcd(int a, int b);
//-----------------------------------------------------------------------------
/**
 * @brief Retorna o n-ésimo número primo, onde n é menor ou igual a 41.
 *
 * Esta função só suporta valores de n até 41 porque o 41º número primo é 179,
 * que é o maior número primo cujo quadrado (179^2 = 32.041) ainda cabe em um inteiro com sinal de 32 bits.
 *
 * @param num A posição do número primo a ser retornado (deve ser ≤ 41).
 * @return O n-ésimo número primo, ou 0 se num > 41.
 */
int isPrime(int num);
//-----------------------------------------------------------------------------
/**
 * @brief Aplica o critério de Eisenstein para determinar se um polinômio é irredutível.
 *
 * Verifica se o polinômio satisfaz as condições do Critério de Eisenstein para irredutibilidade
 * usando o número primo p. De acordo com o critério, todos os coeficientes, exceto o líder,
 * devem ser divisíveis por p, e o termo constante não deve ser divisível por p².
 *
 * Embora o Critério de Eisenstein seja uma ferramenta poderosa, ele tem limitações — por exemplo,
 * o polinômio (x^4 + 1) é irredutível sobre os racionais, apesar de não satisfazer o critério.
 *
 * @param x O polinômio a ser testado.
 * @param p Um número primo usado para o teste de Eisenstein.
 * @return 1 se o polinômio for irredutível pelo critério de Eisenstein,
 *        -1 caso contrário.
 */
int eisenstein(polynomial x, double p);
//-----------------------------------------------------------------------------
/**
 * @brief Imprime todas as raízes complexas não reais de um polinômio, assumindo módulo unitário.
 *
 * Esta função assume que o polinômio tem apenas raízes complexas localizadas no círculo unitário
 * e calcula as raízes complexas na forma exponencial e^{iθ} e e^{-iθ}, onde θ é um ângulo
 * em radianos representado como múltiplo de π. Ela lida com os casos quando o polinômio é
 * monômio ou não e ajusta o número de raízes impressas conforme necessário.
 *
 * As raízes são exibidas usando a fórmula de Euler e mostradas como expressões fatoradas.
 *
 * @param p O polinômio do qual as raízes complexas serão extraídas e impressas.
 */
void complexp(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Gera uma representação simbólica fatorada das raízes de um polinômio do segundo grau.
 *
 * Esta função recebe os coeficientes e o discriminante de uma equação quadrática e retorna
 * uma string alocada dinamicamente contendo a fatoração simbólica das raízes. Ela trata diversos casos:
 *  - Raízes reais e racionais
 *  - Raízes irracionais (com radicais)
 *  - Raízes complexas conjugadas
 *  - Raízes puramente imaginárias
 *
 * Também é compatível com casos onde o polinômio representa um quadrado de binômio,
 * útil em fatorações de polinômios de grau 4.
 *
 * @warning A string retornada é alocada dinamicamente com `malloc()` e deve ser liberada
 * manualmente com `free()` para evitar vazamento de memória.
 *
 * @param aexp Expoente do termo principal do polinômio (normalmente 2 ou 4).
 * @param den Denominador utilizado na fórmula de Bhaskara (geralmente 2 * a).
 * @param b Coeficiente do termo linear do polinômio.
 * @param delta Discriminante da equação quadrática (b² - 4ac).
 * @param rootsPair Vetor de dois doubles contendo as raízes numéricas aproximadas (usado apenas para formatação).
 *
 * @return Ponteiro para uma string alocada dinamicamente contendo a fatoração simbólica.
 *         A responsabilidade de liberar a memória é do chamador da função.
 */
char* bhaskaraSimplify(int aexp, int den, int b, double delta, double* rootsPair);
//-----------------------------------------------------------------------------
/**
 * @brief Calcula e imprime a fatoração de um polinômio de segundo grau
 *        usando a fórmula de Bhaskara (quadrática).

 * Esta função simplifica os resultados dependendo do valor de `b` e de `delta`.
 *
 * @param p Um polinômio de grau 2 ou 4 com três termos na ordem decrescente
 *          de expoente.
 */
void bhaskara(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Aplica a identidade da soma ou diferença de potências (a^n ± b^n).
 *        Para n ímpar, tanto soma quanto diferença podem ser fatoradas.
 *        Para n par, apenas a diferença pode ser fatorada.
 *
 * Assume que o polinômio está na forma x^n ± a e imprime os passos de fatoração com base na fórmula.
 * Lida com expoentes pares e ímpares.
 *
 * @param p O polinômio a ser fatorado.
 */
void powerDiffFactoring(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Tenta fatorar um polinômio usando o Teorema das Raízes Racionais e divisão sintética.
 *
 * Aplica o método de Briot–Ruffini (divisão sintética) para encontrar raízes inteiras e racionais de um polinômio,
 * assim como raízes aproximadas de radicais quando possível. Exibe os termos fatorados lineares e radicais.
 *
 * @param p O polinômio a ser fatorado.
 */
void briot_ruffini(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Imprime a fatoração de um polinômio ou exibe a mensagem apropriada
 *        com base no tipo de fatoração que pode ser realizada nele.
 *
 * Esta função avalia um polinômio para determinar a técnica de fatoração apropriada:
 * - Se o polinômio é irredutível com base no critério de Eisenstein, ela exibe uma mensagem
 *   indicando que não há fatoração não-imaginária ou raízes representáveis por radicais.
 * - Se o polinômio é quadrático, aplica a fatoração por diferença de potências, se aplicável.
 * - Se o polinômio tem grau par e três ou menos termos, aplica o método de Bhaskara.
 * - Se o polinômio não atender a nenhum dos critérios acima, tenta o método Briot-Ruffini
 *   ou realiza fatoração complexa, se necessário.
 *
 * A função lida com polinômios com número variável de termos e coeficientes,
 * incluindo verificação de casos especiais, como coeficientes únicos e avaliação do critério de Eisenstein.
 *
 * @param p O polinômio a ser fatorado.
 */
void printFac(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Ao contrário de `printFac()`, esta função também desaloca a memória dos
 * termos de `p` usando `free(p.terms)`.
 *
 * @param p O polinômio a ser fatorado.
 */
 void pFree(polynomial p);

#endif // POLYNOMIALSZ_H_INCLUDED
