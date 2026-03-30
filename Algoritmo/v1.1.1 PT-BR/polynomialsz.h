#ifndef POLYNOMIALSZ_H_INCLUDED
#define POLYNOMIALSZ_H_INCLUDED

/**
 * @file polynomialsz.h
 * @brief Arquivo de cabeçalho para as funcionalidades de fatoração de polinômios.
 *
 * Este arquivo fornece as estruturas de dados e declarações necessárias para a
 * fatoração de polinômios, incluindo:
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h> /**< Necessário para operações matemáticas. */

#define _USE_MATH_DEFINES /**< Habilita o uso de M_PI e outras constantes em alguns compiladores. */

// polynomialsz.h

/**
 * @brief Variável simbólica usada em expressões polinomiais (padrão: 'x')
 */
extern char var;

/**
 * @brief Expoente do fator x extraído por divideX()
 */
extern int degreeX;

/**
 * @brief Fator MDC extraído por divideGCD()
 */
extern int divider;

/**
 * @brief Flag que indica se a fatoração foi completa (1) ou não (0)
 */
extern int sol;

/**
 * @struct term
 * @brief Representa um único termo em um polinômio.
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
 * @brief Representa um polinômio como um vetor de termos.
 *
 * Os termos não precisam estar ordenados, mas geralmente estão em ordem decrescente de expoentes.
 */
typedef struct spolynomial {
    term *terms;        /**< Vetor dinâmico de termos do polinômio. */
    int numTerms;       /**< Número de termos no polinômio. */
} polynomial;

//-----------------------------------------------------------------------------
/**
 * @brief Cria e define o coeficiente e o expoente de um termo do polinômio.
 *
 * @param coef O coeficiente do termo.
 * @param exp O expoente do termo.
 * @return O termo completo.
 */
term setTerms(int coef, int exp);
//-----------------------------------------------------------------------------
/**
 * @brief Retorna a base de uma potência com expoente igual a ind, que resulta em rad.
 *
 * @param rad O radicando (resultado da potência).
 * @param ind O índice (expoente da potência).
 * @return base A base (resultado do radical).
 */
double nrt(double rad, int ind);
//-----------------------------------------------------------------------------
/**
 * @brief Retorna o expoente (índice) de uma potência que tem base rt e resultado rad.
 *
 * @param rad O radicando (resultado da potência).
 * @param rt A raiz (base da potência).
 * @return ind O índice (expoente da potência).
 */
int indOfRoot(int rad, double rt);
//-----------------------------------------------------------------------------
/**
 * @brief Cria e retorna um polinômio com um número específico de termos.
 *
 * Esta função aloca memória para uma estrutura de polinômio e inicializa
 * seu vetor interno de termos com base no número especificado de termos.
 * É usada como ponto de partida para definir ou manipular polinômios.
 *
 * @param numTerms O número de termos que o polinômio deve conter.
 * @return Uma estrutura de polinômio com espaço alocado para o número especificado de termos.
 */
polynomial pCreate(int numTerms);
//-----------------------------------------------------------------------------
/**
 * @brief Imprime o polinômio.
 *
 * @param p O polinômio.
 */
void pPrint(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Calcula o máximo divisor comum (MDC) de dois inteiros usando o algoritmo de Euclides.
 *
 * @param a O primeiro inteiro.
 * @param b O segundo inteiro.
 * @return O máximo divisor comum de a e b.
 */
int gcd(int a, int b);
//-----------------------------------------------------------------------------
/**
 * @brief Gera uma representação em string da forma fatorada para uma equação
 *        quadrática com raízes inteiras (racionais) e discriminante positivo.
 *
 * @param aexp O expoente do primeiro termo (2 ou 4, indicando fatores lineares ou quadráticos)
 * @param den O valor do denominador (2a) da fórmula de Bhaskara
 * @param b O coeficiente b da equação quadrática
 * @param delta O valor do discriminante (b² - 4ac)
 * @param absDelta O valor absoluto de delta
 * @param bSimplify Buffer de string pré-alocado para armazenar o resultado
 * @param powerRoot String indicando o formato do expoente da raiz ("^2" para quartico, vazio para quadrático)
 * @return Ponteiro para o buffer de string formatado
 */
char* intSrPositiveDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Gera uma representação em string da forma fatorada para uma equação
 *        quadrática com raízes inteiras (racionais) e discriminante negativo.
 *
 * @param aexp O expoente do primeiro termo (2 ou 4, indicando fatores lineares ou quadráticos)
 * @param den O valor do denominador (2a) da fórmula de Bhaskara
 * @param b O coeficiente b da equação quadrática
 * @param delta O valor do discriminante (b² - 4ac)
 * @param absDelta O valor absoluto de delta
 * @param bSimplify Buffer de string pré-alocado para armazenar o resultado
 * @param powerRoot String indicando o formato do expoente da raiz ("^2" para quartico, vazio para quadrático)
 * @return Ponteiro para o buffer de string formatado
 */
char* intSrNegativeDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Gera uma representação em string da forma fatorada para uma equação
 *        quadrática com raízes irracionais e discriminante positivo.
 *
 * @param aexp O expoente do primeiro termo (2 ou 4, indicando fatores lineares ou quadráticos)
 * @param den O valor do denominador (2a) da fórmula de Bhaskara
 * @param b O coeficiente b da equação quadrática
 * @param delta O valor do discriminante (b² - 4ac)
 * @param absDelta O valor absoluto de delta
 * @param bSimplify Buffer de string pré-alocado para armazenar o resultado
 * @param powerRoot String indicando o formato do expoente da raiz ("^2" para quartico, vazio para quadrático)
 * @return Ponteiro para o buffer de string formatado
 */
char* floatSrPositiveDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Gera uma representação em string da forma fatorada para uma equação
 *        quadrática com raízes complexas irracionais e discriminante negativo.
 *
 * @param aexp O expoente do primeiro termo (2 ou 4, indicando fatores lineares ou quadráticos)
 * @param den O valor do denominador (2a) da fórmula de Bhaskara
 * @param b O coeficiente b da equação quadrática
 * @param delta O valor do discriminante (b² - 4ac)
 * @param absDelta O valor absoluto de delta
 * @param bSimplify Buffer de string pré-alocado para armazenar o resultado
 * @param powerRoot String indicando o formato do expoente da raiz ("^2" para quartico, vazio para quadrático)
 * @return Ponteiro para o buffer de string formatado
 */
char* floatSrNegativeDelta(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Determina a função de formatação apropriada para raízes racionais
 *        com base no sinal do discriminante.
 *
 * @param aexp O expoente do primeiro termo (2 ou 4)
 * @param den O valor do denominador (2a) da fórmula de Bhaskara
 * @param b O coeficiente b da equação quadrática
 * @param delta O valor do discriminante (b² - 4ac)
 * @param absDelta O valor absoluto de delta
 * @param bSimplify Buffer de string pré-alocado para armazenar o resultado
 * @param powerRoot String indicando o formato do expoente da raiz
 * @return Ponteiro para o buffer de string formatado
 */
char* rationalSqRoots(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Determina a função de formatação apropriada para raízes irracionais
 *        com base no sinal do discriminante.
 *
 * @param aexp O expoente do primeiro termo (2 ou 4)
 * @param den O valor do denominador (2a) da fórmula de Bhaskara
 * @param b O coeficiente b da equação quadrática
 * @param delta O valor do discriminante (b² - 4ac)
 * @param absDelta O valor absoluto de delta
 * @param bSimplify Buffer de string pré-alocado para armazenar o resultado
 * @param powerRoot String indicando o formato do expoente da raiz
 * @return Ponteiro para o buffer de string formatado
 */
char* irrationalSqRoots(int aexp, int den, int b, double delta, int absDelta, char* bSimplify, char* powerRoot);
//-----------------------------------------------------------------------------
/**
 * @brief Função principal de simplificação que coordena a formatação de raízes
 *        racionais/irracionais com base nas propriedades do discriminante.
 *
 * Esta função verifica se o discriminante é um quadrado perfeito para determinar
 * se as raízes são racionais ou irracionais, então chama as funções de formatação
 * apropriadas.
 *
 * @param aexp O expoente do primeiro termo (2 ou 4)
 * @param den O valor do denominador (2a) da fórmula de Bhaskara
 * @param b O coeficiente b da equação quadrática
 * @param delta O valor do discriminante (b² - 4ac)
 * @param absDelta O valor absoluto de delta
 * @param rootsPair Vetor contendo as duas raízes calculadas
 * @return Ponteiro para o buffer de string formatado (deve ser liberado pelo chamador)
 */
char* bhaskaraSimplify(int aexp, int den, int b, double delta, int absDelta, double* rootsPair);
//-----------------------------------------------------------------------------
/**
 * @brief Calcula e imprime a fatoração de um polinômio de segundo grau
 *        usando a fórmula de Bhaskara.
 *
 * Esta função simplifica o resultado dependendo dos valores de "b" e delta.
 *
 * @param p Um polinômio de grau 2 ou 4 com três termos em ordem decrescente
 *          de expoente.
 */
void bhaskara(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Fatora um polinômio usando o Teorema das Raízes Racionais e divisão sintética.
 *
 * Aplica o método de Briot-Ruffini para encontrar raízes inteiras no intervalo
 * [-√maxNum, √maxNum]. Quando um fator quadrático é encontrado, chama bhaskara()
 * para completar a fatoração (incluindo raízes complexas).
 *
 * @param p Polinômio a ser fatorado (modificado durante o processo)
 * @param maxNum Maior valor absoluto dos coeficientes do polinômio original
 *
 * @note O polinômio de entrada é reduzido à medida que raízes são encontradas.
 * @note A variável global sol indica se a fatoração foi completa.
 */
void briotRuffini(polynomial p, int maxNum);
//-----------------------------------------------------------------------------
/**
 * @brief Função principal de fatoração que coordena todas as estratégias.
 *
 * Orquestra o processo completo de fatoração na seguinte ordem:
 * 1. Extrai fator x (se houver) via divideX() [recursivo]
 * 2. Extrai MDC dos coeficientes via divideGCD()
 * 3. Se for quadrático/biquadrático: usa bhaskara()
 * 4. Se os coeficientes forem sequenciais: usa briotRuffini()
 * 5. Se os coeficientes forem ±1: usa cyclotomicFac()
 * 6. Caso contrário: usa aberth() para aproximação numérica
 *
 * @param p Polinômio a ser fatorado
 *
 * @note A variável global sol indica se a fatoração foi completa
 * @note As variáveis globais divider e degreeX armazenam fatores extraídos
 */
void fac(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Imprime o polinômio seguido de sua fatoração.
 *
 * @param p Polinômio a ser impresso e fatorado
 */
void printFac(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Diferente de `printFac()`, esta alternativa também desaloca a memória dos termos
 * do polinômio usando `free(p.terms)`.
 *
 * @param p O polinômio a ser fatorado.
 */
 void pFree(polynomial p);
 //-----------------------------------------------------------------------------
 /**
 * @brief Extrai o fator x^(minExp) do polinômio e chama fac() recursivamente.
 *
 * Encontra o menor expoente entre todos os termos, divide todos os termos
 * por x^(minExp) e chama fac() no polinômio reduzido. O fator extraído é
 * armazenado na variável global degreeX para impressão posterior.
 *
 * @param p Polinômio original (não modificado diretamente)
 *
 * @note A variável global degreeX armazena o expoente do fator x extraído
 */
void divideX(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Divide todos os coeficientes do polinômio pelo MDC entre eles.
 *
 * Calcula o máximo divisor comum (MDC) de todos os coeficientes do polinômio
 * e divide cada coeficiente por esse valor. O fator extraído é armazenado
 * na variável global divider para impressão posterior.
 *
 * @param p Polinômio a ser simplificado (modificado in-place)
 *
 * @note A variável global divider armazena o MDC extraído
 */
void divideGCD(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Remove termos com coeficiente zero do polinômio e realoca memória.
 *
 * Percorre o vetor de termos, mantendo apenas aqueles com coeficiente != 0.
 * Se todos os termos forem zero, libera completamente o polinômio.
 * Realoca o vetor para o tamanho exato necessário.
 *
 * @param p Ponteiro para o polinômio a ser limpo (modificado in-place)
 */
void removeZeros(polynomial *p);
//-----------------------------------------------------------------------------
/**
 * @brief Fatora polinômios ciclotômicos (coeficientes ±1).
 *
 * Para polinômios da forma x^n + x^(n-1) + ... + x + 1 (ou variações com ±1),
 * encontra as raízes complexas no círculo unitário: exp(2πik/n).
 * Imprime os fatores lineares complexos usando notação exponencial.
 *
 * @param p Polinômio ciclotômico a ser fatorado
 *
 * @note Assume que o polinômio tem coeficientes apenas ±1
 * @note Imprime fatores no formato (x - Exp[i*π/m])
 */
void cyclotomicFac(polynomial p);
//-----------------------------------------------------------------------------
/**
 * @brief Aproxima todas as raízes de um polinômio usando o método de Aberth.
 *
 * Implementa o método de Aberth (ou Aberth–Ehrlich) para aproximação simultânea
 * de todas as raízes (reais e complexas) de um polinômio.
 *
 * Algoritmo:
 * 1. Converte coeficientes para double
 * 2. Estima raio máximo das raízes: R = 1 + max|coef[i]/coef[0]|
 * 3. Inicializa raízes igualmente espaçadas no círculo de raio R
 * 4. Itera a correção de Aberth até convergência:
 *    z_k^(novo) = z_k - [P(z_k)/P'(z_k)] / [1 - (P(z_k)/P'(z_k)) * Σ_{j≠k} 1/(z_k - z_j)]
 *
 * @param p Polinômio cujas raízes serão aproximadas
 *
 * @note Usa números complexos da biblioteca <complex.h>
 * @note Precisão de convergência: 1e-12
 * @note Máximo de 30 iterações
 */
void aberth(polynomial p);

#endif // POLYNOMIALSZ_H_INCLUDED
