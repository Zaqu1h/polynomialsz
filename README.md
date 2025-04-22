# 🧮 Polynomialsz: Um algoritmo em C de fatoração polinomial ![Linguagem](https://img.shields.io/badge/C-99-blue) ![Versão](https://img.shields.io/badge/versão-1.0.0-blue) 
---

## ✍🏽 Descrição
>  Este projeto conta com funções para encontrar raízes e fatorar polinômios. O código é escrito em linguagem C e possui documentações e descrições tanto em inglês quanto em português.
>  O projeto começou, originalmente, apenas como uma atividade complementar para a disciplina de Fundamentos Matemáticos para Ciência da Computação - FMCC -, em 2024. Todavia, após uma série
>  de melhorias, este pode servir como apoio para estudantes de disciplinas relacionáveis, como as de programação, cálculo, álgebra e, é claro, a própria disciplina de FMCC. 

---

## ✔️ Recursos
> O algoritmo tem suporte para a análise e fatoração de polinômios através de diversas técnicas, são elas 💁🏼:
>- Critério de Eisenstein
>- Briot-Ruffini (divisão sintética)
>- Método de Bhaskara (fórmula quadrática)
>- Fatoração por diferença de potências
>- Análise de raízes complexas

#### ⚠️ Observação: *Este algoritmo ainda apresenta incapacidade em fatorar alguns polinômios. Isso ocorre, na maioria das vezes, com polinômios com graus muito elevados ou com raízes complexas... Complexas.*

---

## 📝 Instruções de uso

Este repositório está organizado em versões, com suporte para português e inglês.  
Siga os passos abaixo para clonar o projeto e acessar os códigos-fonte
#### Observação: *Apesar do guia abaixo, você pode simplesmente executar os códigos na IDE de sua preferência, obviamente.*

### 🧬 Passo 1: Clone o repositório

Se você tem o Git instalado, use o comando:

```bash
git clone https://github.com/Zaqu1h/polynomialsz.git
```

### 📁 Passo 2: Acesse o arquivo desejado

> Algoritmo/
>
>├── v1.0.0 PTBR/   → Código e arquivo de cabeçalho com comentários em português
> 
>└── v1.0.0/        → Code and header file with english comments
>

Selecione a pasta da versão e o idioma que quiser.
```bash
cd polynomialsz/Algoritmo/v1.0.0 PTBR
```
### ⚙️ Passo 3: Compile o programa

Dentro da pasta, compile o programa com gcc...
```bash
gcc -o polynomialsz polynomialsz.c -lm
```
...E execute com:
```bash
./polynomialsz
```
### 🧪 Passo 4: Execute os testes
Se você deseja rodar os testes, eles estão disponíveis dentro da pasta Algoritmo/. Para compilar e rodar os testes, use o comando:
```bash
gcc -o teste1 polynomialsz.c teste1.c -lm
./teste1
```
#### Observação: *Há um arquivo teste2.c também!... Embora só tenha dois polinômios sendo testados nele...*

### ⌨️ Passo 5: Comece a digitar!
A utilização do código é mais fácil do que parece:
```bash

polynomial U = pCreate(3);  < Cria um polinômio com 3 termos

U.terms[0] = setTerms(1, 2);
U.terms[1] = setTerms(0, 1);
U.terms[2] = setTerms(4, 0);

printFac(U);  < Mostra o polinômio fatorado (ou pelo menos eu espero que sim)
               e NÃO desaloca a memória dos termos, diferente da função pFree().

  return 0;
```
O polinômio criado em questão foi: x² + 4 (ou x² + 0x + 4). O resultado apresentado, neste caso, será:
```bash
(x - i2)(x + i2)
```

---

## 📰 Documentações

> docs/
>
> ├── v1.0.0 PTBR/   → Documentação em Português
> 
> └── v1.0.0/        → English Documentation

As documentações foram feitas usando Doxygen. Você pode procurar por arquivos chamados "refman.rtf" ou "index.html".

---

## 📋 Planos | To-do List

> - Desenvolver e aplicar materiais e metodologia para que o projeto possa se tornar um bom apoio para aprendizado;
> - Melhorar o algoritmo para suporte a outros tipos de raízes complexas, além de alterar suas representações, se necessário;
> - Melhorar as condicionais do algoritmo para um desempenho mais satisfatório e uma leitura menos cansativa;
> - Adicionar mais métodos e critérios, se necessário, como Cardano ou Ferrari.

---

## ⁉️ Dúvidas, dicas sugestões

Sinta-se livre para me contactar pelo instagram ou pelo meu email, disponíveis a um clique de distância no fim deste README!

---

## ✨ Agradecimentos

A maioria não contribuiu diretamente, mas foram responsáveis pela minha persistência no projeto.

- Alba Silva Passos
- César Alberto Bravo Pariente - (Orientador);
- Esbel Tomás Valero Orellana - [@etvorellana](https://github.com/etvorellana) (Professor de Linguagem de Programação);
- Alex Amaral dos Santos - [@AlexAmaral1](https://github.com/AlexAmaral1)
- Larissa Brito - [@laribrito](https://github.com/laribrito?tab=overview&from=2024-12-01&to=2024-12-31)
- José Fernando Santos Rodrigues Júnior

---

## 🦱 Autor

Desenvolvido por Isaque Silva Passos Ribeiro - [@Zaqu1h](https://github.com/Zaqu1h) -, um estudante de Ciência da Computação na Universidade Estadual de Santa Cruz (UESC).

E-mail - isaquebinha03@hotmail.com

[![Instagram](https://img.shields.io/badge/Instagram-E4405F?style=for-the-badge&logo=instagram&logoColor=white)](https://www.instagram.com/zaquih_/)

