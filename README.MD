# cgSolver - Solver de Sistemas Lineares por Gradiente Conjugado 

**Integrante:** Helias Torres

## Descrição do Projeto 
Programa em C++ que resolve sistemas de equações lineares Ax=b usando o Método do Gradiente Conjugado, com suporte a matrizes esparsas representadas por diagonais.Inclui medição de desempenho e análise de memória.

---

## Estrutura de Arquivos 

| Arquivo | Descrição |
| :--- | :--- |
| **`main.cpp`** | Versão não otimizada do solver (compilada com -00). |
| **`mainoptimized.cpp`** | Versão otimizada com AVX2 e loop unrolling (compilada com -02 -mavx2 -mfma). |
| **`Makefile`** | Compila ambas as versões e gera os binários cgSolver e cgSolverOpt. |
| **`relatorio.pdf`** | Relatório completo com descrição, tabelas de desempenho e análise. |
| **`README.md`** | Este arquivo - instruções de uso e informações do projeto. |
| **`testes/run_tests.sh`** | Script para rodar todos os casos de teste automaticamente. |
| **`testes/resultados_nao_otimizado.txt`** | Saída completa dos testes da versão não otimizada. |
| **`testes/resultados_otimizado.txt`** | Saída completa dos testes da versão otimizada. |

---

## Requisitos

* **Sistema operacional:** Linux ou WSL (Windows Subsystem for Linux).
* **Compilador:** g++ com suporte a C++11 ou superior.
* **Processador:** Suporte a AVX2 (Intel Haswell 2013+ ou AMD Ryzen 1000+).

> **Nota:** Para verificar o suporte AVX2 no Linux, execute o comando:
> ```bash
> grep -o 'avx[^]*' /proc/cpuinfo | sort -u

---

## Como Compilar

Para compilar ambas as versões de uma vez:
> ```bash 
> make all
Para compilar individualmente.
> ```bash
> make cgSolver

>```bash
> make cgSolverOpt
Para remover os binários:
>```bash
>make clean

---

## Como executar
Ambos os programas possuem um menu interativo. Execute e siga as instruções na tela: 

Versão não otimizada:
>```bash
>./cgSolver  
  
Versão otimizada:
>```bash
> ./cgSolverOpt  

O programa solicitará:
1. O tamanho da matriz (1=1024, 2=4096, 3=16384, 4=65536, 5=262144, 6=1048576)
2. O número de bandas (1=7 bandas, 2=27 bandas)
Após a execução, o programa exibe os resultados e encerra.  
  
---

## Parâmetros Fixos

- Seed do gerador aleatório: 42 (garante reprodutibilidade)
- Tolerância de convergência: 1e-10
- Máximo de iterações: 100.000

---

## Progresso do Trabalho
  
- Estrutura da matriz esparsa por diagonais 
- Geração aleatória de matriz simétrica e positiva definida 
- Implementação do Método do Gradiente Conjugado 
- Instrumentalização com medição de tempo (chrono) 
- Análise analítica de memória (RAM alocada, tráfego, largura de banda) 
- Menu interativo para seleção de caso de teste 
- Versão otimizada com AVX2, loop unrolling e flag -O2 
- Testes para todos os 12 casos (6 tamanhos × 2 bandas) 
- Relatório com análise comparativa  

---

## Observações

- A versão otimizada (cgSolverOpt) requer processador com suporte a AVX2. Caso
contrário, o programa não compilará com as flags -mavx2 -mfma.
- Os testes foram realizados sempre na mesma máquina para garantir
comparabilidade dos resultados.
- A seed 42 foi usada em todos os testes. Alterar a seed gera uma matriz diferente
mas igualmente válida.