#include <iostream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <immintrin.h>  // cabecalho para instrucoes AVX2

// ESTRUTURA DA MATRIZ ESPARSA
// identica a versao nao otimizada

struct SparseMatrix {
    int n, num_diags;
    int *offsets;
    double **diags;
};

// TEMPO EM MILISSEGUNDOS

double getTime() {
    auto agora = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(agora.time_since_epoch()).count();
}


// GERACAO DA MATRIZ (igual a versao nao otimizada)
// a geracao nao e o gargalo, entao nao precisa ser otimizada

SparseMatrix* gerarMatriz(int n, int num_diags, unsigned int seed) {
    srand(seed);
    int half = (num_diags - 1) / 2;
    SparseMatrix *mat = new SparseMatrix{n, num_diags, new int[num_diags], new double*[num_diags]};

    for (int d = 0; d < num_diags; d++) {
        mat->offsets[d] = d - half;
        mat->diags[d] = new double[n]();
    }

    for (int d = half + 1; d < num_diags; d++) {
        int offset = mat->offsets[d], mirror = half - (d - half);
        for (int i = 0; i < n - offset; i++) {
            double valor = (double)rand() / RAND_MAX;
            mat->diags[d][i] = valor;
            mat->diags[mirror][i + offset] = valor;
        }
    }

    for (int i = 0; i < n; i++) {
        double soma = 0.0;
        for (int d = 0; d < num_diags; d++) {
            if (d != half && i + mat->offsets[d] >= 0 && i + mat->offsets[d] < n)
                soma += fabs(mat->diags[d][i]);
        }
        mat->diags[half][i] = soma + 1.0;
    }
    return mat;
}

double* gerarVetor(int n, unsigned int seed) {
    srand(seed + 1);
    double *v = new double[n];
    for (int i = 0; i < n; i++) v[i] = (double)rand() / RAND_MAX;
    return v;
}

void liberarMatriz(SparseMatrix *mat) {
    for (int d = 0; d < mat->num_diags; d++) delete[] mat->diags[d];
    delete[] mat->diags; delete[] mat->offsets; delete mat;
}


// OTIMIZACAO 1+2: MATVEC COM AVX2 + LOOP UNROLLING
// esta e a principal otimizacao

void matvec(SparseMatrix *mat, double *v, double *resultado) {
    int n = mat->n;

    // zera o vetor resultado
    for (int i = 0; i < n; i++) resultado[i] = 0.0;

    for (int d = 0; d < mat->num_diags; d++) {
        int offset = mat->offsets[d];
        double *diag = mat->diags[d];

        // calcula os limites seguros para AVX e unrolling
        // AVX processa 4 doubles por vez, entao precisamos de i+3 < n
        // e tambem j = i+offset deve estar entre 0 e n-1
        int i_start = (offset < 0) ? -offset : 0;
        int i_end   = (offset > 0) ? n - offset : n;

        // garante que temos multiplos de 4 para o bloco AVX
        int i_end_avx = i_start + ((i_end - i_start) / 4) * 4;

        // BLOCO AVX2: processa 4 doubles por iteracao
        for (int i = i_start; i < i_end_avx; i += 4) {
            int j = i + offset;

            // carrega 4 elementos de diag[i..i+3]
            // _mm256_loadu_pd: carrega 4 doubles sem exigir alinhamento
            __m256d a = _mm256_loadu_pd(&diag[i]);

            // carrega 4 elementos de v[j..j+3]
            __m256d b = _mm256_loadu_pd(&v[j]);

            // carrega 4 elementos atuais de resultado[i..i+3]
            __m256d r = _mm256_loadu_pd(&resultado[i]);

            // faz: r = r + a * b  (fused multiply-add)
            // _mm256_fmadd_pd: multiplica a*b e soma em r numa unica instrucao
            r = _mm256_fmadd_pd(a, b, r);

            // escreve os 4 resultados de volta em resultado[i..i+3]
            _mm256_storeu_pd(&resultado[i], r);
        }

        // LOOP UNROLLING: processa os elementos restantes (2 a 2)
        // os ultimos elementos que nao completam um grupo de 4
        int i_end_unroll = i_end_avx + ((i_end - i_end_avx) / 2) * 2;
        for (int i = i_end_avx; i < i_end_unroll; i += 2) {
            int j = i + offset;
            resultado[i]   += diag[i]   * v[j];
            resultado[i+1] += diag[i+1] * v[j+1];
        }

        // ESCALAR: processa o ultimo elemento restante (se houver)
        for (int i = i_end_unroll; i < i_end; i++) {
            resultado[i] += diag[i] * v[i + offset];
        }
    }
}

// OTIMIZACAO 3: DOT PRODUCT COM AVX2
// acumula 4 produtos por iteracao usando AVX


double dot(double *a, double *b, int n) {

    // registrador AVX que acumula 4 somas parciais simultaneamente
    __m256d soma_avx = _mm256_setzero_pd();  // inicializa com zeros

    // bloco AVX: processa 4 elementos por vez
    int i_end_avx = (n / 4) * 4;
    for (int i = 0; i < i_end_avx; i += 4) {
        __m256d va = _mm256_loadu_pd(&a[i]);  // carrega a[i..i+3]
        __m256d vb = _mm256_loadu_pd(&b[i]);  // carrega b[i..i+3]

        // soma_avx += va * vb  (4 multiplicacoes simultaneas)
        soma_avx = _mm256_fmadd_pd(va, vb, soma_avx);
    }

    // reduz os 4 valores acumulados em soma_avx para um unico double
    // _mm256_extractf128_pd: extrai metade superior do registrador AVX
    __m128d low  = _mm256_castpd256_pd128(soma_avx);    // metade inferior
    __m128d high = _mm256_extractf128_pd(soma_avx, 1);  // metade superior
    __m128d soma128 = _mm_add_pd(low, high);             // soma as duas metades

    // soma os dois doubles restantes no registrador de 128 bits
    __m128d soma128_shift = _mm_unpackhi_pd(soma128, soma128);
    double soma = _mm_cvtsd_f64(_mm_add_pd(soma128, soma128_shift));

    // processa os elementos restantes que nao completam grupo de 4
    for (int i = i_end_avx; i < n; i++)
        soma += a[i] * b[i];

    return soma;
}

double norma(double *v, int n) { return sqrt(dot(v, v, n)); }

// GRADIENTE CONJUGADO - logica identica, operacoes otimizadas

double* gradienteConjugado(SparseMatrix *mat, double *b, double tol, int max_iter,
                           double &t_matvec, double &t_dot, int &iter) {
    int n = mat->n;
    t_matvec = t_dot = 0.0;

    double *x  = new double[n]();
    double *r  = new double[n]();
    double *p  = new double[n]();
    double *Ap = new double[n]();

    for (int i = 0; i < n; i++) r[i] = p[i] = b[i];

    double t0 = getTime();
    double rtr = dot(r, r, n);
    t_dot += getTime() - t0;

    iter = 0;

    while (iter < max_iter) {
        t0 = getTime(); matvec(mat, p, Ap); t_matvec += getTime() - t0;
        t0 = getTime(); double ptAp = dot(p, Ap, n); t_dot += getTime() - t0;

        double alpha = rtr / ptAp;
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        t0 = getTime(); double res = norma(r, n); t_dot += getTime() - t0;
        iter++;

        if (res < tol) {
            printf("  Convergiu na iteracao %d — residuo: %.4e\n", iter, res);
            break;
        }

        t0 = getTime(); double rtr_novo = dot(r, r, n); t_dot += getTime() - t0;
        double beta = rtr_novo / rtr;
        for (int i = 0; i < n; i++) p[i] = r[i] + beta * p[i];
        rtr = rtr_novo;
    }

    if (iter >= max_iter) printf("  Atingiu maximo de iteracoes (%d).\n", max_iter);

    delete[] r; delete[] p; delete[] Ap;
    return x;
}

// FUNCAO DE TESTE - identica ao main.cpp

void rodarTeste(int n, int num_diags, unsigned int seed) {
    double tol = 1e-10, t_matvec, t_dot;
    int max_iter = 100000, iteracoes;

    double t_inicio = getTime();
    SparseMatrix *mat = gerarMatriz(n, num_diags, seed);
    double *b = gerarVetor(n, seed);
    double t_geracao = getTime() - t_inicio;

    double t_cg_inicio = getTime();
    double *x = gradienteConjugado(mat, b, tol, max_iter, t_matvec, t_dot, iteracoes);
    double t_cg = getTime() - t_cg_inicio;

    // verificacao: ||b - Ax||
    double *Ax = new double[n]();
    for (int d = 0; d < mat->num_diags; d++) {
        int offset = mat->offsets[d];
        for (int i = 0; i < n; i++) {
            int j = i + offset;
            if (j >= 0 && j < n) Ax[i] += mat->diags[d][i] * x[j];
        }
    }
    double residuo = 0.0;
    for (int i = 0; i < n; i++) { double diff = b[i] - Ax[i]; residuo += diff * diff; }
    residuo = sqrt(residuo);

    // calculo de memoria
    double bytes_double = sizeof(double);
    double mem_matriz   = (num_diags * sizeof(int) + num_diags * sizeof(double*) + num_diags * n * bytes_double);
    double mem_vetores  = (5 * n * bytes_double);
    double mem_total    = (mem_matriz + mem_vetores) / (1024.0 * 1024.0);

    double trafego_matvec = ((num_diags + 2) * n * bytes_double) / (1024.0 * 1024.0);
    double trafego_dot    = (2 * n * bytes_double) / (1024.0 * 1024.0);
    double trafego_iter   = trafego_matvec + (3 * trafego_dot) + ((9 * n * bytes_double) / (1024.0 * 1024.0));
    double trafego_total  = trafego_iter * iteracoes;

    printf("\n RESULTADOS:\n");
    printf("  Tamanho       : %d x %d\n", n, n);
    printf("  Bandas        : %d\n", num_diags);
    printf("  Geracao       : %.3f ms\n", t_geracao);
    printf("  Total CG      : %.3f ms\n", t_cg);
    printf("  - Matvec      : %.3f ms  (%.1f%%)\n", t_matvec, 100.0 * t_matvec / t_cg);
    printf("  - Dot         : %.3f ms  (%.1f%%)\n", t_dot,    100.0 * t_dot    / t_cg);
    printf("  Iteracoes     : %d\n", iteracoes);
    printf("  ||b - Ax||    : %.4e  %s\n", residuo, residuo < tol * 10 ? "OK" : "ATENCAO!");

    printf("\n ANALISE DE MEMORIA:\n");
    printf("  RAM Alocada Total : %.2f MB\n", mem_total);
    printf("  Trafego por Matvec: %.2f MB por chamada\n", trafego_matvec);
    printf("  Trafego por Dot   : %.2f MB por chamada\n", trafego_dot);
    printf("  Trafego Total CG  : %.2f MB movimentados\n", trafego_total);
    printf("  Largura de Banda  : %.2f GB/s\n", (trafego_total / 1024.0) / (t_cg / 1000.0));

    liberarMatriz(mat); delete[] b; delete[] x; delete[] Ax;
}

// MENU INTERATIVO

int main() {
    int tamanhos[] = {1024, 4096, 16384, 65536, 262144, 1048576}, opt_n, opt_b;

    printf("========================================\n  SOLVER OTIMIZADO - GRADIENTE CONJUGADO\n========================================\n"
           "\nEscolha o tamanho da matriz:\n  1 - 1024\n  2 - 4096\n  3 - 16384\n  4 - 65536\n  5 - 262144\n  6 - 1048576\n\nOpcao: ");
    if (!(std::cin >> opt_n) || opt_n < 1 || opt_n > 6) return printf("Opcao invalida!\n"), 1;

    printf("\nEscolha o numero de bandas:\n  1 - 7 bandas\n  2 - 27 bandas\n\nOpcao: ");
    if (!(std::cin >> opt_b) || opt_b < 1 || opt_b > 2) return printf("Opcao invalida!\n"), 1;

    int n = tamanhos[opt_n - 1], num_diags = (opt_b == 1) ? 7 : 27;
    printf("\n>>> Executando versao OTIMIZADA: n=%d | bandas=%d\n", n, num_diags);

    rodarTeste(n, num_diags, 42);
    return 0;
}