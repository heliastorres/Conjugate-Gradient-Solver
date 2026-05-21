#include <iostream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cstdio>

// Estrutura compacta da Matriz Esparsa (Ponteiro de Ponteiro)
struct SparseMatrix {
    int n, num_diags;
    int *offsets;
    double **diags;
};

// Retorna o tempo atual em segundos
double getTime() {
    auto agora = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(agora.time_since_epoch()).count();
}

// Aloca e gera a matriz esparsa por dominância diagonal (Simétrica Positiva Definida)
SparseMatrix* gerarMatriz(int n, int num_diags, unsigned int seed) {
    srand(seed);
    int half = (num_diags - 1) / 2;
    SparseMatrix *mat = new SparseMatrix{n, num_diags, new int[num_diags], new double*[num_diags]};
    
    for (int d = 0; d < num_diags; d++) {
        mat->offsets[d] = d - half;
        mat->diags[d] = new double[n]();
    }

    // Preenche a metade superior e espelha na inferior
    for (int d = half + 1; d < num_diags; d++) {
        int offset = mat->offsets[d], mirror = half - (d - half);
        for (int i = 0; i < n - offset; i++) {
            double valor = (double)rand() / RAND_MAX;
            mat->diags[d][i] = valor;
            mat->diags[mirror][i + offset] = valor;
        }
    }

    // Garante a dominância diagonal na diagonal principal
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

// Gera o vetor b aleatório
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

// Multiplicação Matriz-Vetor: resultado = A * v
void matvec(SparseMatrix *mat, double *v, double *resultado) {
    for (int i = 0; i < mat->n; i++) resultado[i] = 0.0;
    for (int d = 0; d < mat->num_diags; d++) {
        int offset = mat->offsets[d];
        for (int i = 0; i < mat->n; i++) {
            int j = i + offset;
            if (j >= 0 && j < mat->n) resultado[i] += mat->diags[d][i] * v[j];
        }
    }
}

// Produto escalar entre dois vetores
double dot(double *a, double *b, int n) {
    double soma = 0.0;
    for (int i = 0; i < n; i++) soma += a[i] * b[i];
    return soma;
}

double norma(double *v, int n) { return sqrt(dot(v, v, n)); }

// Algoritmo do Gradiente Conjugado 
double* gradienteConjugado(SparseMatrix *mat, double *b, double tol, int max_iter,
                                   double &t_matvec, double &t_dot, int &iter) {
    int n = mat->n;
    t_matvec = t_dot = 0.0;
    double *x = new double[n](), *r = new double[n](), *p = new double[n](), *Ap = new double[n]();

    for (int i = 0; i < n; i++) r[i] = p[i] = b[i];

    double t0 = getTime(); double rtr = dot(r, r, n); t_dot += getTime() - t0;
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

// Gerencia o fluxo de teste e validação das respostas
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

    // Verificação do resíduo local: ||b - Ax||
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

    // CÁLCULO ANALÍTICO DE MEMÓRIA 
    double bytes_double = sizeof(double);
    
    // Pegada de memória total alocada (RAM ocupada)
    double mem_matriz = (num_diags * sizeof(int) + num_diags * sizeof(double*) + num_diags * n * bytes_double);
    double mem_vetores = (5 * n * bytes_double); // b, x, r, p, Ap
    double mem_total_alocada = (mem_matriz + mem_vetores) / (1024.0 * 1024.0); // em MB

    // Tráfego de memória por chamada/iteração (Volume de dados)
    double trafego_matvec = ((num_diags + 2) * n * bytes_double) / (1024.0 * 1024.0); // em MB
    double trafego_dot = (2 * n * bytes_double) / (1024.0 * 1024.0); // em MB
    
    // Tráfego total estimado no laço do CG
    // Cada iteração faz: 1 matvec, 3 dots/normas, e 3 atualizações de vetores (3*3 = 9 acessos por n)
    double trafego_por_iteracao = trafego_matvec + (3 * trafego_dot) + ((9 * n * bytes_double) / (1024.0 * 1024.0));
    double trafego_cg_total = trafego_por_iteracao * iteracoes;

    // Print unificado dos resultados finais
    printf("\n RESULTADOS:\n");
    printf("  Tamanho       : %d x %d\n", n, n);
    printf("  Bandas        : %d\n", num_diags);
    printf("  Geracao       : %.3f ms\n", t_geracao);
    printf("  Total CG      : %.3f ms\n", t_cg);
    printf("  - Matvec      : %.3f ms  (%.1f%%)\n", t_matvec, 100.0 * t_matvec / t_cg);
    printf("  - Dot         : %.3f ms  (%.1f%%)\n", t_dot, 100.0 * t_dot / t_cg);
    printf("  Iteracoes     : %d\n", iteracoes);
    printf("  ||b - Ax||    : %.4e  %s\n", residuo, residuo < tol * 10 ? "OK" : "ATENCAO!");
    
    printf("\n ANALISE DE MEMORIA:\n");
    printf("  RAM Alocada Total : %.2f MB\n", mem_total_alocada);
    printf("  Trafego por Matvec: %.2f MB lidos/escritos por chamada\n", trafego_matvec);
    printf("  Trafego por Dot   : %.2f MB lidos por chamada\n", trafego_dot);
    printf("  Trafego Total CG  : %.2f MB movimentados na RAM\n", trafego_cg_total);
    printf("  Largura de Banda  : %.2f GB/s\n", (trafego_cg_total / 1024.0) / (t_cg / 1000.0));

    liberarMatriz(mat); delete[] b; delete[] x; delete[] Ax;
}

int main() {
    int tamanhos[] = {1024, 4096, 16384, 65536, 262144, 1048576}, opt_n, opt_b;

    printf("========================================\n  SOLVER - GRADIENTE CONJUGADO\n========================================\n"
           "\nEscolha o tamanho da matriz:\n  1 - 1024\n  2 - 4096\n  3 - 16384\n  4 - 65536\n  5 - 262144\n  6 - 1048576\n\nOpcao: ");
    if (!(std::cin >> opt_n) || opt_n < 1 || opt_n > 6) return printf("Opcao invalida!\n"), 1;

    printf("\nEscolha o numero de bandas:\n  1 - 7 bandas\n  2 - 27 bandas\n\nOpcao: ");
    if (!(std::cin >> opt_b) || opt_b < 1 || opt_b > 2) return printf("Opcao invalida!\n"), 1;

    int n = tamanhos[opt_n - 1], num_diags = (opt_b == 1) ? 7 : 27;
    printf("\n>>> Executando: n=%d | bandas=%d\n", n, num_diags);
    
    rodarTeste(n, num_diags, 42);
    return 0;
}