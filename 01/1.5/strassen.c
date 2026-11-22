#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXN 1024
#define RUNS 5
#define THRESHOLD 64

static int n;
static double a[MAXN][MAXN], b[MAXN][MAXN];
static double c[MAXN][MAXN], ref[MAXN][MAXN];

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void init(void) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            a[i][j] = (double)rand() / RAND_MAX;
            b[i][j] = (double)rand() / RAND_MAX;
        }
}

/* Эталон: обычное умножение, нужно для проверки и как базовый случай. */
static void naive(double out[MAXN][MAXN]) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0;
            for (int k = 0; k < n; k++)
                s += a[i][k] * b[k][j];
            out[i][j] = s;
        }
}

/* Умножение небольших блоков с произвольным шагом между строками. */
static void multiply_base(const double *A, int lda,
                          const double *B, int ldb,
                          double *C, int ldc, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++)
            C[i * ldc + j] = 0.0;
        for (int k = 0; k < size; k++) {
            const double aik = A[i * lda + k];
            for (int j = 0; j < size; j++)
                C[i * ldc + j] += aik * B[k * ldb + j];
        }
    }
}

/* Сложить или вычесть два блока: D = X + sign * Y. */
static void combine(const double *X, int ldx,
                    const double *Y, int ldy,
                    double *D, int size, int sign) {
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            D[i * size + j] = X[i * ldx + j] + sign * Y[i * ldy + j];
}

static double *alloc_block(int size) {
    double *p = malloc((size_t)size * size * sizeof *p);
    if (!p) {
        fprintf(stderr, "Не удалось выделить память для блока %dx%d\n", size, size);
        exit(EXIT_FAILURE);
    }
    return p;
}

/* Рекурсивное умножение Штрассена: C = A * B. */
static void strassen_rec(const double *A, int lda,
                         const double *B, int ldb,
                         double *C, int ldc, int size) {
    if (size <= THRESHOLD) {
        multiply_base(A, lda, B, ldb, C, ldc, size);
        return;
    }

    const int h = size / 2;
    const double *A11 = A,             *A12 = A + h;
    const double *A21 = A + h * lda,   *A22 = A + h * lda + h;
    const double *B11 = B,             *B12 = B + h;
    const double *B21 = B + h * ldb,   *B22 = B + h * ldb + h;

    double *M[7];
    for (int i = 0; i < 7; i++)
        M[i] = alloc_block(h);
    double *S = alloc_block(h);
    double *T = alloc_block(h);

    /* M1 = (A11 + A22)(B11 + B22) */
    combine(A11, lda, A22, lda, S, h, 1);
    combine(B11, ldb, B22, ldb, T, h, 1);
    strassen_rec(S, h, T, h, M[0], h, h);

    /* M2 = (A21 + A22) B11 */
    combine(A21, lda, A22, lda, S, h, 1);
    strassen_rec(S, h, B11, ldb, M[1], h, h);

    /* M3 = A11 (B12 - B22) */
    combine(B12, ldb, B22, ldb, T, h, -1);
    strassen_rec(A11, lda, T, h, M[2], h, h);

    /* M4 = A22 (B21 - B11) */
    combine(B21, ldb, B11, ldb, T, h, -1);
    strassen_rec(A22, lda, T, h, M[3], h, h);

    /* M5 = (A11 + A12) B22 */
    combine(A11, lda, A12, lda, S, h, 1);
    strassen_rec(S, h, B22, ldb, M[4], h, h);

    /* M6 = (A21 - A11)(B11 + B12) */
    combine(A21, lda, A11, lda, S, h, -1);
    combine(B11, ldb, B12, ldb, T, h, 1);
    strassen_rec(S, h, T, h, M[5], h, h);

    /* M7 = (A12 - A22)(B21 + B22) */
    combine(A12, lda, A22, lda, S, h, -1);
    combine(B21, ldb, B22, ldb, T, h, 1);
    strassen_rec(S, h, T, h, M[6], h, h);

    /* Собираем четыре блока результата. */
    for (int i = 0; i < h; i++)
        for (int j = 0; j < h; j++) {
            int p = i * h + j;
            C[i * ldc + j] = M[0][p] + M[3][p] - M[4][p] + M[6][p];
            C[i * ldc + j + h] = M[2][p] + M[4][p];
            C[(i + h) * ldc + j] = M[1][p] + M[3][p];
            C[(i + h) * ldc + j + h] = M[0][p] - M[1][p] + M[2][p] + M[5][p];
        }

    free(S);
    free(T);
    for (int i = 0; i < 7; i++)
        free(M[i]);
}

/* Результат алгоритма Штрассена записывается в c. */
static void strassen(void) {
    strassen_rec(&a[0][0], MAXN, &b[0][0], MAXN,
                 &c[0][0], MAXN, n);
}

static double maxdiff(void) {
    double d = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double x = fabs(c[i][j] - ref[i][j]);
            if (x > d) d = x;
        }
    return d;
}

int main(void) {
    const int sizes[] = {128, 256, 512, 1024};
    const int nsizes = (int)(sizeof sizes / sizeof sizes[0]);

    printf("n,method,run,seconds,gflops,maxdiff\n");
    for (int s = 0; s < nsizes; s++) {
        n = sizes[s];
        init();
        naive(ref);

        for (int m = 0; m < 2; m++) {
            const char *name = m ? "strassen" : "naive";
            for (int r = 1; r <= RUNS; r++) {
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        c[i][j] = 0.0;

                double t0 = now();
                if (m) strassen();
                else   naive(c);
                double t = now() - t0;

                printf("%d,%s,%d,%.6f,%.3f,%.3e\n",
                       n, name, r, t, 2.0 * n * n * n / t / 1e9, maxdiff());
            }
        }
    }
    return 0;
}
