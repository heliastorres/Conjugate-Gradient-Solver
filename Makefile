CXX      = g++

# -O0: sem otimizacoes (versao nao otimizada)
# -O2: otimizacoes do compilador (versao otimizada)
# -mavx2: habilita instrucoes AVX2
# -mfma: habilita fused multiply-add (necessario para _mm256_fmadd_pd)

# compila versao nao otimizada
cgSolver: main.cpp
	$(CXX) -O0 -Wall -o cgSolver main.cpp

# compila versao otimizada
cgSolverOpt: mainoptimized.cpp
	$(CXX) -O2 -Wall -mavx2 -mfma -o cgSolverOpt mainoptimized.cpp

# compila as duas de uma vez
all: cgSolver cgSolverOpt

clean:
	rm -f cgSolver cgSolverOpt