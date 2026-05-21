#!/bin/bash
# ============================================================
# run_tests.sh — Executa todos os casos de teste automaticamente
# Uso: bash testes/run_tests.sh
# Certifique-se de ter compilado com: make all
# ============================================================

SOLVER="./cgSolver"
SOLVER_OPT="./cgSolverOpt"
OUT_NAO_OTIM="testes/resultados_nao_otimizado.txt"
OUT_OTIM="testes/resultados_otimizado.txt"

# Verifica se os binarios existem
if [ ! -f "$SOLVER" ]; then
  echo "Erro: $SOLVER nao encontrado. Execute 'make all' primeiro."
  exit 1
fi
if [ ! -f "$SOLVER_OPT" ]; then
  echo "Erro: $SOLVER_OPT nao encontrado. Execute 'make all' primeiro."
  exit 1
fi

echo "Iniciando testes..."

# Limpa arquivos de saida anteriores
> "$OUT_NAO_OTIM"
> "$OUT_OTIM"

# Tamanhos: 1=1024 2=4096 3=16384 4=65536 5=262144 6=1048576
# Bandas:   1=7    2=27

for tamanho in 1 2 3 4 5 6; do
  for bandas in 1 2; do

    # Nome amigavel para o log
    case $tamanho in
      1) N=1024    ;; 2) N=4096    ;; 3) N=16384   ;;
      4) N=65536   ;; 5) N=262144  ;; 6) N=1048576 ;;
    esac
    case $bandas in 1) B=7 ;; 2) B=27 ;; esac

    echo "Testando n=$N | bandas=$B ..."

    # Versao nao otimizada
    echo "\n=== n=$N | bandas=$B ===" >> "$OUT_NAO_OTIM"
    printf "%s\n%s\n" "$tamanho" "$bandas" | "$SOLVER" >> "$OUT_NAO_OTIM" 2>&1

    # Versao otimizada
    echo "\n=== n=$N | bandas=$B ===" >> "$OUT_OTIM"
    printf "%s\n%s\n" "$tamanho" "$bandas" | "$SOLVER_OPT" >> "$OUT_OTIM" 2>&1

  done
done

echo ""
echo "Testes concluidos!"
echo "  Resultados nao otimizados: $OUT_NAO_OTIM"
echo "  Resultados otimizados:     $OUT_OTIM"
