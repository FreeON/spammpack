#N=8192
#B=128
N=4096
B=4
REPEAT=5

OPTIONS=(
    "-N ${N}"
    "-b ${B}"
    "-T exp_decay"
    "-l 0.994"
    "-R ${REPEAT}"
    "-c"
    "-v"
    )

THREADS=( 48 1 24 2 4 8 12 16 20 28 32 36 40 )
