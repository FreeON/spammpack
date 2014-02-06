#N=8192
#B=128
#LAMBDA=0.994

N=4096
B=4
LAMBDA=0.994

#N=2048
#B=4
#LAMBDA=0.99

REPEAT=5

OPTIONS=(
    "-N ${N}"
    "-b ${B}"
    "-T exp_decay"
    "-l ${LAMBDA}"
    "-R ${REPEAT}"
    "-c"
    "-v"
    )

THREADS=( 48 1 24 2 4 8 12 16 20 28 32 36 40 )
