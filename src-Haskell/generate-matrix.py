#!/usr/bin/env python

def main ():

    import argparse
    import random

    parser = argparse.ArgumentParser()

    parser.add_argument("-N",
            help = "Matrix size.",
            required = True,
            type = int
            )

    options = parser.parse_args()

    print("%%MatrixMarket matrix coordinate real general")
    print("%d %d %d" % (options.N, options.N, options.N**2))
    for i in range(options.N):
        for j in range(options.N):
            print("%d %d %e" % (i, j, random.random()))

if __name__ == "__main__":
    main()
