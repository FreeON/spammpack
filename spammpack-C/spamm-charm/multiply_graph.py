#!/usr/bin/env python

import networkx
import matplotlib.pyplot as plt

G = networkx.DiGraph()

N = 2

for i in range(N):
  for j in range(N):
    for k in range(N):
      if abs(i-k) <= 1 and abs(k-j) <= 1:
        A = "A{:d}{:d}".format(i, k)
        B = "B{:d}{:d}".format(k, j)
        C = "C{:d}{:d}".format(i, j)
        M = "M{:d}{:d}{:d}".format(i, j, k)
        G.add_node(A)
        G.add_node(B)
        G.add_node(C)
        G.add_node(M)
        G.add_edge(A, M)
        G.add_edge(B, M)
        G.add_edge(M, C)

networkx.draw_spring(G, dim = 8)
plt.show()
