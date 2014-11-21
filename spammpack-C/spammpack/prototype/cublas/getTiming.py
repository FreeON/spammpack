#!/usr/bin/python

import re, sys

# Open stdin for reading.
lines = sys.stdin.readlines()

DMA1 = []
DMA2 = []
blas = []
cublas = []

for line in lines:
  print line.rstrip()
  result = re.compile("([0123456789]+)x([0123456789]+) matrix").search(line)
  if result:
    N = int(result.group(1))

  result = re.compile("DMA1 time:.+= ([0123456789.e+-]+) s per").search(line)
  if result:
    DMA1.append(float(result.group(1)))

  result = re.compile("DMA2 time:.+= ([0123456789.e+-]+) s per").search(line)
  if result:
    DMA2.append(float(result.group(1)))

  result = re.compile("cublas time:.+= ([0123456789.e+-]+) s per").search(line)
  if result:
    cublas.append(float(result.group(1)))

  result = re.compile("blas time:.+= ([0123456789.e+-]+) s per").search(line)
  if result:
    blas.append(float(result.group(1)))

avg_DMA1 = 0.0
avg_DMA2 = 0.0
avg_blas = 0.0
avg_cublas = 0.0

for i in DMA1:
  avg_DMA1 += i

for i in DMA2:
  avg_DMA2 += i

for i in blas:
  avg_blas += i

for i in cublas:
  avg_cublas += i

avg_DMA1 /= float(len(DMA1))
avg_DMA2 /= float(len(DMA2))
avg_blas /= float(len(blas))
avg_cublas /= float(len(cublas))

print "average over %i trials, %ix%i matrix" % (len(DMA1), N, N)
print "avg(DMA1) = %e" % (avg_DMA1)
print "avg(DMA2) = %e" % (avg_DMA2)
print "avg(blas) = %e" % (avg_blas)
print "avg(cublas) = %e" % (avg_cublas)

print "%i %e %e %e %e" % (N, avg_DMA1, avg_DMA2, avg_blas, avg_cublas)
