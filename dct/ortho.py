#!/usr/bin/python
# Build matrix of DCT-II coefficients and check that rows and columns are
# orthonormal.
#
#       k=0  k=1  k=2  k=3  k=4  k=5  k=6  k=7
# n=0 [  c0   c1   c2   c3   c4   c5   c6   c7  ]
# n=1 [  c0   c3   c6  -c7  -c4  -c1  -c2  -c5  ]
# n=2 [  c0   c5  -c6  -c1  -c4   c7   c2   c3  ]
# n=3 [  c0   c7  -c2  -c5   c4   c3  -c6  -c1  ]
# n=4 [  c0  -c7  -c2   c5   c4  -c3  -c6   c1  ]
# n=5 [  c0  -c5  -c6   c1  -c4  -c7   c2  -c3  ]
# n=6 [  c0  -c3   c6   c7  -c4   c1  -c2   c5  ]
# n=7 [  c0  -c1   c2  -c3   c4  -c5   c6  -c7  ]
# 
import numpy as np

N = 8

def show(m):
  for j in range(N):
    print "[",
    for i in range(N):
      print "%7.4f" % (m[j,i]),
    print "]"

print "coefficients:"
coeff = np.zeros((N, N))
for n in range(N):
  for k in range(N):
    coeff[n,k] = np.cos(k * np.pi / N * (n + .5))

if True:
  # Scaling
  coeff[:,0] *= 1. / np.sqrt(2.)
  coeff[:,:] *= np.sqrt(2. / N)

show(coeff)

print "dot product of rows:"
rows = np.zeros((N, N))
for j in range(N):
  for i in range(N):
    rows[j,i] = np.dot( coeff[j,:], coeff[i,:] )
show(rows)

print "dot product of cols:"
cols = np.zeros((N, N))
for j in range(N):
  for i in range(N):
    cols[j,i] = np.dot( coeff[:,j], coeff[:,i] )
show(cols)
