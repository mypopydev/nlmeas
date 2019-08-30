#!/usr/bin/python
# Coefficients of DCT-II as a matrix.
N = 8
denom = 2 * N
print "c_X = cos(pi * X / %d)" % denom
print "%"
print "% coefficient matrix:"
print "\\begin{bmatrix}"
for n in range(N):
  for k in range(N):
    num = k * (2*n + 1)
    print "%-7s"%("c_{%d}" % num),
    if k == N-1:
      print "\\\\"
    else:
      print "&",
print "\\end{bmatrix}"

print "%"
print "% within 90 degrees:"
print "%"
print "\\begin{bmatrix}"
for n in range(N):
  for k in range(N):
    num = k * (2*n + 1)
    while num > denom * 2:
      num -= denom * 2
    if num > denom:
      num = 2 * denom - num
    pre = " "
    if num > denom / 2:
      pre = "-"
      num = denom - num
    print "%-4s"%("%sc_%d" % (pre, num)),
    if k == N-1:
      print "\\\\"
    else:
      print "&",
print "\\end{bmatrix}"
