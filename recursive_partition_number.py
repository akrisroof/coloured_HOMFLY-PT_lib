# http://www.cs.utsa.edu/~wagner/python/fp/part.html
import sys

LIM = 1000000
pmem = [0.0]*(LIM+1)

def insertp(t, loc):
   pmem[loc] = t
   
def p(n):
    if n < 0:
        return 0
    if n == 0:
        return 1
    if pmem[n] != 0:
        return pmem[n]
    sign = 1
    res = 0
    for k in range(1, n+1):
        p1 = n - k*(3*k - 1)/2
        if p1 >= 0:
            term1 = p(p1)
            insertp(term1, p1)
        else:
            term1 = 0
        p2 = n - k*(3*k + 1)/2
        if p2 >= 0:
            term2 = p(p2)
            insertp(term2, p2)
        else:
            term2 = 0
        term = term1 + term2
        res = res + sign*term
        sign = -sign
    return res

text_file = open("Output.txt", "w")

print "LIM:", LIM
n = 0
while n <= LIM:
    res = p(n)
    if n%1 == 0:
        print n, res
        text_file.write("\n %s" % n)
        text_file.write(" %s" % res)
    n = n+1


text_file.close()

print "\n"
