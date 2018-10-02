from collections import Counter
import  math
from sympy import Symbol, factor, expand, Rational
import numpy as np

# Generating All Partitions: A Comparison Of Two Encodings
# https://arxiv.org/abs/0909.2331
def accel_asc(n):
    a = [0 for unused_variable in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

# Symmetries and Groups Michaelmas, Term 2008, Hugh Osborn
# http://www.damtp.cam.ac.uk/user/ho/GNotes.pdf
def casimir(a):
    num = 0
    for i in range(0,len(a)):
        for j in range(0,a[i]):
            num = num + 2*(j - i)
    return num

# This function prints values of Casimir operator on each conjugacy class of symmetric group S_{n}
def printAllCasimirs(n):
    partitions = list(accel_asc(n))
    for i in range(0,len(partitions)):
        print('CasimirValue(',partitions[i],') = ',casimir(partitions[i]))
        
# This function returns total number of boxes in Young diagram
def NormOfYoungDiagram(a):
    sum = 0
    for i in range(0,len(a)):
        sum = sum + a[i]
    return sum

def printAllNormOfYoungDiagram(n):
    partitions = list(accel_asc(n))
    for i in range(0,len(partitions)):
        print('Norm(',partitions[i],') = ',NormOfYoungDiagram(partitions[i]))

# ON THE HECKE ALGEBRAS AND THE COLORED HOMFLY POLYNOMIAL
# https://arxiv.org/pdf/math/0601267.pdf
def firstTermOfHOMPLY(k,r,setOfColours):
    t_power = 0
    nu_power = k * (r - 1)/2
    for i in range(0, len(setOfColours)):
        t_power = t_power + casimir(setOfColours[i])            
    t_power = t_power * k * r/2
    print('t^{',t_power, '} *', '{/nu}^{',nu_power,'}')
    
def denominatorSizeOfConjClass(youngDiag):
    a = Counter(youngDiag)
    con = 1
    for j in range(0,NormOfYoungDiagram(youngDiag)+1):
        con = con * ((j**(a[j])) * math.factorial(a[j]))
    return con

def sizeOfConjClass(youngDiag):
    res = math.factorial(NormOfYoungDiagram(youngDiag))//denominatorSizeOfConjClass(youngDiag)
    return res

def printSizesOfConjClasses(number):
    qqq = list(accel_asc(number))
    summ = 0
    for n in range(0,len(qqq)):
        summ = summ + sizeOfConjClass(qqq[n])
        print("ConjClassSize(",qqq[n],") =",sizeOfConjClass(qqq[n]))  
    print(summ," -vs- ",math.factorial(number))
    
def productTermInSchurLatex(youngDiag):
    con='Z'
    for i in range(0,len(youngDiag)):
        if con == 'Z':
            con = '\\frac{q^{'+str(youngDiag[i])+'/2}-q^{-'+str(youngDiag[i])+'/2} }{t^{'+str(youngDiag[i])+'/2}-t^{-'+str(youngDiag[i])+'/2} }'
        else:    
            con = con+'\\cdot\\frac{q^{'+str(youngDiag[i])+'/2}-q^{-'+str(youngDiag[i])+'/2} }{t^{'+str(youngDiag[i])+'/2}-t^{-'+str(youngDiag[i])+'/2} }'
    return con #it returns answer in latex format

def getNumberBeforeSymbol(starting_index_in_file,symbol,charater_data):
    count=''
    i = starting_index_in_file
    while charater_data[i] != symbol:
        count = count + charater_data[i]
        i = i + 1
    return [int(count), i]

def loadOfCharacterTableToMemory(number):
    qqq = list(accel_asc(number))
    f = open('CharacterTable/characterTable'+str(number)+'.txt', 'r+')
    charater_data = f.read()
    f.close()
    #a = [[0 for unused_variable in range(0, len(qqq))] for unused_variable in range(0, len(qqq))]
    a = np.zeros(shape=(len(qqq), len(qqq)), dtype=np.int64)
    index_in_file = 1
    while charater_data[index_in_file] != '}':
        i_index_info = getNumberBeforeSymbol(index_in_file,',', charater_data)
        i_index = i_index_info[0]
        index_in_file = i_index_info[1] + 1
        j_index_info = getNumberBeforeSymbol(index_in_file,':', charater_data)
        j_index = j_index_info[0]
        index_in_file = j_index_info[1] + 1
        a_ij_info = getNumberBeforeSymbol(index_in_file, ',', charater_data)
        i_index = i_index - 1
        j_index = j_index - 1
        a[i_index][j_index] = a_ij_info[0]
        index_in_file = a_ij_info[1]
        index_in_file += 1
    print('finally')
    return a

def productTermInSchur(youngDiag):
    con = 1
    t = Symbol('t')
    q = Symbol('n')
    for i in range(0,len(youngDiag)):
        con = con * (q**(Rational(youngDiag[i], 2))-q**(-Rational(youngDiag[i],2)))/(t**(Rational(youngDiag[i], 2))-t**(- Rational(youngDiag[i], 2)))
    res = factor(con)
    return res #it returns answer

def quasiHadamarProduct(table, r):
    res = [0 for unused_variable in range(0, len(table))]
    for i in range(0,len(res)):
       res[i] = table[i] * r
    return res

def specialSchurPolynomial(number, index):
    count = 0
    qqq = list(accel_asc(number))
    charaters = loadOfCharacterTableToMemory(number)
    for i in range(0,len(qqq)):
        count = count + Rational(charaters[len(qqq)-index-1][len(qqq)-i-1], denominatorSizeOfConjClass(qqq[i]))*productTermInSchur(qqq[i])
    return count
