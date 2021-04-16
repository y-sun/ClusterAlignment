#!/usr/bin/env python3

import MD
import sys
import pylab as plt

if(len(sys.argv) != 3):
    print("exe file(only read second conlum) factor")
    exit()

factor=float(sys.argv[2])

fin=open(sys.argv[1],"r")
sc=[]
for line in fin:
    ll=line.split()
    sc.append(float(ll[-1])/factor)
fin.close()



q,f= MD.freq(sc, 0.005, 0, max(sc)+0.005)
fout=open("score-dis.dat","w+")
print("#score(rc="+sys.argv[2]+") freq",file=fout)


sss=0
for k in range(len(q)):
    sss += f[k]
    print(q[k],f[k]/sum(f), sss/sum(f),file=fout)
fout.close()

plt.plot(q,f)
plt.grid()
plt.savefig("his.png")
