#!/usr/bin/env python3

import matplotlib.pyplot as plt

den_cut=16
delta_cut=2
draw_cut=1.0

fin=open("density-delta.dat","r")
fout=open("extracted.xyz","w+")
den=[]
delta=[]
p=[]
dis=[]
sid=[]
cn=0
for line in fin:
    ll=line.split()
    cn+=1
    if(float(ll[1]) > 0):
        den.append(float(ll[0]))
        delta.append(float(ll[1]))
        x=[]
        for k in range(3):
            x.append(float(ll[2+k]))
        p.append(x)
        dis.append(float(ll[5]))
        sid.append(cn)
fin.close()
print("Searching high density point(HDP) with Density >",den_cut, " Delta > ", delta_cut)
#fig, ax = plt.subplots()
#ax.scatter(den, delta,".")
#for i , txt in enumerate(sid):
#    ax.annotate( txt, (den[i],delta[i]))
#plt.savefig("hdp.png")

# output extracted.xyz
count=0
pool=[]
flag=[0 for i in range(len(den))]
for i in range(len(den)):
    if( den[i]>den_cut and delta[i]> delta_cut):
        add=1
        for pp in pool:
            dis=(pp[0]-p[i][0])**2+(pp[1]-p[i][1])**2+(pp[2]-p[i][2])**2
            if(dis < draw_cut**2):
                add=0
        if(add==1):
            count +=1
            flag[i]=1
            pool.append(p[i])

print(count,file=fout)
print("den >",den_cut, "delta >",delta_cut,file=fout)
for i in range(len(den)):
    if( flag[i]==1):
        print("T", p[i][0],p[i][1],p[i][2],file=fout)
fout.close()
print("Got", count, "HDP postions, written in extracted.xyz")

# merge with CHGCAR-test
fin=open("CHGCAR-test","r")
fout=open("CHG-extracted.vasp","w+")
for i in range(5):
    print(fin.readline().strip("\n"),file=fout)
print("Al",file=fout)
print(count,file=fout)
print("Cartesian",file=fout)
center=7.5
for i in range(len(den)):
    if( flag[i]==1):
        print(p[i][0]+center,p[i][1]+center,p[i][2]+center,file=fout)
for k in range(3):
    fin.readline()
for line in fin:
    print(line.strip("\n"),file=fout)
fout.close()
print("Merge HDP position with CHGCAR-test, written in CHG-extracted.vasp")
