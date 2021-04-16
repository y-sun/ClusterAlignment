#!/usr/bin/env python3

import os
import subprocess
import sys

ntotal = int(sys.argv[1])
cutoff= float(sys.argv[2])
change_type=0

def output_aligned(clique,ncluster):
    fp=open("position.xyz","r")
    natom=int(fp.readline().split()[0])
    fp.seek(0)

    jump=int((clique[0]-1)*(2*ncluster-clique[0])/2+0.5)*(natom+2)
    for i in range(jump):
        fp.readline()
    fout=open("position-maxclq.dat","w+")
    print((maxn-1)*natom,file=fout)
    print(maxn, clique[0],file=fout)
    j=1
    for i in range(ncluster-clique[0]):
        if(j >= maxn):break
        fp.readline()
        ll=fp.readline().split()
        if((int(ll[2])+1)==clique[0] and (int(ll[3])+1)==clique[j]):
            for k in range(natom):
                line=fp.readline()
                if(change_type==0):
                   print(line.strip("\n"),file=fout)
                else:
                   ll=line.split()
                   if(ll[0]=="1"):
                      print("Fe",ll[1],ll[2],ll[3],file=fout)
                   elif(ll[0]=="2"):
                    print("O",ll[1],ll[2],ll[3],file=fout)
                   else:
                      print("Wrong element!")
            j += 1
        else:
            for k in range(natom):
                fp.readline()
    fout.close()
    return



# get pair-wise score
matrix=[[0 for i in range(ntotal)] for j in range(ntotal)]
fpw=open("score.dat","r")
for line in fpw:
    ll=line.split()
    i=int(ll[0])
    j=int(ll[1])
    score=float(ll[2])
    matrix[i-1][j-1]=score
    matrix[j-1][i-1]=score
fpw.close()
#--

ite=0
pool=[i+1 for i in range(ntotal)]
pool_size=ntotal
pool_limit=5
while pool_size >pool_limit:
    # make graph
    fout=open("graph.txt","w+")
    print("p dege",pool_size,file=fout)
    space=0
    for i in range(pool_size):
        for j in range(pool_size):
            space+=1
            if(matrix[pool[i]-1][pool[j]-1] < cutoff):
                print("e",  i+1, j+1,file=fout)
    fout.close()
    proc=subprocess.Popen(["/home/yangsun/bin/mcqd.x graph.txt"],  stdout=subprocess.PIPE, shell=True)
    (output,err) = proc.communicate()
    #os.system("cp graph.txt graph-"+str(ite)+".txt")
    os.remove("graph.txt");

    fin=open("cliques.dat","r")
    maxn=0
    maxclq="0"
    for line in fin:
        ll=line.split()
        if(int(ll[0]) > maxn):
            maxn=int(ll[0])
            maxclq=line
    fin.close()
    fout=open("cliques-"+str(ite)+".dat","w+")
    print(maxn,end=" ",file=fout)
    ll=maxclq.split()
    for i in range(maxn):
        print(pool[int(ll[i+1])-1],end=" ",file=fout)
    fout.close()
    clique=[]
    for i in range(maxn):
        clique.append(pool[int(ll[i+1])-1])
    if(maxn>pool_limit):
        print(ite, maxn, pool[int(ll[1])-1])
        output_aligned(clique,ntotal)
        os.system("mv position-maxclq.dat position_clq_"+str(ite)+".xyz")
    else:
        os.system("rm cliques-"+str(ite)+".dat graph-"+str(ite)+".txt")
    for i in range(maxn):
        b=pool.index(clique[i])
        del pool[b]
    pool_size=len(pool)
    ite += 1
print(" ")
