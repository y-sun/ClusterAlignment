#!/usr/bin/env python3

import sys

fp=open("../pairwise-"+sys.argv[1]+"/para.in","r")

for line in fp:
    if("INPUT_FILE" in line):
        path=line.split()[1]

fxyz=open("../pairwise-"+sys.argv[1]+"/"+path,"r")
lines=fxyz.readlines()
natom=int(lines[0].split()[0])
ncluster=int(len(lines)/(natom+2))
clusters=[]
for i in range(ncluster):
    atoms=[]
    for j in range(natom+2):
        atoms.append(lines[ (natom+2)*i+j])
    clusters.append(atoms)

fin=open("cliques-out.dat","r")
for line in fin:
    if("Clique#" in line):
        break
ft=open("results.dat","w+")
print("Clique#  clique_size  first_cluster#  crystal# ", file=ft)
for line in fin:
    ll=line.split()
    print(line.strip("\n"), end=" ", file=ft)
    if(len(ll)==0):
        break
    cid=int(ll[2])
    fout=open(ll[0]+".xyz","w+")
    for k in range(natom+2):
        print(clusters[cid-1][k].strip("\n"),file=fout)
    print(clusters[cid-1][1].strip("\n"),file=ft)
    fout.close()

