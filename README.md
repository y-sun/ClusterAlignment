# ClusterAlignment
Source code of Cluster Alignment method

## Compile
    make -f Makefile

## Run
    align.x para.in

## Parameters
Input parameters includes

```TASK_TYPE
     FCC, HCP, BCC, ICO(icosahedral), OCT(octahedral), TET(tetrahedral), 3661, 16661, 15551, other, PW
     If choosing "other", template information will be read from "template.dat".
     If choosing "PW", pair-wise alignment will be performed.  
     With "PW", "EFFECT_NATOM" can be used to limit the atom number during alignment. 
INPUT_FILE      
     The name of file containing clusters to be aligned. 
NCLUSTER        
     Number of clusters to be aligned.
TRANSLATE       
     0/1 : turn off/on the cluster translation
RESIZE          
     0/1 : turn off/on the cluster scaler 
POSITION_WRITE  
     0/1 : turn off/on to write the aligned position
MIRROR  
     0/1 : turn off/on to check the mirror image during the alignment 
ELEMENT  
     0/1 : turn off/on to check the chemical element 
RANDOM_TIME    
     It is highly suggested to optimize the random rotation time before starting calculation. 
     This factor will affect the computational time significantly. It mainly depends on 
     the symmetry of the motif(template) and the accuracy you want to get. 




## Reference
Please cite Sun et al, Sci. Rep. 6,23734 (2016) if the resutls are used in your paper.
The first Cluster Alignment code was written with Fortran by X.W. Fang in Phys. Rev. B 82,184204(2010)
