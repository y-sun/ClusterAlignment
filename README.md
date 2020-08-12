# ClusterAlignment
Source code of Cluster Alignment method

## Compile
    make -f Makefile

## Run
    align.x para.in

## Parameters
Input parameters should be written into a parameter file, e.g. `para.in` here. It includes,

* TASK_TYPE {FCC, HCP, BCC, ICO(icosahedral), OCT(octahedral), TET(tetrahedral), 3661, 16661, 15551, other, PW}

If choosing "other", template information will be read from "template.dat". If choosing "PW", pair-wise alignment will be performed. With "PW", "EFFECT_NATOM" can be used to limit the atom number during alignment. 
     
* INPUT_FILE {string} : The name of file containing clusters to be aligned. 

* NCLUSTER {interge number} : Number of clusters to be aligned.
* TRANSLATE {0 or 1} : turn off/on the cluster translation
* RESIZE {0 or 1 or # > 1} : turn off/on the cluster scaler. If > 1, the template will be scaled multiple times and alignment will be redo with each templates.
* POSITION_WRITE {0 or 1} : turn off/on to write the aligned position
* MIRROR {0 or 1} : turn off/on to check the mirror image during the alignment 
* ELEMENT {0 or 1} : turn off/on to check the chemical element 
* RANDOM_TIME {interge number} : It is highly suggested to optimize the random rotation time before starting calculation. This factor directly scales the computational time linearly. Larger `RANDOM_TIME` leads to higher accuracy. The Lower symmetry template usually requires larger `RANDOM_TIME`.


## Reference
Please cite Sun *et al.*, Sci. Rep. 6,23734 (2016) if the resutls are used in your paper.

    @article{CA2016,
    title={‘Crystal genes’ in metallic liquids and glasses},
    author={Sun, Yang and Zhang, Feng and Ye, Zhuo and Zhang, Yue and Fang, Xiaowei and Ding, Zejun and Wang, Cai-Zhuang and Mendelev, Mikhail I and Ott, Ryan T and Kramer, Matthew J and others},
    journal={Scientific reports},
    volume={6},
    pages={23734},
    year={2016},
    publisher={Nature Publishing Group}

*Little hisgory* The idea of cluster alignment was initiated by [Kai-Ming Ho](https://scholar.google.com/citations?user=cGlRoOAAAAAJ&hl=en) and [Cai-Zhuang Wang](https://scholar.google.com/citations?user=9r-VpcgAAAAJ&hl=en) at [Ames laboraty](https://www.ameslab.gov/). The first code of cluster-template alignment was written with Fortran by [X.W. Fang](https://www.linkedin.com/in/%E5%B0%8F%E4%BC%9F-%E6%96%B9-0b9613b0/) *et al.* in 2009, published in [Phys. Rev. B 82,184204(2010)](https://doi.org/10.1103/PhysRevB.82.184204). [Y. Sun](https://scholar.google.com/citations?user=91yBLrMAAAAJ&hl=en) and [F. Zhang](https://scholar.google.com/citations?user=uL51e5oAAAAJ&hl=en) at Ames lab added pairwise alignment method as published in [Sun *et al.*, Sci. Rep. 6,23734 (2016)](https://doi.org/10.1038/srep23734). Later, the code was rewrote with C++ in 2016 to improve the optimization algorithm and efficiencies. Now the codes are mainly maintained by Y. Sun (yangsun017@gmail.com).
