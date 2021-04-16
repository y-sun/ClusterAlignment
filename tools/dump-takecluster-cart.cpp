#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>

using namespace std;

typedef struct{
   int id;
   int type;
   float x;
   float y;
   float z;
}povec;


int main(int argc, char *argv[])
{
    // exe filename natom nstep cn
    if(argc != 7)
    {
       cout << "exe filename(with x y z) natom nstep cn njump nevery" << endl;
       return 0;
    }
    int nstep=atoi(argv[3]);
    int natom=atoi(argv[2]);
    int cn=atoi(argv[4]);
    int njump=atoi(argv[5]);
    int nevery=atoi(argv[6]);
    double cut=15;
   
    string fn1=string(argv[1]);
    ifstream fdump(fn1.c_str());
    string ss;
    double box[3];
    povec *pa=new povec[natom];
    double *diff=new double[natom];
    ofstream fout("cluster-todo.dat");

    for (int i=0; i<njump*(natom+9); i++)
       getline(fdump,ss);
    //nstep = int(nstep/nevery);

    for (int istep=0; istep<nstep; istep++)
    {
       if( istep%nevery != 0)
           for (int k=0; k<(natom+9); k++)
               getline(fdump,ss);
       else
       {
           getline(fdump,ss);
           int timestep;
           fdump >> timestep ; getline(fdump,ss);

           for ( int k=0; k<3; k++)
              getline(fdump, ss);
           double bl, bh;
           int aa, bb;
           for (int k=0; k<3; k++)
           {
              fdump >> bl >> bh ; getline(fdump,ss);
              box[k] = bh -bl;
           }
           getline(fdump,ss);
           for (int k=0; k<natom; k++)
           {
              fdump >> pa[k].id >> pa[k].type >> pa[k].x >> pa[k].y >> pa[k].z; 
              getline(fdump,ss);
           }

           for (int iatom=0; iatom<natom; iatom++)
           {
              int nbs[natom], count=0;
              double distance[natom];
              povec pdir[natom];
              for (int jatom=0; jatom<natom; jatom++)
              {
                 double dx=pa[jatom].x-pa[iatom].x,
                        dy=pa[jatom].y-pa[iatom].y,
                        dz=pa[jatom].z-pa[iatom].z;
                 if(dx/box[0]> 0.5) dx-=box[0];
                 if(dx/box[0]<-0.5) dx+=box[0];
                 if(dy/box[1]> 0.5) dy-=box[1];
                 if(dy/box[1]<-0.5) dy+=box[1];
                 if(dz/box[2]> 0.5) dz-=box[2];
                 if(dz/box[2]<-0.5) dz+=box[2];
                 //dx *= box[0]; dy *= box[1]; dz *= box[2];
                 double dis=dx*dx+dy*dy+dz*dz;
                 if( dis < cut*cut)
                 {
                    nbs[count]=jatom;
                    distance[count]=dis;
                    pdir[count].x=dx; pdir[count].y=dy; pdir[count].z=dz;
                    pdir[count].id=pa[jatom].id; pdir[count].type=pa[jatom].type;
                    count += 1;
                 }
              }
              if(count < cn) cout << "neighbor cutoff too small! " << count << endl;
              // sort
              for (int kk=0; kk<count-1; kk++)
                 for (int jj=kk; jj < count; jj++)
                 {
                    if(distance[kk] > distance[jj])
                    {
                       double dswap=distance[jj];
                       distance[jj]=distance[kk];
                       distance[kk]=dswap;
                       
                       int idx=nbs[jj];
                       nbs[jj]=nbs[kk];
                       nbs[kk]=idx;
                       
                       povec pp;
                       pp.x=pdir[jj].x; pp.y=pdir[jj].y; pp.z=pdir[jj].z; pp.id=pdir[jj].id; pp.type=pdir[jj].type;
                       pdir[jj].x=pdir[kk].x; pdir[jj].y=pdir[kk].y; pdir[jj].z=pdir[kk].z; pdir[jj].id=pdir[kk].id; pdir[jj].type=pdir[kk].type;
                       pdir[kk].x=pp.x; pdir[kk].y=pp.y; pdir[kk].z=pp.z;pdir[kk].id=pp.id; pdir[kk].type=pp.type;
                    }
                 }
              fout << cn << endl;
              fout << timestep << " " << iatom+1 << endl;
              cout << timestep << " " << iatom+1 << endl;
              for (int kk=0; kk<cn; kk++)
              {
                //fout << "AA " << pdir[kk].x << " " << pdir[kk].y <<  " " << pdir[kk].z << " " << pdir[kk].id << endl;
                 if(pdir[kk].type==1)
                    fout << "AA " << pdir[kk].x << " " << pdir[kk].y <<  " " << pdir[kk].z << " " << pdir[kk].id << endl;
                 else if(pdir[kk].type==2)
                    fout << "BB " << pdir[kk].x << " " << pdir[kk].y <<  " " << pdir[kk].z << " " << pdir[kk].id << endl;
              }
           }
       }
    }
    fout.close(); 
}
