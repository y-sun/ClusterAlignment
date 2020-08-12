/* ==========ALIGNMENT================ *
 * Written by Yang Sun at Ames Lab in  *
 * Oct. 2016.    yangsun017@gmail.com  *
 * Please cite Sci. Rep. 6,23734(2016) *
 * ----------------------------------- *
 * If anything is changed or needs to  *
 * be changed, plese let me know.      *
 * ----------------------------------- *
 * First ALIGNMENT code is written by  *
 * X.W. Fang in PRB 82,184204(2010)    *
 * =================================== *
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include "mpi.h"

#include "cluster.h"
#include "io.h"
#include "alignment.h"
#include "motif.h"
#include "version.h"

using namespace std;

int main(int argc, char* argv[]) {
    // mpi
    int numprocs, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    
    if(argc != 2 ) {
        if(myid==0)
            cout << "command: align.x para.in" << endl;
        return 0 ;
    }
    string filename=read_para(argv);
    extern int user_jobtype;
    extern int user_ncluster;
    extern int user_natom;
    extern string user_motif;
    extern bool user_position;
    extern int ran_ini;
    
    Cluster *clusters= new Cluster [user_ncluster];
    read_data(filename, user_ncluster, clusters);
    if(myid==0)
        mkdir("TMP", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(myid==0) {
        cout << "Alignment " << VERSION << " | Please cite [Sci. Rep. 6,23734] as reference" << endl;
        cout << "initial_random_number: " << ran_ini << endl;
    }
    
    switch (user_jobtype) {
        case 0:  // cluster-template alignment
        {
            // get template
            Cluster *temp=new Cluster [1];
            if(user_motif == "other" )
                read_data("template.dat", 1, temp);
            else
                motif(user_motif, temp[0]);

            // set two array for mpi
            double* score=new double [user_ncluster];
            double* score_out=new double [user_ncluster];
            for (int i=0; i<user_ncluster; i++) {score[i]=0; score_out[i]=0;}
            MPI_Barrier(MPI_COMM_WORLD);

            // alignment
            for (int i=0; i<user_ncluster ; i++) {
                if(myid == i%numprocs){
                    string comment=to_string((long long)i)+"  "+user_motif;
                    score[i] = align( clusters[i], temp[0], to_string((long long)i), comment );
                    cout << "cluster " << (i+1) << " to " << user_motif << " : " << score[i] <<endl;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // collect score from each cpu
            MPI_Reduce(score, score_out, user_ncluster, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            // output data 
            if(myid==0) {
                write_data(score_out);
                if(user_position) 
                    file_merge(user_ncluster);
            }

            // finish
            delete [] temp;
            break;
        }

        case 1:  // pair-wise alignment
        {
            int ntotal =  user_ncluster*(user_ncluster-1)/2;
            double* score=new double [ntotal];
            double* score_out=new double [ntotal];
            for (int i=0; i<ntotal; i++) {score[i]=0; score_out[i]=0;}
            MPI_Barrier(MPI_COMM_WORLD);
            
            int mpi_k=0;
            for (int i=0; i<user_ncluster-1; i++) 
                for (int j=i+1; j<user_ncluster; j++) {
                    if(myid == mpi_k%numprocs) {
                        int k = (2*user_ncluster - i-1)*i*0.5+(j-i) -1;
                        string comment=to_string((long long)i) + " " +to_string((long long)j)+" pw";
                        if( user_natom > 0) {  
                            // only user part of cluster to align
                            Cluster cluster_effect=clusters[i];
                            cluster_effect.cut(user_natom); 
                            comment += "_effect_atom_"+to_string((long long)user_natom);
                            score[k]= align( clusters[j] , cluster_effect, to_string((long long)k), comment ) ;
                        }
                        else
                            score[k]= align( clusters[j] , clusters[i], to_string((long long)k), comment ) ;
                        cout << "cluster " << (i+1) << " to cluster " << (j+1) <<  " : " << score[k] <<endl;
                    }
                    mpi_k += 1;
                }
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Reduce(score, score_out, ntotal, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myid==0) {
                write_data(score_out);
                if(user_position)
                    file_merge(ntotal);
            }
            break;
        }
    }
    delete [] clusters;
    MPI_Finalize();
}
