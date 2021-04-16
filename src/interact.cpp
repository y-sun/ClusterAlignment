#include "cluster.h"

namespace lj {
    const double rcut=3.0;
    const double sigma=1.0;
    const double epsfactor=1.0;
    const double shift2=0.5;

    const double rcut2=rcut*rcut;
    const double sigma2=sigma*sigma;
    const double encut=-4.0*epsfactor*pow(sigma2/rcut2+shift2,3);
}

double Cluster::score (const Cluster & temp)  {
//   score between this Cluster and Cluster temp                                      
    for (int i=0; i<m_size; i++)
        m_atoms[i].iscore=999;
    if( m_size < temp.m_size ) {                                                      
        std::cout << "Error: too few atoms in the cluster" << std::endl;              
        return 9999;                                                                  
    }

    m_atoms[0].pair = 0; // center atom only choose template center
    for (int i=1; i<m_size; i++) {
        m_atoms[i].pair = -1;
    }
    

    int cal_size;
    if( m_size > (int)(1.5*temp.m_size) )
        // reduce calculating time for big clusters. Bond-sorting is required!
        cal_size = (int)(1.5*temp.m_size);
    else
        cal_size = m_size;
    
    extern bool user_element;

    int const narray=cal_size*temp.m_size;                                              
    double sc[narray];                                                                
    for (int i=1; i<cal_size; i++) // ignore center atoms                               
        for (int j=1; j<temp.m_size; j++) {                                           
            int p_sc=i*temp.m_size+j;                                                      
            double dx=m_atoms[i].x-temp.m_atoms[j].x,                                 
                   dy=m_atoms[i].y-temp.m_atoms[j].y,                                 
                   dz=m_atoms[i].z-temp.m_atoms[j].z;                                 
            sc[p_sc] = (dx*dx + dy*dy+ dz*dz)/temp.m_atoms[j].bond2;                                          
            if( user_element and (m_atoms[i].type != temp.m_atoms[j].type) ) 
                sc[p_sc] += 0.10;
        }                                                                             
    std::vector<bool> select_C(cal_size), select_T(temp.m_size);                        
    std::fill(select_C.begin(), select_C.end(), 0);                                   
    std::fill(select_T.begin(), select_T.end(), 0);                                   
    
    int Tstatus=0, flag=1;                                                            
    double sct=0;                                                                     
    do {
        double scmin=9999;                                                            
        int minT=0, minC=0;
        for (int j=1; j<temp.m_size; j++) {                                           
            if(select_T[j] == 1)                                                      
                continue;                                                             
            for (int i=1; i<cal_size; i++) {                                            
                if(select_C[i] == 1)
                    continue;
                int p_sc=i*temp.m_size+j;
                if( sc[p_sc] < scmin ) {
                    scmin = sc[p_sc];                                                 
                    minC = i;                                                         
                    minT = j;                                                         
                }
            }
        }

        select_C[minC]=1;
        select_T[minT]=1;
        sct += scmin;
        Tstatus += 1;
        m_atoms[minC].pair = minT;
        m_atoms[minC].iscore = scmin;
    } while (Tstatus < (temp.m_size-1) );
    double alignment_score = sqrt(sct/(temp.m_size-1));

    return sqrt(sct/(temp.m_size-1));
}


double Cluster::paired_score (const Cluster & temp, double* force, double* torque)  {
    int cal_size;
    for (int i=0; i<m_size; i++)
        m_atoms[i].iscore=999;
    
    if( m_size > (int)(1.5*temp.m_size) )
        // reduce calculating time for big clusters. Bond-sorting is required!
        cal_size = (int)(1.5*temp.m_size);
    else
        cal_size = m_size;

    extern bool user_element;
    double sc=0;
    for(int i=0; i<3; i++) { force[i]=0; torque[i]=0;}
    int count = 0;
    for (int i=0; i<cal_size; i++) {
        if(m_atoms[i].pair <= 0) 
            continue;
        count += 1;
        int j=m_atoms[i].pair;
        double dx=m_atoms[i].x-temp.m_atoms[j].x,                                 
               dy=m_atoms[i].y-temp.m_atoms[j].y,                                 
               dz=m_atoms[i].z-temp.m_atoms[j].z;                                
        double isc=(dx*dx + dy*dy+ dz*dz)/temp.m_atoms[j].bond2;
        if( user_element and (m_atoms[i].type != temp.m_atoms[j].type) ) 
            isc += 0.10;
        sc += isc;
        m_atoms[i].iscore=isc;
        double fa[3]={-1*dx, -1*dy, -1*dz};
        for (int k=0; k<3; k++)
            force[k] += fa[k];
        torque[0] += (m_atoms[i].y-temp.m_atoms[0].y)*fa[2]-(m_atoms[i].z-temp.m_atoms[0].z)*fa[1];
        torque[1] += (m_atoms[i].z-temp.m_atoms[0].z)*fa[0]-(m_atoms[i].x-temp.m_atoms[0].x)*fa[2];
        torque[2] += (m_atoms[i].x-temp.m_atoms[0].x)*fa[1]-(m_atoms[i].y-temp.m_atoms[0].y)*fa[0];
    }
    if( count != (temp.m_size -1 ))
        std::cout << "Error: number of aligned atoms!";
    return sqrt(sc/(temp.m_size-1));
}

double Cluster::bond2_compare (const Cluster & temp) {
    // using same number of bonds to compare
    double bonds2=0;
    for (int i=0; i< temp.m_size; i++)
        bonds2 += m_atoms[i].bond2;
    return temp.m_bonds2/bonds2;
}


