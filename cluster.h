#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using std::string;

class Atom {
public:
    double x;
    double y;
    double z;
    double bond2; // square of bond length to (0,0,0)
    string type;
    void val(string ot, double ox, double oy, double oz);
    int pair;
    double iscore;
};


class Cluster {
private:
    int m_size;
    Atom *m_atoms;
    string m_comment;
    double m_bonds2;  // square summation of all the bond length
public:
    Cluster();
    Cluster (int s, string words);
    Cluster (const Cluster & st);  // deep copy
    Cluster & operator=(const Cluster & st); //deep copy
    ~Cluster ();
    void update(int s, string words);
    void read_atom ( int elem, Atom val );
    void bond_sum() ;
    void print();
    void writefile( std::ofstream &fout);
    void comment(string words);
    
    void rotate(double ux, double uy, double uz, double theta);
    void translate(double ux, double uy, double uz, double dis);
    void resize(double scale);
    void recenter();
    void cut(int natom);
    void mirror();
    
    //interact
    double score (const Cluster & temp) ;
    double paired_score (const Cluster & temp, double *force, double *torque);
    double bond2_compare (const Cluster & temp);
};

#endif
