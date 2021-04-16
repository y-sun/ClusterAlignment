#include "cluster.h"

void Atom::val(string ot, double ox, double oy, double oz) {
    x=ox;
    y=oy;
    z=oz;
    type=ot;
}


Cluster::Cluster() {
    m_size=1;
    m_atoms = new Atom [m_size];
    m_comment = "";
    m_bonds2 = 0;
}

Cluster::Cluster(int s, string words) {
    m_size = s;
    m_atoms = new Atom [m_size];
    m_comment = words;
}

Cluster::Cluster(const Cluster & st) { // deep copy
    m_size = st.m_size;
    m_atoms = new Atom [m_size];
    m_comment = st.m_comment;
    m_bonds2 = st.m_bonds2;
    for (int i=0; i<m_size; i++)
        read_atom(i, st.m_atoms[i] );
}

Cluster & Cluster::operator=(const Cluster & st) { // operator "=" for deep copy
    m_size = st.m_size;
    if (m_atoms) delete [] m_atoms;  // ADD THIS unless you want to explode your memory
    m_atoms = new Atom [m_size];
    m_comment = st.m_comment;
    m_bonds2 = st.m_bonds2;
    for (int i=0; i<m_size; i++)
        read_atom(i, st.m_atoms[i] );
    return *this;
}

Cluster::~Cluster() {
    if(m_atoms)
        delete [] m_atoms;
}

void Cluster::update(int s, string words) {
    m_size = s;
    m_atoms = new Atom [m_size];
    m_comment = words;
}


void Cluster::read_atom(int elem, Atom val) {
    m_atoms[elem].x = val.x;
    m_atoms[elem].y = val.y;
    m_atoms[elem].z = val.z;
    m_atoms[elem].bond2 = val.bond2;
    m_atoms[elem].type = val.type;
    m_atoms[elem].pair = val.pair;
    m_atoms[elem].iscore = val.iscore;
}

void Cluster::bond_sum() {
    m_bonds2=0;
    for (int i=0; i<m_size; i++) {
        m_atoms[i].bond2 = m_atoms[i].x*m_atoms[i].x +
                           m_atoms[i].y*m_atoms[i].y +
                           m_atoms[i].z*m_atoms[i].z;
        m_bonds2 +=  m_atoms[i].bond2;
    }
}

void Cluster::print() {
    std::cout << m_size << std::endl;
    std::cout << m_comment << std::endl;
    for (int i=0; i < m_size ; i++) {
        std::cout << m_atoms[i].type << " "
             << m_atoms[i].x << " "
             << m_atoms[i].y << " "
             << m_atoms[i].z << std::endl;
    }
}

void Cluster::writefile( std::ofstream &fout) {
    fout << m_size << std::endl;
    fout << m_comment << std::endl;
    for (int i=0; i < m_size ; i++) {
        fout << m_atoms[i].type << " "
             << m_atoms[i].x << " "
             << m_atoms[i].y << " "
             << m_atoms[i].z << " "
             << m_atoms[i].pair << " "
             << m_atoms[i].iscore << std::endl;
    }
}

void Cluster::comment( string words) {
    m_comment = m_comment + " " + words;
}


void Cluster::rotate(double ux, double uy, double uz, double theta) {
   // normalize 
   double u=sqrt(ux*ux+uy*uy+uz*uz);
   ux = ux/u; uy = uy/u; uz = uz/u;
   theta = theta / 180 * M_PI;

   // rotation matrix from axis and angle
   // https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
   double R[][3] = {
       { cos(theta)+ux*ux*(1-cos(theta)), ux*uy*(1-cos(theta))-uz*sin(theta), ux*uz*(1-cos(theta))+uy*sin(theta) },
       { uy*ux*(1-cos(theta))+uz*sin(theta), cos(theta)+uy*uy*(1-cos(theta)), uy*uz*(1-cos(theta))-ux*sin(theta) },
       { uz*ux*(1-cos(theta))-uy*sin(theta), uz*uy*(1-cos(theta))+ux*sin(theta), cos(theta)+uz*uz*(1-cos(theta)) } };

   for (int i=0; i<m_size; i++) {
       double x = m_atoms[i].x*R[0][0] + m_atoms[i].y*R[0][1] + m_atoms[i].z*R[0][2],
              y = m_atoms[i].x*R[1][0] + m_atoms[i].y*R[1][1] + m_atoms[i].z*R[1][2],
              z = m_atoms[i].x*R[2][0] + m_atoms[i].y*R[2][1] + m_atoms[i].z*R[2][2];
       m_atoms[i].x = x;
       m_atoms[i].y = y;
       m_atoms[i].z = z;
   }
}

void Cluster::translate(double ux, double uy, double uz, double dis) {
    // normalize
    double u=sqrt(ux*ux+uy*uy+uz*uz);
    ux = ux/u; uy = uy/u; uz = uz/u;
    for (int i=0; i<m_size; i++) {
        m_atoms[i].x += ux*dis;
        m_atoms[i].y += uy*dis;
        m_atoms[i].z += uz*dis;
    }
}

void Cluster::resize(double scale) {
    for (int i=0; i<m_size; i++) {
        m_atoms[i].x *= scale;
        m_atoms[i].y *= scale;
        m_atoms[i].z *= scale;
        m_atoms[i].bond2 = scale*scale*m_atoms[i].bond2;
    }
    m_bonds2 = scale*scale*m_bonds2;
}

void Cluster::recenter() {
    double x = m_atoms[0].x, y = m_atoms[0].y, z = m_atoms[0].z;
    for (int i=0; i<m_size; i++) {
        m_atoms[i].x -= x;
        m_atoms[i].y -= y;
        m_atoms[i].z -= z;
    }
}

void Cluster::cut(int natom) {
    m_size = natom;
    bond_sum();
}

void Cluster::mirror() {
    for (int i=0; i<m_size; i++) 
        m_atoms[i].x *= (-1.0);
}
