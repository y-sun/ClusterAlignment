#include "motif.h"

void motif(string name, Cluster &A) {
    if( name == "FCC" ) {
        A.update(13,"FCC");
        Atom a;
        a.val("T", 0.0, 0.0, 0.0); A.read_atom(0,a);
        a.val("T",-2.0, 0.0,-2.0); A.read_atom(1,a);
        a.val("T",-2.0, 0.0, 2.0); A.read_atom(2,a);
        a.val("T",-2.0, 2.0, 0.0); A.read_atom(3,a);
        a.val("T", 0.0,-2.0,-2.0); A.read_atom(4,a);
        a.val("T", 0.0,-2.0, 2.0); A.read_atom(5,a);
        a.val("T", 2.0,-2.0, 0.0); A.read_atom(6,a);
        a.val("T", 0.0, 2.0,-2.0); A.read_atom(7,a);
        a.val("T", 2.0, 0.0,-2.0); A.read_atom(8,a);
        a.val("T",-2.0,-2.0, 0.0); A.read_atom(9,a);
        a.val("T", 0.0, 2.0, 2.0); A.read_atom(10,a);
        a.val("T", 2.0, 0.0, 2.0); A.read_atom(11,a);
        a.val("T", 2.0, 2.0, 0.0); A.read_atom(12,a);
        A.bond_sum();
    }
    else if( name == "HCP" ) {
        A.update(13,"HCP");
        Atom a;
        a.val("T", 0.000, 0.000000, 0.00000); A.read_atom(0,a);
        a.val("T", 0.000, 1.483790,-2.08528); A.read_atom(1,a);
        a.val("T", 0.000, 1.483790, 2.08528); A.read_atom(2,a);
        a.val("T", 1.285,-0.741895,-2.08528); A.read_atom(3,a);
        a.val("T", 1.285,-0.741895, 2.08528); A.read_atom(4,a);
        a.val("T",-1.285,-0.741895, 2.08528); A.read_atom(5,a);
        a.val("T",-1.285,-0.741895,-2.08528); A.read_atom(6,a);
        a.val("T", 1.285,-2.225685, 0.00000); A.read_atom(7,a);
        a.val("T",-1.285, 2.225685, 0.00000); A.read_atom(8,a);
        a.val("T",-1.285,-2.225685, 0.00000); A.read_atom(9,a);
        a.val("T", 1.285, 2.225685, 0.00000); A.read_atom(10,a);
        a.val("T", 2.570, 0.000000, 0.00000); A.read_atom(11,a);
        a.val("T",-2.570, 0.000000, 0.00000); A.read_atom(12,a);
        A.bond_sum();
    }
    else if( name == "BCC") {
        A.update(15,"BCC");
        Atom a;
        a.val("T", 0.00000, 0.00000, 0.00000); A.read_atom(0,a);
        a.val("T",-1.74826,-1.74826,-1.74826); A.read_atom(1,a);
        a.val("T",-1.74826,-1.74826, 1.74826); A.read_atom(2,a);
        a.val("T",-1.74826, 1.74826,-1.74826); A.read_atom(3,a);
        a.val("T",-1.74826, 1.74826, 1.74826); A.read_atom(4,a);
        a.val("T", 1.74826,-1.74826,-1.74826); A.read_atom(5,a);
        a.val("T", 1.74826,-1.74826, 1.74826); A.read_atom(6,a);
        a.val("T", 1.74826, 1.74826,-1.74826); A.read_atom(7,a);
        a.val("T", 1.74826, 1.74826, 1.74826); A.read_atom(8,a);
        a.val("T", 0.00000,-3.49652, 0.00000); A.read_atom(9,a);
        a.val("T", 0.00000, 0.00000,-3.49652); A.read_atom(10,a);
        a.val("T", 3.49652, 0.00000, 0.00000); A.read_atom(11,a);
        a.val("T", 0.00000, 0.00000, 3.49652); A.read_atom(12,a);
        a.val("T",-3.49652, 0.00000, 0.00000); A.read_atom(13,a);
        a.val("T", 0.00000, 3.49652, 0.00000); A.read_atom(14,a);
        A.bond_sum();
    }
    else if( name == "ICO") {
        A.update(13,"ICO");
        Atom a;
        double phi=(1+sqrt(5))*3/4;
        a.val("T", 0.00, 0.00, 0.00); A.read_atom(0,a);
        a.val("T", 0.00, 1.50,  phi); A.read_atom(1,a);
        a.val("T", 0.00, 1.50,0-phi); A.read_atom(2,a);
        a.val("T", 0.00,-1.50,  phi); A.read_atom(3,a);
        a.val("T", 0.00,-1.50,0-phi); A.read_atom(4,a);
        a.val("T", 1.50,  phi, 0.00); A.read_atom(5,a);
        a.val("T", 1.50,0-phi, 0.00); A.read_atom(6,a);
        a.val("T",-1.50,  phi, 0.00); A.read_atom(7,a);
        a.val("T",-1.50,0-phi, 0.00); A.read_atom(8,a);
        a.val("T",  phi, 0.00, 1.50); A.read_atom(9,a);
        a.val("T",0-phi, 0.00, 1.50); A.read_atom(10,a);
        a.val("T",  phi, 0.00,-1.50); A.read_atom(11,a);
        a.val("T",0-phi, 0.00,-1.50); A.read_atom(12,a);
        A.bond_sum();
    }
    else if( name == "TET" ){
        A.update(5,"TET");
        Atom a;
        a.val("T", 0.00, 0.00, 0.00); A.read_atom(0,a);
        a.val("T",-1.50,-1.50, 1.50); A.read_atom(1,a);
        a.val("T",-1.50, 1.50,-1.50); A.read_atom(2,a);
        a.val("T", 1.50,-1.50,-1.50); A.read_atom(3,a);
        a.val("T", 1.50, 1.50, 1.50); A.read_atom(4,a);
        A.bond_sum();
    }
    else if( name == "OCT") {
        A.update(7,"OCT");
        Atom a;
        a.val("T", 0.0, 0.0, 0.0); A.read_atom(0,a);
        a.val("T", 3.0, 0.0, 0.0); A.read_atom(1,a);
        a.val("T",-3.0, 0.0, 0.0); A.read_atom(2,a);
        a.val("T", 0.0, 3.0, 0.0); A.read_atom(3,a);
        a.val("T", 0.0,-3.0, 0.0); A.read_atom(4,a);
        a.val("T", 0.0, 0.0, 3.0); A.read_atom(5,a);
        a.val("T", 0.0, 0.0,-3.0); A.read_atom(6,a);
        A.bond_sum();
    }
    else if( name == "3661") {
        A.update(17,"3661");
        Atom a;
        a.val("T", 0.0, 0.0, 0.0); A.read_atom(0,a);
        a.val("T", 1.67682641, 0.09667704, 2.66250680); A.read_atom(1,a);
        a.val("T", 0.04534018,-3.05087948, 0.80650901); A.read_atom(2,a);
        a.val("T",-2.77714569,-1.32162396, 0.71984724); A.read_atom(3,a);
        a.val("T", 2.66127668,-1.59840304, 0.66326486); A.read_atom(4,a);
        a.val("T", 0.27451485, 3.14912045, 0.34652943); A.read_atom(5,a);
        a.val("T",-0.80199249, 1.63890407, 2.60779653); A.read_atom(6,a);
        a.val("T", 2.78875095, 1.47267791, 0.47222708); A.read_atom(7,a);
        a.val("T",-2.71214597,-0.00166187,-1.73022540); A.read_atom(8,a);
        a.val("T",-2.64967143, 1.74945698, 0.52880946); A.read_atom(9,a);
        a.val("T",-1.24153159, 2.30988989,-1.89241518); A.read_atom(10,a);
        a.val("T",-1.45476038,-2.44670257,-1.53643066); A.read_atom(11,a);
        a.val("T",-0.89485752,-1.34651662, 2.80987142); A.read_atom(12,a);
        a.val("T", 2.72627640,-0.27844095,-1.78680779); A.read_atom(13,a);
        a.val("T", 1.46918280, 2.16713407,-1.98136085); A.read_atom(14,a);
        a.val("T", 1.25627103,-2.58999153,-1.62493312); A.read_atom(15,a);
        a.val("T", 0.01440817,-0.27995239,-3.51732100); A.read_atom(16,a);
        A.bond_sum();
    }
    else if( name == "15551") {
        A.update(18,"15551");
        Atom a;
        a.val("T", 0.0, 0.0, 0.0); A.read_atom(0,a);
        a.val("T",  0.00000,  3.30000,  0.00000 ); A.read_atom(1,a);
        a.val("T",  3.13849,  1.01976,  0.00000); A.read_atom(2,a);
        a.val("T", -3.13849,  1.01976,  0.00000); A.read_atom(3,a);
        a.val("T",  1.93969, -2.66976,  0.00000); A.read_atom(4,a);
        a.val("T", -1.93969, -2.66976,  0.00000); A.read_atom(5,a);
        a.val("T",  0.00000, -2.52795,  2.12120); A.read_atom(6,a);
        a.val("T",  2.40422, -0.78118,  2.12120); A.read_atom(7,a);
        a.val("T", -2.40422, -0.78118,  2.12120); A.read_atom(8,a);
        a.val("T",  1.48589,  2.04515,  2.12120); A.read_atom(9,a);
        a.val("T", -1.48589,  2.04515,  2.12120); A.read_atom(10,a);
        a.val("T",  0.00000, -2.52795, -2.12120); A.read_atom(11,a);
        a.val("T",  2.40422, -0.78118, -2.12120); A.read_atom(12,a);
        a.val("T", -2.40422, -0.78118, -2.12120); A.read_atom(13,a);
        a.val("T",  1.48589,  2.04515, -2.12120); A.read_atom(14,a);
        a.val("T", -1.48589,  2.04515, -2.12120); A.read_atom(15,a);
        a.val("T",  0.00000,  0.00000, -3.50000); A.read_atom(16,a);
        a.val("T",  0.00000,  0.00000,  3.50000); A.read_atom(17,a);
        A.bond_sum();
    }
    else if( name == "16661") {
        A.update(21,"16661");
        Atom a;
        a.val("T", 0.0, 0.0, 0.0);           A.read_atom(0,a);
        a.val("T",  2.39,  -1.32,     2.39); A.read_atom(1,a);
        a.val("T", -2.39,  -1.32,    -2.39); A.read_atom(2,a);
        a.val("T",  0   ,  -2.58,    -2.00); A.read_atom(3,a);
        a.val("T",  0   ,  -2.58,     2.00); A.read_atom(4,a);
        a.val("T", -2.39,  -1.32,     2.39); A.read_atom(5,a);
        a.val("T",  2.39,  -1.32,    -2.39); A.read_atom(6,a);
        a.val("T",  0   ,   0   ,    -3.50); A.read_atom(7,a);
        a.val("T",  0   ,   0   ,     3.50); A.read_atom(8,a);
        a.val("T",  3.50,   0   ,     0   ); A.read_atom(9,a);
        a.val("T", -3.50,   0   ,     0   ); A.read_atom(10,a);
        a.val("T", -2.00,  -2.58,     0   ); A.read_atom(11,a);
        a.val("T",  2.00,  -2.58,     0   ); A.read_atom(12,a);
        a.val("T", -2.39,   1.32,    -2.39); A.read_atom(13,a);
        a.val("T", -2.00,   2.58,     0   ); A.read_atom(14,a);
        a.val("T", -2.39,   1.32,     2.39); A.read_atom(15,a);
        a.val("T",  0   ,   2.58,    -2.00); A.read_atom(16,a);
        a.val("T",  0   ,   2.58,     2.00); A.read_atom(17,a);
        a.val("T",  2.39,   1.32,    -2.39); A.read_atom(18,a);
        a.val("T",  2.00,   2.58,     0   ); A.read_atom(19,a);
        a.val("T",  2.39,   1.32,     2.39); A.read_atom(20,a);
        A.bond_sum();
    }
 



    else
        std::cout << "Cannot find information about Motif " 
                  << name << "." << std::endl;
    return;
}
