#include "io.h"

int user_jobtype;
int user_ncluster;
int user_natom=-999;
int user_rant=200;
string user_motif;
bool user_position=1;
bool user_translate=1;
int user_resize=0;
bool user_mirror=0;
bool user_element=0;
int ran_ini;

string read_para(char* argv[]) {
    std::ifstream fp(argv[1]);
    string ss,filename;
    string job;
    while ( fp >> ss ) {
        if(ss=="TASK_TYPE")
            fp >> user_motif;
        if(ss=="EFFECT_NATOM")
            fp >> user_natom; 
        if(ss=="INPUT_FILE")
            fp >> filename;
        if(ss=="NCLUSTER")
            fp >> user_ncluster;
        if(ss=="TRANSLATE")
            fp >> user_translate;
        if(ss=="RESIZE")
            fp >> user_resize;
        if(ss=="RANDOM_TIME")
            fp >> user_rant;
        if(ss=="POSITION_WRITE")
            fp >> user_position;
        if(ss=="MIRROR")
            fp >> user_mirror;
        if(ss=="ELEMENT")
            fp >> user_element;
    }
    fp.close();
    if( user_motif == "FCC"  or
        user_motif == "HCP"  or
        user_motif == "BCC"  or
        user_motif == "ICO"  or
        user_motif == "OCT"  or
        user_motif == "TET"  or
        user_motif == "3661" or
        user_motif == "other")
        user_jobtype = 0;
    else if( user_motif == "PW")
        user_jobtype = 1;
    //int rand_init=time(NULL);
    //ran_ini = rand_init  ;
    //ran_ini =  (rand_init%1000)*10+ (rand_init%2)*10000;
    //std::cout << "initial_random_number: " << rand_init << std::endl;
    //srand(ran_ini);
    
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, 99999); // define the range

    ran_ini = distr(gen) ;

    return filename;
}

void read_data(std::string filename, int nc, Cluster *clusters) {
    std::ifstream fin;
    std::string ss;
    Atom a;
    int nt;
    fin.open(filename.c_str());
    for (int i=0; i<nc; i++) {
        fin >> nt ; getline(fin,ss);
        getline(fin,ss);
        clusters[i].update(nt,ss);
        for (int k=0; k<nt; k++){
            fin >> a.type >> a.x >> a.y >> a.z; getline(fin,ss);
            clusters[i].read_atom(k,a);
        }
        clusters[i].recenter();
        clusters[i].bond_sum();
    }
    fin.close();
    return;
}

void write_data(double* score) {
    std::ofstream fout("score.dat");
    switch (user_jobtype) {
        case 0:
        {
            for (int i=0; i<user_ncluster; i++) 
                fout << i+1 << "  " 
                     << std::fixed << std::setprecision(6) 
                     << score[i] << std::endl;
            break;
        }
        case 1:
        {
            for (int i=0; i<user_ncluster-1; i++)
                for(int j=i+1; j<user_ncluster; j++) {
                    int k = (2*user_ncluster - i-1)*i*0.5+(j-i) -1;
                    fout << i+1 << "  " << j+1 << "  " 
                        << std::fixed << std::setprecision(6) 
                        << score[k] << std::endl;
                }
            break;
        }
    }
    fout.close();
    return;
}

void file_merge(int nfiles) {
    std::ofstream fout("position.xyz");
    for (int i=0; i<nfiles; i++) {
        std::string filename="TMP/"+std::to_string((long long)i) + ".xyz";
        std::ifstream fin(filename.c_str());
        std::string ss;
        while ( std::getline(fin, ss) )
            fout << ss << std::endl;
        fin.close();
        std::remove(filename.c_str());
    }
    fout.close();
}
