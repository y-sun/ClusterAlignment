#include "alignment.h"

namespace para {
    const int nran=20000; // steps of random rotate
    const int nran_conv=500; // if score not change in this steps, break the random rotation
    const int nlo=1000;   // steps of local rotate
    const int nmn=100;   // steps of minimization
    const double rdeg=5.0; // degree
    
    const double ratio=0.7; // tuning factor for converge
    const double LSTEPSTART=0.1; // angstrom
    const int LMAX=200;     
    const double converge=1e-4;
    const double LEPSO=1e-8;

    const int opt=10;
    const double opt_start=0.9;
    const double opt_end=1.1;
    const double opt_delta=0.02;
}

double align ( Cluster & A_init, const Cluster & B_init, string tag, string words) {
    // initialization
    extern bool user_translate;
    extern int user_resize;
    extern bool user_mirror;
    extern bool user_position;
    extern int  user_rant;
    string filename="TMP/"+tag+".xyz";
    int nran_conv;
    if( user_rant > para::nran_conv) 
        nran_conv =  user_rant;
    else
        nran_conv = para::nran_conv;
    
    clock_t t0, td, t_rot=0, t_tran=0;

    //resize, if no resize, use one-element array.
    int nopt;
    //double delta=(para::opt_end - para::opt_start)/para::opt;
    double delta=para::opt_delta;
    double start=1-delta*int(user_resize/2);
    Cluster* B;
    
    if( user_mirror) {
        if( user_resize ==0 ){
            B = new Cluster [2];
            nopt = 2;
            B[0] = B_init;
            B[1] = B_init;
            B[1].mirror();
        }
        else if( user_resize > 0) {
            nopt = user_resize * 2;
            B = new Cluster [nopt];
            double fac = 1.0/sqrt( A_init.bond2_compare(B_init) );
            for (int i=0; i< user_resize; i++) {
                B[i*2] = B_init;
                B[i*2].resize( fac * (start+delta*i) );
                B[i*2+1] = B[i*2];
                B[i*2+1].mirror();
            }
        }
    }
    else {
        if( user_resize ==0 ){
            B =  new Cluster [1];
            nopt = 1;
            B[0] = B_init;
        }
        else if (user_resize > 0) {
            nopt = user_resize;
            B = new Cluster [nopt];
            double fac = 1.0/sqrt( A_init.bond2_compare(B_init) );
            for (int i=0; i< nopt; i++) {
                B[i] = B_init;
                B[i].resize( fac * (start+delta*i) );
            }
        }
    }
    
    // ALIGNMENT: 
    //      Random_rotate -> Translate (optional) -> Local_rotate
    //      Directions of translation & local rotation depende on force & torque.
    int iconv;
    double sc_f=9999;
    int best_opt;
    Cluster A_f;
    for (int iopt=0; iopt<nopt; iopt++) {
        double sc_ran=9999;
        Cluster A_ran_min;
        for (int i=0; i<user_rant; i++) {
            Cluster A_copy=A_init;
            // random rotation
            random_rotate(A_copy);
            double sc_min=A_copy.score(B[iopt]);
            
            // minimizations
            if(user_translate) {       
                for (int imn=0; imn < para::nmn; imn++) {
                    // local rotation
                    double sc = local_rotate(A_copy,B[iopt],sc_min);
                    // translation
                    sc=translate(A_copy,B[iopt], sc);
                    // converge 
                    sc_min=sc;
                    if(2*std::abs(sc-sc_min) < para::converge*(std::abs(sc)+std::abs(sc_min)+para::LEPSO))
                        break;
                }

            }
            else {
                double sc=local_rotate(A_copy, B[iopt], sc_min);
                sc_min = sc ;
            }
            if( sc_min < sc_ran ){ 
                sc_ran = sc_min ;
                A_ran_min = A_copy;
                iconv=0;
            }
            else
                iconv += 1;
            //if(iconv >= nran_conv)
                //break;
        }
        if( sc_ran < sc_f) {
            sc_f = sc_ran;
            best_opt =  iopt;
            A_f = A_ran_min;
        }
    }
    delete[] B;
    A_f.comment(words+" "+std::to_string((long double)sc_f)+" opt  "+std::to_string((long long)best_opt));
    if(user_position) { 
        std::ofstream fout;
        fout.open(filename.c_str());
        A_f.writefile(fout);
        fout.close();
    }
    return sc_f;
}

double local_rotate(Cluster &A, const Cluster & B, double sc) {
    double degree = para::rdeg;
    int count = 0;
    double vfor[3],vtor[3];
    double sc_now=A.paired_score(B,vfor,vtor);
    for (int j=0; j < para::nlo; j++) {
        count += 1;
        Cluster A_backup=A;
        A.rotate(vtor[0],vtor[1],vtor[2],degree);
        double cfor[3],ctor[3];
        sc_now=A.paired_score(B,cfor,ctor);
        if(2*std::abs(sc_now-sc) < para::converge*(std::abs(sc)+std::abs(sc_now)+para::LEPSO) ){
            sc = sc_now;
            break;
        }
        else {
            if(sc_now < sc) {
                sc = sc_now;
                for (int m=0; m<3; m++)
                    vtor[m]=ctor[m];
            }
            else {
                degree = degree *para::ratio; 
                A = A_backup;
            }
        }
    }
    return sc;
}

void random_rotate(Cluster & A) {
    double ux = 0.5 - ((double) rand() / (RAND_MAX)),
           uy = 0.5 - ((double) rand() / (RAND_MAX)),
           uz = 0.5 - ((double) rand() / (RAND_MAX)),
           degree = 360 * ((double) rand() / (RAND_MAX));
    A.rotate(ux,uy,uz,degree);
    return;
}

double translate(Cluster & A, const Cluster& B, double sc) {
    double vfor[3], vtor[3];
    double sc_now=A.paired_score(B,vfor,vtor);
    double step=para::LSTEPSTART;
    for (int j=0; j<para::LMAX; j++) {
        Cluster A_backup = A;
        A.translate(vfor[0],vfor[1],vfor[2],step);
        double cfor[3],ctor[3];
        sc_now=A.paired_score(B,cfor,ctor);
        if( 2*std::abs(sc_now-sc) < para::converge*(std::abs(sc)+std::abs(sc_now)+para::LEPSO)){
            sc = sc_now;
            break;
        }
        else {
            if( sc_now < sc ) {
                sc = sc_now;
                for (int m=0; m<3; m++)
                    vfor[m]=cfor[m];
            }
            else {
                step = step * para::ratio;
                A = A_backup;         
            }
        }
    }
    return sc;
}
