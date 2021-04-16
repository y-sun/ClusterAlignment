#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <ctime>
#include <stdlib.h>
#include <fstream>
#include "cluster.h"
#include "io.h"


double align(Cluster & A, const Cluster &B, string tag, string words);
void random_rotate(Cluster &A);
double translate(Cluster & A, const Cluster &B, double sc) ;
double local_rotate(Cluster & A, const Cluster &B, double sc);


#endif
