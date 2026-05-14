///#-- ========================================================================================
//#--               Modelling and Simulation Project #13
//#-- ========================================================================================
// Particles in confinement: between two parallel plates 
//An interesting and important avenue to alter phase behavior is to change the boundary conditions by, 
// e.g. confining particles between walls aka plates.  
// In fact, the exact structures formed can depend on the distance between the plates.   
 
//Goal: Examine what happens when hard spheres are confined between two plates of different separation.  An intriguing additional exploration is to see what happens when they are not parallel. 

#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Simulation parameters */
const int mc_steps = 10000;
const int output_steps = 100;

const double diameter = 1.0;
const double delta = 0.001;

/* Choose initial configuration */
const char* init_filename =
    "square/trilat_H1.5_eta0.40.dat";

/* Variables */
int n_particles = 0;

double radius;

double r[N][NDIM];
double box[NDIM];

/* Read initial configuration */
void read_data(void){

    FILE *fp = fopen(init_filename, "r");

    if (fp == NULL){
        printf("Error: could not open file.\n");
        return;
    }

    fscanf(fp, "%d", &n_particles);

    for (int d = 0; d < NDIM; ++d) {

        double min, max;

        fscanf(fp, "%lf %lf", &min, &max);

        box[d] = max;
    }

    for (int n = 0; n < n_particles; ++n) {

        for (int d = 0; d < NDIM; ++d) {

            fscanf(fp, "%lf", &r[n][d]);
        }

        double skip;

        fscanf(fp, "%lf", &skip);
    }

    fclose(fp);

    printf("Particles: %d\n", n_particles);

    printf("Box: %lf %lf %lf\n",
           box[0],
           box[1],
           box[2]);
}

/* Monte Carlo move */
int move_particle(void){

    int n =
        (int)(dsfmt_genrand() * n_particles);

    double newpos[NDIM];

    /* random displacement */
    for(int d = 0; d < NDIM; d++){

        double displacement =
            (2.0 * dsfmt_genrand() - 1.0)
            * delta;

        newpos[d] =
            r[n][d] + displacement;
    }

    /* periodic boundaries in x,y */
    for(int d = 0; d < 2; d++){

        if(newpos[d] < 0.0)
            newpos[d] += box[d];

        else if(newpos[d] >= box[d])
            newpos[d] -= box[d];
    }

    /* hard walls in z */
    if(newpos[2] < radius ||
       newpos[2] > box[2] - radius){

        return 0;
    }

    /* overlap check */
    for(int m = 0; m < n_particles; m++){

        if(m == n)
            continue;

        double dist = 0.0;

        for(int d = 0; d < NDIM; d++){

            double dx =
                newpos[d] - r[m][d];

            /* nearest image convention */
            if(d < 2){

                if(dx > 0.5 * box[d])
                    dx -= box[d];

                if(dx < -0.5 * box[d])
                    dx += box[d];
            }

            dist += dx * dx;
        }

        if(dist < diameter * diameter){

            return 0;
        }
    }

    /* accept move */
    for(int d = 0; d < NDIM; d++){

        r[n][d] = newpos[d];
    }

    return 1;
}

/* Write output */
void write_data(int step){

    char buffer[128];

    sprintf(buffer,
            "coords_step%07d.dat",
            step);

    FILE* fp = fopen(buffer, "w");

    fprintf(fp, "%d\n", n_particles);

    for(int d = 0; d < NDIM; ++d){

        fprintf(fp,
                "%lf %lf\n",
                0.0,
                box[d]);
    }

    for(int n = 0; n < n_particles; ++n){

        for(int d = 0; d < NDIM; ++d){

            fprintf(fp,
                    "%f\t",
                    r[n][d]);
        }

        fprintf(fp,
                "%lf\n",
                diameter / 2.0); //writing radius for Ovito file
    }

    fclose(fp);
}

int main(void){

    radius = 0.5 * diameter;

    read_data();

    if(n_particles == 0){

        printf("Error: no particles.\n");

        return 0;
    }

    dsfmt_seed(time(NULL));

    int accepted = 0;

    /* Monte Carlo loop */
    for(int step = 0;
        step < mc_steps;
        ++step){

        for(int n = 0;
            n < n_particles;
            ++n){

            accepted += move_particle();
        }

        if(step % output_steps == 0){

            printf(
                "Step %d | Acceptance: %lf\n",
                step,
                (double)accepted /
                (n_particles * output_steps)
            );

            accepted = 0;

            write_data(step);
        }
    }

    return 0;
}