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
#include <direct.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

#define N_BINS 200
const double H_list[] = {2.0}; //{1.0, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
const double eta_list[] = {0.40}; //{0.30, 0.35, 0.40, 0.45, 0.50, 0.55};

/* Initialization variables */
const int mc_steps = 50000; // production steps per density
const int eq_steps = 50000; // Added equilibration steps before sampling g(r)
const int output_steps = 1000;
const int sample_every = 10; // Sample g(r) every this many steps
const double diameter = 1.0;

/* Simulation variables */
int n_particles;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];
double box_min[NDIM], box_max[NDIM]; // Added these to store box limits globally

/* Compression/Expansion stage */
const int ENABLE_COMPRESSION = 1; // Set to 1 to enable compression stage
const double eta_target = 0.50;   // Target packing fraction for compression stage
const int compression_steps = 25000; // Number of steps for compression stage
const int post_compression_eq_steps = 50000; // Additional equilibration steps after compression

int wall_hit_count = 0;

/* Functions */
void read_data(const char *filename){
    FILE *fp = fopen(filename, "r");
    assert(fp != NULL);

    // Read number of particles
    assert(fscanf(fp, "%d", &n_particles) == 1);

    assert(n_particles <= N);

    // Read box limits (min max for each dimension)
    for(int d = 0; d < NDIM; d++){
        assert(fscanf(fp, "%lf %lf", &box_min[d], &box_max[d]) == 2);
        box[d] = box_max[d] - box_min[d];
    }

    // Read particle coordinates
    for(int i = 0; i < n_particles; i++){

        double dummy;  // fourth column (ignored)

        assert(fscanf(fp, "%lf %lf %lf %lf",
                      &r[i][0],
                      &r[i][1],
                      &r[i][2],
                      &dummy) == 4);
    }

    fclose(fp);

    // debug check
    printf("Read %d particles\n", n_particles);

    for(int d = 0; d < NDIM; d++)
        printf("Box[%d] = %f\n", d, box[d]);

    for(int i = 0; i < 3; i++){
        printf("Particle %d: %f %f %f\n",
               i, r[i][0], r[i][1], r[i][2]);
    }
}


int move_particle(double delta){
    // Pick random particle
    int i = (int)(dsfmt_genrand() * n_particles);

    // Store old position
    double old_pos[NDIM];
    for (int d = 0; d < NDIM; d++) {
        old_pos[d] = r[i][d];
    }

    // Trial move
    // x and y: periodic boundaries
    for (int d = 0; d < 2; d++) {
        r[i][d] += (2.0 * dsfmt_genrand() - 1.0) * delta;;

        // Periodic boundaries
        if (r[i][d] < box_min[d]) {
            r[i][d] += box[d];
        }
        if (r[i][d] >= box_max[d]) {
            r[i][d] -= box[d];
        }
    }
    // z: hard walls at box_min[2] and box_max[2]
    r[i][2] += (2.0 * dsfmt_genrand() - 1.0) * delta;
    if (r[i][2] < box_min[2] + radius || r[i][2] > box_max[2] - radius) {
        // Reject move if wall is hit
        for (int d = 0; d < NDIM; d++) {
            r[i][d] = old_pos[d];
        }
        wall_hit_count++;
        return 0;
    }

    // Overlap checking
    for (int j = 0; j < n_particles; j++) {

        if (i == j) continue; // Do not check distance to the same particle we are moving

        double dist2 = 0.0;

        for (int d = 0; d < 2; d++) { // Check in x and y
            double dx = r[i][d] - r[j][d];

            // Ensure periodic boundaries with nearest image convention
            if (dx > 0.5 * box[d]) dx -= box[d];
            if (dx < -0.5 * box[d]) dx += box[d];

            dist2 += dx * dx; // Add distance in dimension to total distance
        }

        double dz = r[i][2] - r[j][2]; // Check in z without periodic boundaries
        dist2 += dz * dz;

        // If distance < diameter -> overlap
        if (dist2 < diameter * diameter - 1e-12) {
            // Move is rejected: dimension values are returned to old position
            for (int d = 0; d < NDIM; d++) {
                r[i][d] = old_pos[d];
            }

            return 0;
        }
    }

    // No overlap -> accept move
    return 1;
}

void write_data(double H, double eta, int step) {
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "run2/data/coords/square_steps/H%.1f_eta%.2f/coords_step%07d.dat", H, eta, step);            // ============ Make sure file pathing is correct for run ================
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n\n", n_particles);
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter/2.0);
    }
    fclose(fp);
}

void write_data_compress(double H, double eta, int step) {
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "run2/data/compress/square_steps/H%.1f_eta%.2f-%.2f/compress_step%07d.dat", H, eta, eta_target, step);            // ============ Make sure file pathing is correct for run ================
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n\n", n_particles);
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter/2.0);
    }
    fclose(fp);
}

// Scaling box for compression/expansion stage to reach target packing fraction eta_target
void scale_box_step(double scale_xy) {
    box[0] *= scale_xy; box_max[0] = box_min[0] + box[0];
    box[1] *= scale_xy; box_max[1] = box_min[1] + box[1];
    for (int i = 0; i < n_particles; i++) {
        r[i][0] = box_min[0] + (r[i][0] - box_min[0]) * scale_xy;
        r[i][1] = box_min[1] + (r[i][1] - box_min[1]) * scale_xy;
    }
}

// g(r) functions

/* Loops over all unique pairs (i < j), computes seperation using minimum image convention
Use r_max = box_min / 2 to never cross minimum image boundaries
hist[b] counts pairs whose distance falls in [b*dr, (b+1)*dr]
Add 2 because its a pair, representing both directions*/

// In-plane g(r) for confined system (only x and y dimentsions)
void accumulate_gr_inplane(double *hist, double r_max, double dr) {
    for (int i = 0; i < n_particles - 1; i++) {
        for (int j = i + 1; j < n_particles; j++) {
            double dx = r[i][0] - r[j][0];
            double dy = r[i][1] - r[j][1];
            // Minimum image convention for x and y
            if (dx > 0.5 * box[0]) dx -= box[0];
            if (dx < -0.5 * box[0]) dx += box[0];
            if (dy > 0.5 * box[1]) dy -= box[1];
            if (dy < -0.5 * box[1]) dy += box[1];
            double dist = sqrt(dx * dx + dy * dy);
            if (dist < r_max) {
                int b = (int)(dist / dr);
                if (b < N_BINS) hist[b] += 2.0; // Each pair contributes to g(r) twice (i,j and j,i)
            }
        }
    }
}


/* Normalises histogram and writes it to file
Normalisation (Eq. 8.2):
    g(r + 0.5*dr) = hist[b] / (N * n_samples * n_ideal(b))
    
2D ideal gas:
    n_ideal(b) = rho2d * pi * ((r_high)^2 - (r_low^2))
    where r_high = (b+1)*dr and r_low = b*dr
*/
void write_gr_inplane(const char *fname, double *hist, long n_samples, double rho2d, double r_max, double dr) {
    FILE *fp = fopen(fname, "w");
    assert(fp != NULL);
    fprintf(fp, "# r/sigma g(r) (rho2d=%.4f, n=%d, samples=%ld)\n", rho2d, n_particles, n_samples);

    double norm = (double)n_particles * (double)n_samples;

    for (int b = 0; b < N_BINS; b++) {
        double r_low = b * dr;
        double r_high = r_low + dr;
        double r_mid = r_low + 0.5 * dr;
        // 2D ideal gas
        double n_ideal = rho2d * M_PI * (r_high * r_high - r_low * r_low); // Area of annulus
        double gr = (n_ideal > 0.0) ? hist[b] / (norm * n_ideal) : 0.0; // Avoid division by zero
        fprintf(fp, "%.6f %.6f\n", r_mid, gr);
    }
    fclose(fp);
    printf("g(r) written to %s\n", fname);
}

// Density profile rho(z)
#define N_ZBINS 200
void accumulate_rho_z(double *hist, double dz) {
    for (int i = 0; i < n_particles; i++) {
        int b = (int)((r[i][2] - box_min[2]) / dz);
        if (b >= 0 && b < N_ZBINS) hist[b] += 1.0;
    }
}

void write_rhoz(const char *fname, double *zhist, long n_samples, double H, double dz, double area) {
    FILE *fp = fopen(fname, "w");
    assert(fp != NULL);
    fprintf(fp, "# z\trho(z)\n");
    for (int b = 0; b < N_ZBINS; b++) {
        double z_mid = box_min[2] + (b + 0.5) * dz;
        double rhoz = zhist[b] / (n_samples * area * dz);
        fprintf(fp, "%.6f %.6f\n", z_mid, rhoz);
    }
    fclose(fp);
    printf("Density profile rho(z) written to %s\n", fname);
}

// Main
int main(int argc, char* argv[]){

    assert(diameter > 0.0);

    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    dsfmt_seed(time(NULL));

    // Loop over H and eta
    for (int H_idx = 0; H_idx < sizeof(H_list)/sizeof(H_list[0]); H_idx++) {
        for (int eta_idx = 0; eta_idx < sizeof(eta_list)/sizeof(eta_list[0]); eta_idx++){

            double H = H_list[H_idx];
            double eta = eta_list[eta_idx];
            printf("\n===| H = %.2f, eta = %.2f |===\n", H, eta);

            // Set init_frame based on current H and eta
            char init_frame[64];
            snprintf(init_frame, sizeof(init_frame), "square/sqlat_H%.1f_eta%.2f.dat", H, eta);            // ============ Make sure file pathing is correct for run ================

            read_data(init_frame); // Re-read initial configuration for each density to start from the same state

            if(n_particles == 0){
            printf("Error: Number of particles, n_particles = 0.\n");
            return 0;
            }

            // Tune delta to get reasonable acceptance (0.3 - 0.5)
            double delta_local = 0.05 * diameter;
            assert(delta_local > 0.0);

            // number density rho = N / V
            double area = box[0] * box[1];
            double rho2d = n_particles / area;
            printf("rho2d = %.6f, area = %.4f\n", rho2d, area);

            // r_max = half the shortest box side
            double r_max = fmin(box[0], box[1]) * 0.5;
            double dr = r_max / N_BINS;

            // zero histogram
            double hist[N_BINS];
            for(int b = 0; b < N_BINS; b++) hist[b] = 0.0;

            // zero z histogram
            double zhist[N_ZBINS];
            for(int b = 0; b < N_ZBINS; b++) zhist[b] = 0.0;
            double dz = H / N_ZBINS;

            // equilibration
            int accepted = 0;
            for(int step = 0; step < eq_steps; step++) {
                for(int n = 0; n < n_particles; n++){
                    accepted += move_particle(delta_local);
                }

                if ((step + 1) % 1000 == 0) {
                    double rate = (double)accepted / (1000.0 * n_particles);
                    if (rate < 0.20) delta_local *= 0.9;
                    if (rate > 0.40) delta_local *= 1.1;
                    delta_local = fmax(1e-4, fmin(delta_local, 0.5 * diameter));
                    accepted = 0;
                }
            }
            printf("Tuned delta = %.4f, final acceptance = %.2f%%\n", delta_local, 100.0 * accepted / (double)(1000.0 * n_particles));

            // Compression/Expansion stage
            if (ENABLE_COMPRESSION && eta_target != eta) {
                double scale_total = sqrt(eta / eta_target); // Total scaling factor to reach target packing fraction
                double scale_per_step = pow(scale_total, 1.0 / compression_steps); // Scaling factor per step
                accepted = 0;

                char folder_name[128];
                snprintf(folder_name, sizeof(folder_name), "run2/data/compress/square_steps/H%.1f_eta%.2f-%.2f", H, eta, eta_target);           // ============ Make sure file pathing is correct for run ================
                if (_mkdir(folder_name) == 0) {
                    printf("Folder created successfully!\n");
                } else {
                    printf("Unable to create folder. (It may already exist)\n");
                }

                for (int step = 0; step < compression_steps; step++) {
                    scale_box_step(scale_per_step);

                    for (int n = 0; n < n_particles; n++) {
                        accepted += move_particle(delta_local);
                    }

                    if ((step + 1) % 1000 == 0) {
                        double rate = (double)accepted / (1000.0 * n_particles);
                        if (rate < 0.20) delta_local *= 0.9;
                        if (rate > 0.40) delta_local *= 1.1;
                        delta_local = fmax(1e-4, fmin(delta_local, 0.5 * diameter));
                        accepted = 0;
                    }

                    if(step % output_steps == 0){
                        printf("Step %d. Move acceptance = %.2f%%\tWall hits = %d\n", step, 100.0 * accepted / (double)(output_steps * n_particles), wall_hit_count);
                        write_data_compress(H, eta, step);
                        accepted = 0;
                        wall_hit_count = 0;
                    }
                }
                //Recompute all box-derived quantities for production stage
                area   = box[0] * box[1];
                rho2d  = n_particles / area;
                r_max  = fmin(box[0], box[1]) * 0.5;
                dr     = r_max / N_BINS;
                printf("Compression/Expansion done. Final eta = %.4f, rho2d = %.6f\n", n_particles * particle_volume / (area * H), rho2d);

                // Post-compression/expansion equilibration
                accepted = 0;
                for (int step = 0; step < post_compression_eq_steps; step++) {
                    for (int n = 0; n < n_particles; n++) {
                        accepted += move_particle(delta_local);
                    }

                    if ((step + 1) % 1000 == 0) {
                        double rate = (double)accepted / (1000.0 * n_particles);
                        if (rate < 0.20) delta_local *= 0.9;
                        if (rate > 0.40) delta_local *= 1.1;
                        delta_local = fmax(1e-4, fmin(delta_local, 0.5 * diameter));
                        accepted = 0;
                    }
                }
                printf("Post-compression equilibration done. Tuned delta = %.4f, final acceptance = %.2f%%\n", delta_local, 100.0 * accepted / (double)(1000.0 * n_particles));
            }


            // production + g(r) sampling
            accepted = 0;
            long n_samples = 0;

            char folder_name[128];
            snprintf(folder_name, sizeof(folder_name), "run2/data/coords/square_steps/H%.1f_eta%.2f", H, (ENABLE_COMPRESSION ? eta_target : eta));           // ============ Make sure file pathing is correct for run ================
            if (_mkdir(folder_name) == 0) {
                printf("Folder created successfully!\n");
            } else {
                printf("Unable to create folder. (It may already exist)\n");
            }

            for(int step = 0; step < mc_steps; step++){
                for(int n = 0; n < n_particles; n++){
                    accepted += move_particle(delta_local);
                }

                if(step % sample_every == 0){
                    accumulate_gr_inplane(hist, r_max, dr);
                    accumulate_rho_z(zhist, dz);
                    n_samples++;
                }

                if(step % output_steps == 0){
                    printf("Step %d. Move acceptance = %.2f%%\tWall hits = %d\n", step, 100.0 * accepted / (double)(output_steps * n_particles), wall_hit_count);
                    write_data(H, eta, step);
                    accepted = 0;
                    wall_hit_count = 0;
                }
            }

            printf("Production done. Samples = %ld\n", n_samples);

            // Write g(r) for this density
            char fname[64];
            snprintf(fname, sizeof(fname), "run2/data/gr/square_gr/gr_H%.2f_eta%.2f.dat", H, (ENABLE_COMPRESSION ? eta_target : eta));               // ============ Make sure file pathing is correct for run ================
            write_gr_inplane(fname, hist, n_samples, rho2d, r_max, dr);
            snprintf(fname, sizeof(fname), "run2/data/rhoz/square_rhoz/rhoz_H%.2f_eta%.2f.dat", H, (ENABLE_COMPRESSION ? eta_target : eta));         // ============ Make sure file pathing is correct for run ================
            write_rhoz(fname, zhist, n_samples, H, dz, area);
        }
    }

    printf("\nDone.\n");
    return 0;
}
