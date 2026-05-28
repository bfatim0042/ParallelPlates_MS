#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double H_list[] = {1.1, 1.2, 1.3};//{1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
const double eta_list[] = {0.30, 0.35, 0.40, 0.45, 0.50, 0.55};

int main(void) {
    for (int H_idx = 0; H_idx < sizeof(H_list)/sizeof(H_list[0]); H_idx++) {
        for (int eta_idx = 0; eta_idx < sizeof(eta_list)/sizeof(eta_list[0]); eta_idx++) {
            // Parameters
            int n_layers = (int)floor(H_list[H_idx]);   // Number of layers in z direction
            int Nx = 10;        // Number of particles along x
            int Ny = 10;        // Number of particles along y
            const double diameter = 1.0;
            double H = H_list[H_idx];    // Target box height
            double eta = eta_list[eta_idx]; // Target packing fraction

            int N_per_layer = Nx * Ny;
            int N = N_per_layer * n_layers; // Total number of particles

            // Solve for H from eta = N * V_sphere / (Lx * Ly * H)
            double V_sphere = (4.0/3.0) * M_PI * pow(diameter/2.0, 3);

            double A_needed = (N * V_sphere) / (eta * H);
            double ax = sqrt(A_needed / (Nx * Ny)) + 1e-10;
            if (ax < diameter) {
                ax = diameter + 1e-10;
            }
            double ay = ax;

            // Box dimensions
            double Lx = Nx * ax;
            double Ly = Ny * ay;

            double actual_eta = (N * V_sphere) / (Lx * Ly * H);

            // Layer z-positions
            double z_low = diameter / 2.0; // First layer at z = radius
            double z_high = H - diameter / 2.0; // Last layer at z = H - radius
            double z_spacing = (n_layers > 1) ? (z_high - z_low) / (n_layers - 1) : 0.0; // Spacing between layers when more than one layer

            // Check if H is large enough to fit n_layers without overlap
            if (n_layers > 1 && z_spacing < diameter) {
                printf("Error: H = %.4f is too small for %d layers at eta = %.2f.\n", H, n_layers, eta);
                printf("Minimum H needed: %.4f\n", diameter * (n_layers - 1) + diameter);
                return 1;
            }

            // Print parameters
            printf("N = %d particles, %d layers of %d\n", N, n_layers, N_per_layer);
            printf("Box: Lx = %.4f, Ly = %.4f, H = %.4f\n", Lx, Ly, H);
            //printf("ax = %.4f, diameter = %.4f\n", ax, diameter);
            if (n_layers > 1) {
                printf("Layer z-spacing: %.4f (min allowed: %.4f)\n", z_spacing, diameter);
            } else {
                printf("Single layer at z = %.4f\n", H / 2.0);
            }
            printf("Packing fraction eta = %.4f (target: %.4f)\n", actual_eta, eta);

            // Placing particles
            double x[N], y[N], z[N];
            int idx = 0;
            for (int layer = 0; layer < n_layers; layer++) {
                double z_pos = (n_layers == 1) ? H / 2.0 : z_low + layer * z_spacing; // H/2 for single layer, otherwise spaced layers
                double x_offset = 0.0;
                double y_offset = 0.0;

                for (int j = 0; j < Ny; j++) {
                    double row_shift = 0.0;
                    for (int i = 0; i < Nx; i++) {
                        x[idx] = i * ax + x_offset + row_shift;
                        y[idx] = j * ay + y_offset;
                        z[idx] = z_pos;

                        // Wrap around box boundaries [0, L)
                        x[idx] = fmod(x[idx], Lx);
                        y[idx] = fmod(y[idx], Ly);
                        if (x[idx] < 0) x[idx] += Lx;
                        if (y[idx] < 0) y[idx] += Ly;

                        idx++;
                    }
                }
            }

            // Overlap check
            int overlaps = 0;
            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {
                    double dx = x[i] - x[j];
                    double dy = y[i] - y[j];
                    double dz = z[i] - z[j];

                    // Apply periodic boundaries in x and y
                    if (dx >  0.5 * Lx) dx -= Lx;
                    if (dx < -0.5 * Lx) dx += Lx;
                    if (dy >  0.5 * Ly) dy -= Ly;
                    if (dy < -0.5 * Ly) dy += Ly;

                    double dist2 = dx*dx + dy*dy + dz*dz;
                    if (dist2 < diameter * diameter - 1e-10) overlaps++;
                }
            }
            printf("Overlap check: %d overlaps found.\n", overlaps);
            if (overlaps > 0) {
                return 1;
            }

            //Write to file
            char filename[256];
            sprintf(filename, "square/sqlat_H%.1f_eta%.2f.dat", H, eta);
            FILE *fp = fopen(filename, "w");
            if (!fp) {
                printf("Error opening file for writing.\n");
                return 1;
            }

            fprintf(fp, "%d\n", N);
            //fprintf(fp, "\n"); // Blank line for Ovito format
            fprintf(fp, "%f %f\n", 0.0, Lx); // x: [0, Lx] // Remove these for Ovito format
            fprintf(fp, "%f %f\n", 0.0, Ly); // y: [0, Ly]
            fprintf(fp, "%f %f\n", 0.0, H);  // z: [0, H]

            for (int i = 0; i < N; i++) {
                fprintf(fp, "%f %f %f %f\n", x[i], y[i], z[i], diameter/2.0);
            }

            fclose(fp);
            printf("Written to %s\n\n", filename);
        }
    }
    return 0;
}