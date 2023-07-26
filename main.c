// Author : Joe McNease
// Date   : 07/2023
// Contact: joemcnease11@gmail.com
//
//
// Finite-Difference Time-Domain (FDTD) solver for the isotropic elastic (P-SV) wave
// equation. Optional Convolutional Perfectly Matched Layers (C-PML) implemented
// at the boundaries to simulate infinite media.
//
// Code based on classical velocity-stress staggered grid formulation developed
// by Madriaga (1976) and Vireaux (1984, 1986).
//
// 
// C-PML formulation is from Martin and Komatitsch et al. (2008), i.e. the
// unsplit C-PML improved at grazing incidence. You can find the original
// C-PML codes that Komatitsch mostly wrote at his github:
//      https://github.com/komatits/seismic_cpml/
//
// In particular, I used this code as a reference:
//      https://github.com/komatits/seismic_cpml/blob/master/
//          seismic_CPML_2D_isotropic_second_order.f90
//
// Use it however you would like, but please cite one (or all) of Komatitsch,
// Roden, Martin, or Gedney's papers on the topic if you use it.
//
//
// To be very clear, below is a reference diagram for the staggered grid in
// space and time.
//                           ^
//                          /
//                t(k+1/2) +---------------------->
//                        /|          |
//                       / |          |
//                      /  |          |
//                     /   |          |
//                    /    |(i,j+1/2) |(i+1/2,j+1/2)
//                   /     +----------+----------->
//                  /      |        vz|
//                 /       |          |
//                /        |          |
//               /         |          |
//              /   vx,rho |(i,j)     |(i+1/2,j)
//             / lambda,mu +----------+----------->
//            /            |          |
//           /             |          |
//          /              V          V
//         /
//   t(k) +---------------------->
//        |          |
//        |          |
//        |          |
//        |          |
//        |(i,j+1/2) |(i+1/2,j+1/2)
//    txx +----------+----------->
//        |          |
//        |          |
//        |          |
//        |          |
//        |(i,j)     |(i+1/2,j)
//        +----------+-----------> X
//        |      txx |
//        |      tzz |
//        V          V
//        Z
//
//
// Since centered finite-differences are used, the formulation used here is
// second-order accurate in space and time. The spatial accuracy can easily be
// increased by adding higher order terms to the difference operators, such as
// in Levander's (1988) fourth-order accurate in space formulation.
//
//


#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>


void fd2d_pml(float *p, const int nx, const int nz, const float dx, const float dz,
              const int nt, const float dt, const float *stf, const int sx, const int sz,
              const float *c, bool save);

void print_array(const float *array, const int size);

void save_pressure(const float *p, const int nx, const int nz, const char *filename);

void fill(float *array, const int size, const float val);

void print_array(const float *array, const int size) {
    for (int i=0; i<size; i++) {
        printf("array[%i] = %10.5f \n", i, array[i]);
    }
}


void save_pressure(const float *p, const int nx, const int nz, const char *filename) {
    FILE *fptr;
    fptr = fopen(filename, "w");

    for (int i=0; i<nz-1; ++i) {
        for (int j=0; j<nx-1; ++j) {
            fprintf(fptr, "%6.3f,", p[i*nx + j]);
        }
        fputs("\n", fptr);
    }
    fclose(fptr);
}


void fill(float *array, const int size, const float val) {
    for (int i=0; i<size-1; i++) {
        array[i] = val;
    }
}


void fd2d_pml(float *p, const int nx, const int nz, const float dx, const float dz,
          const int nt, const float dt, const float *stf, const int sx, const int sz,
          const float *c, bool save) {

    bool USE_PML_LEFT = true;
    bool USE_PML_RIGHT = true;
    bool USE_PML_TOP = false;
    bool USE_PML_BOTTOM = true;

    bool USE_FREE_SURFACE_BC = true;

    int NPOINTS_PML = 30;
    int NPOINTS_PML_X = NPOINTS_PML;
    int NPOINTS_PML_Z = NPOINTS_PML;

    float f0 = 100;
    float NPOWER = 2.0;
    float K_MAX_PML = 1.0;
    float ALPHA_MAX_PML = 2.0*M_PI*(f0/2);

    // These should be passed to function
    float cp = 3000;
    float cs = cp / 1.732;
    float density = 2800;

    float vx[nx*nz];
    float vz[nx*nz];
    float txx[nx*nz];
    float tzz[nx*nz];
    float txz[nx*nz];
    float lambda[nx*nz];
    float rho[nx*nz];
    float mu[nx*nz];
    float pressure[nx*nz];

    // memory variables
    // 
    // TODO: We need to change this so that the memory variables
    // are only stored in the PML regions. The current implementation
    // wastes a ton of memory... (nx*nz)*8*4 bytes currently when it
    // should only require NPOINTS_PML*4*4 (if all sides used and symmetric)
    // bytes. This needs to be changed. It is not difficult, just some extra
    // control flow statements in the time stepping section.
    float memory_dvx_dx[nx*nz];
    float memory_dvx_dz[nx*nz];
    float memory_dvz_dx[nx*nz];
    float memory_dvz_dz[nx*nz];
    float memory_dtxx_dx[nx*nz];
    float memory_dtzz_dz[nx*nz];
    float memory_dtxz_dx[nx*nz];
    float memory_dtxz_dz[nx*nz];

    // damping profile
    float d_x[nx];
    float K_x[nx];
    float alpha_x[nx];
    float a_x[nx];
    float b_x[nx];
    float d_x_half[nx];
    float K_x_half[nx];
    float alpha_x_half[nx];
    float a_x_half[nx];
    float b_x_half[nx];

    float d_z[nz];
    float K_z[nz];
    float alpha_z[nz];
    float a_z[nz];
    float b_z[nz];
    float d_z_half[nz];
    float K_z_half[nz];
    float alpha_z_half[nz];
    float a_z_half[nz];
    float b_z_half[nz];

    float thickness_PML_x = NPOINTS_PML * dx;
    float thickness_PML_z = NPOINTS_PML * dz;
    float Rcoef = 0.001;

    float d0_x = -(NPOWER + 1)*cp*logf(Rcoef) / (2.*thickness_PML_x);
    float d0_z = -(NPOWER + 1)*cp*logf(Rcoef) / (2.*thickness_PML_z);

    for (int i=0; i<nz-1; i++) {
        d_z[i] = 0.;
        d_z_half[i] = 0.;
        K_z[i] = 1.;
        K_z_half[i] = 1.;
        alpha_z[i] = 0.;
        alpha_z_half[i] = 0.;
        a_z[i] = 0.;
        a_z_half[i] = 0.;
    }
    for (int i=0; i<nx-1; i++) {
        d_x[i] = 0.;
        d_x_half[i] = 0.;
        K_x[i] = 1.;
        K_x_half[i] = 1.;
        alpha_x[i] = 0.;
        alpha_x_half[i] = 0.;
        a_x[i] = 0.;
        a_x_half[i] = 0.;
    }

    // damping in the x direction
    float xoriginleft = thickness_PML_x;
    float xoriginright = (nx-1)*dx - thickness_PML_x;

    // Calculate damping profile
    for (int i=0; i<nx-1; i++) {
        float xval = i*dx;

        if (USE_PML_LEFT) {
            float abscissa_in_PML = xoriginleft - xval;
            if (abscissa_in_PML >= 0) {
                float abscissa_normalized = abscissa_in_PML / thickness_PML_x;
                d_x[i] = d0_x * powf(abscissa_normalized, NPOWER);
                K_x[i] = 1.0 + (K_MAX_PML - 1.0) * powf(abscissa_normalized, NPOWER);
                alpha_x[i] = ALPHA_MAX_PML * (1.0 - abscissa_normalized);
            }

            abscissa_in_PML = xoriginleft - (xval + dx/2.0);
            if (abscissa_in_PML >= 0) {
                float abscissa_normalized = abscissa_in_PML / thickness_PML_x;
                d_x_half[i] = d0_x * powf(abscissa_normalized, NPOWER);
                K_x_half[i] = 1.0 + (K_MAX_PML - 1.0) * powf(abscissa_normalized, NPOWER);
                alpha_x_half[i] = ALPHA_MAX_PML * (1.0 - abscissa_normalized);
            }
        }

        if (USE_PML_RIGHT) {
            float abscissa_in_PML = xval - xoriginright;
            if (abscissa_in_PML >= 0) {
                float abscissa_normalized = abscissa_in_PML / thickness_PML_x;
                d_x[i] = d0_x * powf(abscissa_normalized, NPOWER);
                K_x[i] = 1.0 + (K_MAX_PML - 1.0) * powf(abscissa_normalized, NPOWER);
                alpha_x[i] = ALPHA_MAX_PML * (1.0 - abscissa_normalized);
            }

            abscissa_in_PML = (xval + dx/2.0) - xoriginright;
            if (abscissa_in_PML >= 0) {
                float abscissa_normalized = abscissa_in_PML / thickness_PML_x;
                d_x_half[i] = d0_x * powf(abscissa_normalized, NPOWER);
                K_x_half[i] = 1.0 + (K_MAX_PML - 1.0) * powf(abscissa_normalized, NPOWER);
                alpha_x_half[i] = ALPHA_MAX_PML * (1.0 - abscissa_normalized);
            }
        }

        if (alpha_x[i] < 0) { alpha_x[i] = 0.; }
        if (alpha_x_half[i] < 0) { alpha_x_half[i] = 0.; }

        b_x[i] = expf(-(d_x[i]/K_x[i] + alpha_x[i]) * dt);
        b_x_half[i] = expf(-(d_x_half[i]/K_x_half[i] + alpha_x_half[i]) * dt);

        if (fabs(d_x[i]) > 1.e-6) {
            a_x[i] = d_x[i] * (b_x[i] - 1.0) / (K_x[i] * (d_x[i] + K_x[i] * alpha_x[i]));
        }
        if (fabs(d_x_half[i]) > 1.e-6) {
            a_x_half[i] = d_x_half[i] * (b_x_half[i] - 1.0) / (K_x_half[i] * (d_x_half[i] + K_x_half[i] * alpha_x_half[i]));
        }
    }

    float zoriginbottom = thickness_PML_z;
    float zorigintop = (nz-1)*dz - thickness_PML_z;

    for (int j=0; j<nz-1; j++) {
        float zval = j*dz;
        if (USE_PML_TOP) {
            float abscissa_in_PML = zoriginbottom - zval;
            if (abscissa_in_PML >= 0) {
                float abscissa_normalized = abscissa_in_PML / thickness_PML_z;
                d_z[j] = d0_z * powf(abscissa_normalized, NPOWER);
                K_z[j] = 1.0 + (K_MAX_PML - 1.0) * powf(abscissa_normalized, NPOWER);
                alpha_z[j] = ALPHA_MAX_PML * (1.0 - abscissa_normalized);
            }

            abscissa_in_PML = zoriginbottom - (zval + dz/2.);
            if (abscissa_in_PML >= 0) {
                float abscissa_normalized = abscissa_in_PML / thickness_PML_z;
                d_z_half[j] = d0_z * powf(abscissa_normalized, NPOWER);
                K_z_half[j] = 1.0 + (K_MAX_PML - 1.0) * powf(abscissa_normalized, NPOWER);
                alpha_z_half[j] = ALPHA_MAX_PML * (1.0 - abscissa_normalized);
            }
        }

        if (USE_PML_BOTTOM) {
            float abscissa_in_PML = zval - zorigintop;
            if (abscissa_in_PML >= 0) {
                float abscissa_normalized = abscissa_in_PML / thickness_PML_z;
                d_z[j] = d0_z * powf(abscissa_normalized, NPOWER);
                K_z[j] = 1.0 + (K_MAX_PML - 1.0) * powf(abscissa_normalized, NPOWER);
                alpha_z[j] = ALPHA_MAX_PML * (1.0 - abscissa_normalized);
            }

            abscissa_in_PML = (zval + dz/2.) - zorigintop;
            if (abscissa_in_PML >= 0) {
                float abscissa_normalized = abscissa_in_PML / thickness_PML_z;
                d_z_half[j] = d0_z * powf(abscissa_normalized, NPOWER);
                K_z_half[j] = 1.0 + (K_MAX_PML - 1.0) * powf(abscissa_normalized, NPOWER);
                alpha_z_half[j] = ALPHA_MAX_PML * (1.0 - abscissa_normalized);
            }
        }

        b_z[j] = expf(-(d_z[j] / K_z[j] + alpha_z[j]) * dt);
        b_z_half[j] = expf(-(d_z_half[j] / K_z_half[j] + alpha_z_half[j]) * dt);

        if (fabs(d_z[j]) > 1.e-6) {
            a_z[j] = d_z[j] * (b_z[j] - 1.0) / (K_z[j] * (d_z[j] + K_z[j]*alpha_z[j]));
        }
        if (fabs(d_z_half[j]) > 1.e-6) {
            a_z_half[j] = d_z_half[j] * (b_z_half[j] - 1.0) / (K_z_half[j] * (d_z_half[j] + K_z_half[j]*alpha_z_half[j]));
        }
    }

    for (int i=0; i<nz-1; i++) {
        for (int j=0; j<nx-1; j++) {
            int idx = i*nx + j;

            rho[idx] = density;
            mu[idx] = density*cs*cs; // Make scalar with mu = 0;
            lambda[idx] = density*(cp*cp - 2.*cs*cs);
        }
    }
    
    // Check Courant stability condition
    float courant_number = cp * dt * sqrtf(1.0/(dx*dx) + 1.0/(dz*dz));
    printf("Courant Number = %f \n", courant_number);
    if (courant_number > 1.) {
        printf("Time step is too large, simulation is unstable. \n");
        exit(1);
    }

    for (int i=0; i<nz-1; i++) {
        for (int j=0; j<nx-1; j++) {
            int idx = i*nx + j;

            vx[idx] = 0.;
            vz[idx] = 0.;
            txx[idx] = 0.;
            tzz[idx] = 0.;
            txz[idx] = 0.;
            pressure[idx] = 0.;

            memory_dvx_dx[idx] = 0.;
            memory_dvx_dz[idx] = 0.;
            memory_dvz_dx[idx] = 0.;
            memory_dvz_dz[idx] = 0.;
            memory_dtxx_dx[idx] = 0.;
            memory_dtzz_dz[idx] = 0.;
            memory_dtxz_dx[idx] = 0.;
            memory_dtxz_dz[idx] = 0.;
        }
    }

    for (int it=0; it<nt-1; it++) {

        // Compute stress (txx, tzz)
        for (int i=1; i<nz-1; i++) {
            for (int j=0; j<nx-2; j++) {
                int idx = i*nx + j;

                float lambda_half_x = 0.5 * (lambda[idx+1] + lambda[idx]);
                float mu_half_x = 0.5 * (mu[idx+1] + mu[idx]);

                float lambda_plus_2mu_half_x = lambda_half_x + 2.0*mu_half_x;

                float value_dvx_dx = (vx[idx+1] - vx[idx]) / dx;
                float value_dvz_dz = (vz[idx] - vz[idx-nx]) / dz;

                memory_dvx_dx[idx] = b_x_half[j]*memory_dvx_dx[idx] + a_x_half[j]*value_dvx_dx;
                memory_dvz_dz[idx] = b_z[i]*memory_dvz_dz[idx] + a_z[i]*value_dvz_dz;

                value_dvx_dx = value_dvx_dx / K_x_half[j] + memory_dvx_dx[idx];
                value_dvz_dz = value_dvz_dz / K_z[i] + memory_dvz_dz[idx];

                txx[idx] = txx[idx] + (lambda_plus_2mu_half_x*value_dvx_dx + lambda_half_x*value_dvz_dz)*dt;
                tzz[idx] = tzz[idx] + (lambda_half_x*value_dvx_dx + lambda_plus_2mu_half_x*value_dvz_dz)*dt;
            }
        }

        // Compute stress (txz)
        for (int i=0; i<nz-2; i++) {
            for (int j=1; j<nx-1; j++) {
                int idx = i*nx + j;

                float mu_half_z = 0.5 * (mu[idx+nx] + mu[idx]);

                float value_dvz_dx = (vz[idx] - vz[idx-1]) / dx;
                float value_dvx_dz = (vx[idx+nx] - vx[idx]) / dz;

                memory_dvz_dx[idx] = b_x[j]*memory_dvz_dx[idx] + a_x[j]*value_dvz_dx;
                memory_dvx_dz[idx] = b_z_half[i]*memory_dvx_dz[idx] + a_z_half[i]*value_dvx_dz;

                value_dvz_dx = value_dvz_dx / K_x[j] + memory_dvz_dx[idx];
                value_dvx_dz = value_dvx_dz / K_z_half[i] + memory_dvx_dz[idx];

                // If in PML region, add convolutional term
                if (i <= NPOINTS_PML_Z) {
                }
                if (i >= ((nz-1)-NPOINTS_PML_Z)) {
                }
                if (j <= NPOINTS_PML_X) {
                }
                if (j >= ((nx-1)-NPOINTS_PML_X)) {
                }

                txz[idx] = txz[idx] + mu_half_z*(value_dvz_dx + value_dvx_dz)*dt;
            }
        }

        for (int i=1; i<nz-1; i++) {
            for (int j=1; j<nx-1; j++) {
                int idx = i*nx + j;

                float value_dsigma_xx_dx = (txx[idx] - txx[idx-1]) / dx;
                float value_dsigma_xz_dz = (txz[idx] - txz[idx-nx]) / dz;

                memory_dtxx_dx[idx] = b_x[j]*memory_dtxx_dx[idx] + a_x[j]*value_dsigma_xx_dx;
                memory_dtxz_dz[idx] = b_z[i]*memory_dtxz_dz[idx] + a_z[i]*value_dsigma_xz_dz;

                value_dsigma_xx_dx = value_dsigma_xx_dx / K_x[j] + memory_dtxx_dx[idx];
                value_dsigma_xz_dz = value_dsigma_xz_dz / K_z[i] + memory_dtxz_dz[idx];

                vx[idx] = vx[idx] + (value_dsigma_xx_dx + value_dsigma_xz_dz)*dt/rho[idx];
            }
        }

        for (int i=0; i<nz-2; i++) {
            for (int j=0; j<nx-2; j++) {
                int idx = i*nx + j;

                float rho_half_x_half_z = 0.25 * (rho[idx] + rho[idx+1] + rho[idx+nx+1] + rho[idx+nx]);

                float value_dsigma_xz_dx = (txz[idx+1] - txz[idx]) / dx;
                float value_dsigma_zz_dz = (tzz[idx+nx] - tzz[idx]) / dz;

                memory_dtxz_dx[idx] = b_x_half[j]*memory_dtxz_dx[idx] + a_x_half[j]*value_dsigma_xz_dx;
                memory_dtzz_dz[idx] = b_z_half[i]*memory_dtzz_dz[idx] + a_z_half[i]*value_dsigma_zz_dz;

                value_dsigma_xz_dx = value_dsigma_xz_dx / K_x_half[j] + memory_dtxz_dx[idx];
                value_dsigma_zz_dz = value_dsigma_zz_dz / K_z_half[i] + memory_dtzz_dz[idx];

                vz[idx] = vz[idx] + (value_dsigma_xz_dx + value_dsigma_zz_dz)*dt/rho_half_x_half_z;
            }
        }

        // Add source term (directional force)
        // float ANGLE_FORCE = 0;
        // float a = M_PI*M_PI*f0*f0;
        // float t = it*dt;
        // float t0 = 0.01;
        // float factor = 1.e7;
        // float source_term = -factor * 2.0*a*(t-t0)*expf(-a*powf((t-t0), 2));
        // float force_x = cosf(ANGLE_FORCE*(M_PI/180)) * source_term;
        // float force_z = sinf(ANGLE_FORCE*(M_PI/180)) * source_term;
        // int sidx = sz*nx + sx;
        // 
        // float rho_half_x_half_z = 0.25 * (rho[sidx] + rho[sidx+nx] + rho[sidx+nx+1] + rho[sidx+1]);

        // vx[sidx] = vx[sidx] + force_x*dt / rho[sidx];
        // vz[sidx] = vz[sidx] + force_z*dt / rho_half_x_half_z;

        float a = M_PI*M_PI*f0*f0;
        float t = it*dt;
        float t0 = 0.01;
        float factor = 1.e6;
        float source_term = -factor * 2.0*a*(t-t0)*expf(-a*powf((t-t0), 2));
        int sidx = sz*nx + sx;
        float rho_half_x_half_z = 0.25 * (rho[sidx] + rho[sidx+nx] + rho[sidx+nx+1] + rho[sidx+1]);
        txx[sidx] = txx[sidx] + source_term;
        tzz[sidx] = tzz[sidx] + source_term;


        // Calculate pressure
        for (int i=0; i<nz-1; i++) {
            for (int j=0; j<nx-1; j++) {
                int idx = i*nx + j;
                pressure[idx] = 0.5 * (tzz[idx] + txx[idx]);
            }
        }

        for (int i=0; i<nz-1; i++) {
            int idx = i*nx;
            vx[idx] = 0.;
            vz[idx] = 0.;
            vx[idx+nx-1] = 0.;
            vz[idx+nx-1] = 0.;
        }
        for (int i=0; i<nx-1; i++) {
            int idx = i;
            if (USE_FREE_SURFACE_BC) {
                txz[idx] = 0.;
                tzz[idx] = 0.;
            }
            else {
                vx[idx] = 0.;
                vz[idx] = 0.;
            }
            vx[(nz-1)*nx + idx] = 0.;
            vz[(nz-1)*nx + idx] = 0.;
        }

        if (save && (it%10==0)) {
            char fn[50];
            // Save velocity
            // Currently testing pressure
            sprintf(fn, "vz/vz%08d", it);
            // save_pressure(vz, nx, nz, fn);
            save_pressure(pressure, nx, nz, fn);

            // sprintf(fn, "vx/vx%08d", it);
            // save_pressure(vx, nx, nz, fn);
        }

        printf("Time step: %d", it);
        printf("\n");
    }
}


int main() {

#define nx 301
#define nz 201
#define dx 1.0
#define dz 1.0
#define sx floor(nx/2)
#define sz floor(nz/2)
#define nt 1000
#define dt 0.0001
#define t0 0.01
#define f0 100.0
#define c0 3000.0


    printf("Starting program \n");

    int res = mkdir("vz", 0777);

    float c[nx*nz];
    float p[nx*nz];
    float stf[nt];

    fill(c, nx*nz, c0);
    fill(p, nx*nz, 0.);
    fill(stf, nt, 0.);

    for (int i=0; i<nt-1; i++) {
        stf[i] = exp(-(f0*f0)*(i*dt-t0)*(i*dt-t0));
    }

    clock_t start, end;
    double elapsed;

    start = clock();
    fd2d_pml(p, nx, nz, dx, dz, nt, dt, stf, sx, sz, c, true);
    end = clock();

    elapsed = (double)(end-start) / CLOCKS_PER_SEC;
    printf("Time to compute: %f \n", elapsed);

    return 0;
}
