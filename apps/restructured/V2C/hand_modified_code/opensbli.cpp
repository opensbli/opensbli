#include <stdlib.h> 
#include <string.h> 
#include <math.h>
#include <complex.h>
double Delta0block0;
double Delta1block0;
int block0np0;
int block0np1;
double rkold[3];
double rknew[3];
double inv_0;
double inv_1;
double dt;
int niter;
double rc2;
double rc3;

#define OPS_2D
#include "ops_seq.h"
#include "opensbliblock00_kernels.h"

const int nx_sb = 16;
double hamming_window[nx_sb];

int main(int argc, char **argv) 
{
Delta0block0=0.00325221;
Delta1block0=0.00561589;
block0np0=1999;
block0np1=1348;

rkold[0]=1.0/4.0;
rkold[1]=3.0/20.0;
rkold[2]=3.0/5.0;
rknew[0]=2.0/3.0;
rknew[1]=5.0/12.0;
rknew[2]=3.0/5.0;
inv_0=1.0/Delta0block0;
inv_1=1.0/Delta1block0;
dt=1;
niter=1;
rc2=2.0/3.0;
rc3=1.0/12.0;

for (int j=0; j<nx_sb; j++){
  hamming_window[j] = (0.54 - 0.46*cos(2*M_PI*j/(nx_sb)))/0.54;
}

// Initializing OPS 
ops_init(argc,argv,1);
ops_decl_const("Delta0block0" , 1, "double", &Delta0block0);
ops_decl_const("Delta1block0" , 1, "double", &Delta1block0);
ops_decl_const("block0np0" , 1, "int", &block0np0);
ops_decl_const("block0np1" , 1, "int", &block0np1);
ops_decl_const("inv_0" , 1, "double", &inv_0);
ops_decl_const("inv_1" , 1, "double", &inv_1);
ops_decl_const("dt" , 1, "double", &dt);
ops_decl_const("niter" , 1, "int", &niter);
ops_decl_const("rc2" , 1, "double", &rc2);
ops_decl_const("rc3" , 1, "double", &rc3);

// Define and Declare OPS Block
ops_block opensbliblock00 = ops_decl_block(2, "opensbliblock00");

ops_reduction N2sum = ops_decl_reduction_handle(sizeof(double), "double", "reduction_N2sum");
ops_reduction N4sum_r = ops_decl_reduction_handle(sizeof(double), "double", "reduction_N4sum_r");
ops_reduction N4sum_i = ops_decl_reduction_handle(sizeof(double), "double", "reduction_N4sum_i");
ops_reduction N8sum_r = ops_decl_reduction_handle(sizeof(double), "double", "reduction_N8sum_r");
ops_reduction N8sum_i = ops_decl_reduction_handle(sizeof(double), "double", "reduction_N8sum_i");
ops_reduction sbx0 = ops_decl_reduction_handle(sizeof(double), "double", "reduction_sbx0");
ops_reduction sbx1 = ops_decl_reduction_handle(sizeof(double), "double", "reduction_sbx1");

#include "defdec_data_set.h"
// Define and declare stencils
int stencil_0_00temp[] = {0, 0};
ops_stencil stencil_0_00 = ops_decl_stencil(2,1,stencil_0_00temp,"stencil_0_00temp");
int stencil_0_02temp[] = {0, -2, 0, -1, 0, 1, 0, 2};
ops_stencil stencil_0_02 = ops_decl_stencil(2,4,stencil_0_02temp,"stencil_0_02temp");
int stencil_0_01temp[] = {-2, 0, -1, 0, 1, 0, 2, 0};
ops_stencil stencil_0_01 = ops_decl_stencil(2,4,stencil_0_01temp,"stencil_0_01temp");
#include "bc_exchanges.h"
// Init OPS partition
ops_partition("");

int iteration_range_26[] = {-2, block0np0 + 2, -2, block0np1 + 2};
ops_par_loop(opensbliblock00Kernel026, "Grid_based_initialisation0", opensbliblock00, 2, iteration_range_26,
ops_arg_dat(rhou1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(velocity0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(velocity1_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(vorticity0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(vorticity1_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_19[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel019, "MetricsEquation CD x1_B0 xi1 ", opensbliblock00, 2, iteration_range_19,
ops_arg_dat(x1_B0, 1, stencil_0_02, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_20[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel020, "MetricsEquation CD x0_B0 xi1 ", opensbliblock00, 2, iteration_range_20,
ops_arg_dat(x0_B0, 1, stencil_0_02, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_21[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel021, "MetricsEquation CD x1_B0 xi0 ", opensbliblock00, 2, iteration_range_21,
ops_arg_dat(x1_B0, 1, stencil_0_01, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_22[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel022, "MetricsEquation CD x0_B0 xi0 ", opensbliblock00, 2, iteration_range_22,
ops_arg_dat(x0_B0, 1, stencil_0_01, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_25[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel025, "MetricsEquation evaluation", opensbliblock00, 2, iteration_range_25,
ops_arg_dat(wk3_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(D00_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(D01_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(D11_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(D10_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(detJ_B0, 1, stencil_0_00, "double", OPS_WRITE));


double cpu_start0, elapsed_start0;
ops_timers(&cpu_start0, &elapsed_start0);

//ops_halo_transfer(exchange15);
//ops_halo_transfer(exchange16);
//ops_halo_transfer(exchange17);
//ops_halo_transfer(exchange18);
int iteration_range_30[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel030, "Save equations", opensbliblock00, 2, iteration_range_30,
ops_arg_dat(vorticity0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(vorticity1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(vorticity1_old_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(vorticity0_old_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_5[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel005, "Zeroing residuals", opensbliblock00, 2, iteration_range_5,
ops_arg_dat(Residual0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(Residual1_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_8[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel008, "Convective CD velocity0_B0 xi0 ", opensbliblock00, 2, iteration_range_8,
ops_arg_dat(velocity0_B0, 1, stencil_0_01, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_10[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel010, "Convective CD velocity1_B0 xi0 ", opensbliblock00, 2, iteration_range_10,
ops_arg_dat(velocity1_B0, 1, stencil_0_01, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_12[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel012, "Convective CD velocity0_B0 xi1 ", opensbliblock00, 2, iteration_range_12,
ops_arg_dat(velocity0_B0, 1, stencil_0_02, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_13[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel013, "Convective CD velocity1_B0 xi1 ", opensbliblock00, 2, iteration_range_13,
ops_arg_dat(velocity1_B0, 1, stencil_0_02, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_14[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel014, "Convective residual ", opensbliblock00, 2, iteration_range_14,
ops_arg_dat(D01_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(D00_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(D11_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(D10_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(Residual0_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(Residual1_B0, 1, stencil_0_00, "double", OPS_RW));


int iteration_range_32[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel032, "Sub stage advancement", opensbliblock00, 2, iteration_range_32,
ops_arg_dat(vorticity1_old_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(vorticity0_old_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(Residual0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(Residual1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(vorticity0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(vorticity1_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_gbl(&rknew[0], 1, "double", OPS_READ));


int iteration_range_31[] = {0, block0np0, 0, block0np1};
ops_par_loop(opensbliblock00Kernel031, "Temporal solution advancement", opensbliblock00, 2, iteration_range_31,
ops_arg_dat(Residual0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(Residual1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(vorticity1_old_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(vorticity0_old_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_gbl(&rkold[0], 1, "double", OPS_READ));

// Error indicator.

// For more details see: C. T. Jacobs, M. Zauner, N. De Tullio, Satya P. Jammy, David J. Lusher, N. D. Sandham (Submitted). "An error indicator for finite difference methods using spectral techniques with application to aerofoil simulation". Submitted to the ParCFD 2017 Special Issue of Computers & Fluids.

// For a 16^3 subset of a DNS (1) work out A2,A4,A8 in each direction, averaged over 16^2 lines (2) repeat in each direction. Gives a directional error indicator (can also do for each Q). Free stream uniform flow may need attention (do we need to subtract uniform flow from this first?).

// Let r=test exponent (i.e. we'd like the solution to go down faster than k^r - depends on the variable we are looking at).
double r = -0.5;
double A2, A4, A8;
double A2temp, A4temp, A8temp;

double sbx0_sum, sbx1_sum; // The approximate mid-point location of each sub-block.

int iter_l[4];
int start_idx;

// Loop over each nx_sb^2 subblock.
for (int sbi=0; sbi<2*(block0np0/nx_sb)-1; sbi++)
{
    for (int sbj=0; sbj<2*(block0np1/nx_sb)-1; sbj++)
    {
        sbx0_sum = 0;
        sbx1_sum = 0;
        ops_printf("Sub-block: (%d,%d)\n", sbi, sbj);
        for (int ldir=0; ldir<2; ldir++)
        {
            // Reset the maximums.
            A2 = 0.0;
            A4 = 0.0;
            A8 = 0.0;
            
            

            // Loop over nx_sb lines.
            for (int li=0; li<nx_sb; li++)
            {
                if (ldir == 0)
                {
                    // x lines
                    iter_l[0] = sbi*(0.5*nx_sb);
                    iter_l[1] = sbi*(0.5*nx_sb) + nx_sb;
                    iter_l[2] = sbj*(0.5*nx_sb) + li;
                    iter_l[3] = sbj*(0.5*nx_sb) + (li+1);
                    start_idx = sbi*(0.5*nx_sb);
                }
                else
                {
                    if (ldir == 1)
                    {
                        // y lines
                        iter_l[0] = sbi*(0.5*nx_sb) + li;
                        iter_l[1] = sbi*(0.5*nx_sb) + (li+1);
                        iter_l[2] = sbj*(0.5*nx_sb);
                        iter_l[3] = sbj*(0.5*nx_sb) + nx_sb;
                        start_idx = sbj*(0.5*nx_sb);
                    }
                }
                ops_par_loop(opensbliblock00KernelErrorIndicator, "Error indicator", opensbliblock00, 2, iter_l,
                                        ops_arg_dat(vorticity0_B0, 1, stencil_0_00, "double", OPS_READ),
                                        ops_arg_dat(x0_B0, 1, stencil_0_00, "double", OPS_READ),
                                        ops_arg_dat(x1_B0, 1, stencil_0_00, "double", OPS_READ),
                                        ops_arg_gbl(hamming_window, 1, "double", OPS_READ),
                                        ops_arg_gbl(&ldir, 1, "int", OPS_READ),
                                        ops_arg_gbl(&start_idx, 1, "int", OPS_READ),
                                        ops_arg_idx(),
                                        ops_arg_reduce(N2sum, 1, "double", OPS_INC),
                                        ops_arg_reduce(N4sum_r, 1, "double", OPS_INC),
                                        ops_arg_reduce(N4sum_i, 1, "double", OPS_INC),
                                        ops_arg_reduce(N8sum_r, 1, "double", OPS_INC),
                                        ops_arg_reduce(N8sum_i, 1, "double", OPS_INC),
                                        ops_arg_reduce(sbx0, 1, "double", OPS_INC),
                                        ops_arg_reduce(sbx1, 1, "double", OPS_INC));

                double N2sum_reduction = 0.0; 
                double N4sum_r_reduction = 0.0; 
                double N4sum_i_reduction = 0.0; 
                double N8sum_r_reduction = 0.0; 
                double N8sum_i_reduction = 0.0; 
                ops_reduction_result(N2sum, &N2sum_reduction);
                ops_reduction_result(N4sum_r, &N4sum_r_reduction);
                ops_reduction_result(N4sum_i, &N4sum_i_reduction);
                ops_reduction_result(N8sum_r, &N8sum_r_reduction);
                ops_reduction_result(N8sum_i, &N8sum_i_reduction);
                //ops_printf("N2sum: %lf\n", N2sum_reduction);
                //ops_printf("N4sum_r: %lf\n", N4sum_r_reduction);
                //ops_printf("N4sum_i: %lf\n", N4sum_i_reduction);
                //ops_printf("N8sum_r: %lf\n", N8sum_r_reduction);
                //ops_printf("N8sum_i: %lf\n", N8sum_i_reduction);

                double N2sum_final = N2sum_reduction;
                double __complex__ N4sum_final = N4sum_r_reduction + N4sum_i_reduction*I;
                double __complex__ N8sum_final = N8sum_r_reduction + N8sum_i_reduction*I;

                A2temp = fabs(N2sum_final/nx_sb)*pow(2, -2*r);
                if(A2temp > A2)
                {
                    A2 = A2temp;
                }
                A4temp = sqrt(pow(__real__ N4sum_final*2/nx_sb, 2) + pow(__imag__ N4sum_final*2/nx_sb, 2))*pow(2, -r);
                if(A4temp > A4)
                {
                    A4 = A4temp;
                }
                A8temp = sqrt(pow(__real__ N8sum_final*2/nx_sb, 2) + pow(__imag__ N8sum_final*2/nx_sb, 2));
                if(A8temp > A8)
                {
                    A8 = A8temp;
                }
                
                double sbx0_reduction = 0.0; 
                double sbx1_reduction = 0.0; 
                ops_reduction_result(sbx0, &sbx0_reduction);
                ops_reduction_result(sbx1, &sbx1_reduction);
                
                sbx0_sum += sbx0_reduction;
                sbx1_sum += sbx1_reduction;
                
            }
            
            
                    
            double eps = 3.0E-2; // Avoids problems in uniform flow
            double ferrest = log(1 + (int)(A2/(A4+eps)) + (int)(A4/(A8+eps)) + (int)(A2/(A8+eps))); // Expressed as a float. (high=bad, 0=good)
            int ierrest = (int)(A2>(A4+eps)) + (int)(A4>(A8+eps)) + (int)(A2>(A8+eps)); // An integer 0, 1, 2, or 3 (3=worst, 0=best)
            ops_printf("ferrest for block(%d,%d) in direction(%d): %lf\n", sbi, sbj, ldir, ferrest);
            ops_printf("ierrest for block(%d,%d) in direction(%d): %d\n", sbi, sbj, ldir, ierrest);

            
        }
        sbx0_sum /= (nx_sb*2);
            sbx1_sum /= (nx_sb*2);
        ops_printf("sbx0: %f\n", sbx0_sum);
            ops_printf("sbx1: %f\n", sbx1_sum);

    }
    
}

//ops_halo_transfer(exchange15);
//ops_halo_transfer(exchange16);
//ops_halo_transfer(exchange17);
//ops_halo_transfer(exchange18);

double cpu_end0, elapsed_end0;
ops_timers(&cpu_end0, &elapsed_end0);
ops_printf("\nTimings are:\n");
ops_printf("-----------------------------------------\n");
ops_printf("Total Wall time %lf\n",elapsed_end0-elapsed_start0);
ops_exit();
//Main program end 
}
