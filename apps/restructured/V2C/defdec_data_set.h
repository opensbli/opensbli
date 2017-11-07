ops_dat vorticity1_old_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
vorticity1_old_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "vorticity1_old_B0");
}
ops_dat velocity1_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
velocity1_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "velocity1_B0");
}
ops_dat vorticity0_old_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
vorticity0_old_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "vorticity0_old_B0");
}
ops_dat wk1_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
wk1_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "wk1_B0");
}
ops_dat rho_B0;
{
rho_B0 = ops_decl_dat_hdf5(opensbliblock00, 1, "double", "rho_B0", "data.h5");
}
ops_dat wk4_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
wk4_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "wk4_B0");
}
ops_dat vorticity0_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
vorticity0_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "vorticity0_B0");
}
ops_dat D01_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
D01_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "D01_B0");
}
ops_dat D10_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
D10_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "D10_B0");
}
ops_dat Residual1_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
Residual1_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "Residual1_B0");
}
ops_dat rhou1_B0;
{
rhou1_B0 = ops_decl_dat_hdf5(opensbliblock00, 1, "double", "rhou1_B0", "data.h5");
}
ops_dat wk0_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
wk0_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "wk0_B0");
}
ops_dat Residual0_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
Residual0_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "Residual0_B0");
}
ops_dat vorticity1_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
vorticity1_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "vorticity1_B0");
}
ops_dat D00_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
D00_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "D00_B0");
}
ops_dat wk3_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
wk3_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "wk3_B0");
}
ops_dat rhou0_B0;
{
rhou0_B0 = ops_decl_dat_hdf5(opensbliblock00, 1, "double", "rhou0_B0", "data.h5");
}
ops_dat detJ_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
detJ_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "detJ_B0");
}
ops_dat D11_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
D11_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "D11_B0");
}
ops_dat x1_B0;
{
x1_B0 = ops_decl_dat_hdf5(opensbliblock00, 1, "double", "x1_B0", "data.h5");
}
ops_dat velocity0_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
velocity0_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "velocity0_B0");
}
ops_dat x0_B0;
{
x0_B0 = ops_decl_dat_hdf5(opensbliblock00, 1, "double", "x0_B0", "data.h5");
}
ops_dat wk2_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
wk2_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "wk2_B0");
}
ops_dat wk5_B0;
{
int halo_p[] = {5, 5};
int halo_m[] = {-5, -5};
int size[] = {block0np0, block0np1};
int base[] = {0, 0};
double* value = NULL;
wk5_B0 = ops_decl_dat(opensbliblock00, 1, size, base, halo_m, halo_p, value, "double", "wk5_B0");
}