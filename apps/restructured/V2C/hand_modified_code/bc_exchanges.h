// Boundary condition exchange code on opensbliblock00 direction 0 left
ops_halo_group exchange15 ;
{
int halo_iter[] = {2, block0np1 + 4};
int from_base[] = {0, -2};
int to_base[] = {block0np0, -2};
int dir[] = {1, 2};
ops_halo halo0 = ops_decl_halo(vorticity0_B0, vorticity0_B0, halo_iter, from_base, to_base, dir, dir);
ops_halo halo1 = ops_decl_halo(vorticity1_B0, vorticity1_B0, halo_iter, from_base, to_base, dir, dir);
ops_halo grp[] = {halo0,halo1};
exchange15 = ops_decl_halo_group(2,grp);
}
// Boundary condition exchange code on opensbliblock00 direction 0 right
ops_halo_group exchange16 ;
{
int halo_iter[] = {2, block0np1 + 4};
int from_base[] = {block0np0 - 2, -2};
int to_base[] = {-2, -2};
int dir[] = {1, 2};
ops_halo halo0 = ops_decl_halo(vorticity0_B0, vorticity0_B0, halo_iter, from_base, to_base, dir, dir);
ops_halo halo1 = ops_decl_halo(vorticity1_B0, vorticity1_B0, halo_iter, from_base, to_base, dir, dir);
ops_halo grp[] = {halo0,halo1};
exchange16 = ops_decl_halo_group(2,grp);
}
// Boundary condition exchange code on opensbliblock00 direction 1 left
ops_halo_group exchange17 ;
{
int halo_iter[] = {block0np0 + 4, 2};
int from_base[] = {-2, 0};
int to_base[] = {-2, block0np1};
int dir[] = {1, 2};
ops_halo halo0 = ops_decl_halo(vorticity0_B0, vorticity0_B0, halo_iter, from_base, to_base, dir, dir);
ops_halo halo1 = ops_decl_halo(vorticity1_B0, vorticity1_B0, halo_iter, from_base, to_base, dir, dir);
ops_halo grp[] = {halo0,halo1};
exchange17 = ops_decl_halo_group(2,grp);
}
// Boundary condition exchange code on opensbliblock00 direction 1 right
ops_halo_group exchange18 ;
{
int halo_iter[] = {block0np0 + 4, 2};
int from_base[] = {-2, block0np1 - 2};
int to_base[] = {-2, -2};
int dir[] = {1, 2};
ops_halo halo0 = ops_decl_halo(vorticity0_B0, vorticity0_B0, halo_iter, from_base, to_base, dir, dir);
ops_halo halo1 = ops_decl_halo(vorticity1_B0, vorticity1_B0, halo_iter, from_base, to_base, dir, dir);
ops_halo grp[] = {halo0,halo1};
exchange18 = ops_decl_halo_group(2,grp);
}