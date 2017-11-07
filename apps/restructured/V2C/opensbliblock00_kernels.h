#ifndef OPENSBLIBLOCK00_KERNEL_H
#define OPENSBLIBLOCK00_KERNEL_H
 void opensbliblock00Kernel026(const double *rhou1_B0,const double *rhou0_B0,const double *rho_B0,double
*velocity0_B0,double *velocity1_B0,double *vorticity0_B0,double *vorticity1_B0)
{
   if(rho_B0[OPS_ACC2(0,0)] != 0.0)
   {
       velocity0_B0[OPS_ACC3(0,0)] = rhou0_B0[OPS_ACC1(0,0)]/rho_B0[OPS_ACC2(0,0)];
       velocity1_B0[OPS_ACC4(0,0)] = rhou1_B0[OPS_ACC0(0,0)]/rho_B0[OPS_ACC2(0,0)];
   }
   else
   {
       velocity0_B0[OPS_ACC3(0,0)] = 0.0;
       velocity1_B0[OPS_ACC4(0,0)] = 0.0;
   }
   
   vorticity0_B0[OPS_ACC5(0,0)] = 0;
   vorticity1_B0[OPS_ACC6(0,0)] = 0;

}

void opensbliblock00Kernel019(const double *x1_B0,double *wk0_B0, const int *idx)
{

wk0_B0[OPS_ACC1(0,0)] = ((idx[1] == 0) ? (
   (3.00000000000002*x1_B0[OPS_ACC0(0,1)] - 1.83333333333334*x1_B0[OPS_ACC0(0,0)] + 1.06910884386911e-15*x1_B0[OPS_ACC0(0,5)] + 0.333333333333356*x1_B0[OPS_ACC0(0,3)] - 8.34617916606957e-15*x1_B0[OPS_ACC0(0,4)] - 1.50000000000003*x1_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 1) ? (
   (0.719443173328855*x1_B0[OPS_ACC0(0,1)] - 0.322484932882161*x1_B0[OPS_ACC0(0,0)] - 0.376283677513354*x1_B0[OPS_ACC0(0,-1)] - 0.0658051057710389*x1_B0[OPS_ACC0(0,3)] + 0.00571369039775442*x1_B0[OPS_ACC0(0,4)] + 0.0394168524399448*x1_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 2) ? (
   (0.521455851089587*x1_B0[OPS_ACC0(0,1)] + 0.197184333887745*x1_B0[OPS_ACC0(0,0)] - 0.791245592765872*x1_B0[OPS_ACC0(0,-1)] + 0.113446470384241*x1_B0[OPS_ACC0(0,-2)] - 0.00412637789557492*x1_B0[OPS_ACC0(0,3)] - 0.0367146847001262*x1_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 3) ? (
   (0.652141084861241*x1_B0[OPS_ACC0(0,1)] + 0.0451033223343881*x1_B0[OPS_ACC0(0,0)] - 0.727822147724592*x1_B0[OPS_ACC0(0,-1)] + 0.121937153224065*x1_B0[OPS_ACC0(0,-2)] - 0.082033432844602*x1_B0[OPS_ACC0(0,2)] - 0.00932597985049999*x1_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 1) ? (
   -(-1.83333333333334*x1_B0[OPS_ACC0(0,0)] - 1.50000000000003*x1_B0[OPS_ACC0(0,-2)] + 3.00000000000002*x1_B0[OPS_ACC0(0,-1)] + 1.06910884386911e-15*x1_B0[OPS_ACC0(0,-5)] - 8.34617916606957e-15*x1_B0[OPS_ACC0(0,-4)] + 0.333333333333356*x1_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 2) ? (
   -(-0.376283677513354*x1_B0[OPS_ACC0(0,1)] - 0.322484932882161*x1_B0[OPS_ACC0(0,0)] + 0.0394168524399448*x1_B0[OPS_ACC0(0,-2)] + 0.719443173328855*x1_B0[OPS_ACC0(0,-1)] + 0.00571369039775442*x1_B0[OPS_ACC0(0,-4)] - 0.0658051057710389*x1_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 3) ? (
   -(-0.791245592765872*x1_B0[OPS_ACC0(0,1)] + 0.197184333887745*x1_B0[OPS_ACC0(0,0)] - 0.0367146847001262*x1_B0[OPS_ACC0(0,-2)] + 0.521455851089587*x1_B0[OPS_ACC0(0,-1)] + 0.113446470384241*x1_B0[OPS_ACC0(0,2)] - 0.00412637789557492*x1_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 4) ? (
   -(-0.727822147724592*x1_B0[OPS_ACC0(0,1)] + 0.0451033223343881*x1_B0[OPS_ACC0(0,0)] - 0.082033432844602*x1_B0[OPS_ACC0(0,-2)] + 0.652141084861241*x1_B0[OPS_ACC0(0,-1)] - 0.00932597985049999*x1_B0[OPS_ACC0(0,3)] + 0.121937153224065*x1_B0[OPS_ACC0(0,2)])/Delta1block0
)
: (
   (rc2)*x1_B0[OPS_ACC0(0,1)]/Delta1block0 - rc2*x1_B0[OPS_ACC0(0,-1)]/Delta1block0 + (rc3)*x1_B0[OPS_ACC0(0,-2)]/Delta1block0 - rc3*x1_B0[OPS_ACC0(0,2)]/Delta1block0
)))))))));

//    wk0_B0[OPS_ACC1(0,0)] = inv_1*((rc2)*x1_B0[OPS_ACC0(0,1)] + (rc3)*x1_B0[OPS_ACC0(0,-2)] - rc2*x1_B0[OPS_ACC0(0,-1)]
  //    - rc3*x1_B0[OPS_ACC0(0,2)]);

}

void opensbliblock00Kernel020(const double *x0_B0,double *wk1_B0, const int *idx)
{

wk1_B0[OPS_ACC1(0,0)] = ((idx[1] == 0) ? (
   (3.00000000000002*x0_B0[OPS_ACC0(0,1)] - 1.83333333333334*x0_B0[OPS_ACC0(0,0)] + 1.06910884386911e-15*x0_B0[OPS_ACC0(0,5)] + 0.333333333333356*x0_B0[OPS_ACC0(0,3)] - 8.34617916606957e-15*x0_B0[OPS_ACC0(0,4)] - 1.50000000000003*x0_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 1) ? (
   (0.719443173328855*x0_B0[OPS_ACC0(0,1)] - 0.322484932882161*x0_B0[OPS_ACC0(0,0)] - 0.376283677513354*x0_B0[OPS_ACC0(0,-1)] - 0.0658051057710389*x0_B0[OPS_ACC0(0,3)] + 0.00571369039775442*x0_B0[OPS_ACC0(0,4)] + 0.0394168524399448*x0_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 2) ? (
   (0.521455851089587*x0_B0[OPS_ACC0(0,1)] + 0.197184333887745*x0_B0[OPS_ACC0(0,0)] - 0.791245592765872*x0_B0[OPS_ACC0(0,-1)] + 0.113446470384241*x0_B0[OPS_ACC0(0,-2)] - 0.00412637789557492*x0_B0[OPS_ACC0(0,3)] - 0.0367146847001262*x0_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 3) ? (
   (0.652141084861241*x0_B0[OPS_ACC0(0,1)] + 0.0451033223343881*x0_B0[OPS_ACC0(0,0)] - 0.727822147724592*x0_B0[OPS_ACC0(0,-1)] + 0.121937153224065*x0_B0[OPS_ACC0(0,-2)] - 0.082033432844602*x0_B0[OPS_ACC0(0,2)] - 0.00932597985049999*x0_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 1) ? (
   -(-1.83333333333334*x0_B0[OPS_ACC0(0,0)] - 1.50000000000003*x0_B0[OPS_ACC0(0,-2)] + 3.00000000000002*x0_B0[OPS_ACC0(0,-1)] + 1.06910884386911e-15*x0_B0[OPS_ACC0(0,-5)] - 8.34617916606957e-15*x0_B0[OPS_ACC0(0,-4)] + 0.333333333333356*x0_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 2) ? (
   -(-0.376283677513354*x0_B0[OPS_ACC0(0,1)] - 0.322484932882161*x0_B0[OPS_ACC0(0,0)] + 0.0394168524399448*x0_B0[OPS_ACC0(0,-2)] + 0.719443173328855*x0_B0[OPS_ACC0(0,-1)] + 0.00571369039775442*x0_B0[OPS_ACC0(0,-4)] - 0.0658051057710389*x0_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 3) ? (
   -(-0.791245592765872*x0_B0[OPS_ACC0(0,1)] + 0.197184333887745*x0_B0[OPS_ACC0(0,0)] - 0.0367146847001262*x0_B0[OPS_ACC0(0,-2)] + 0.521455851089587*x0_B0[OPS_ACC0(0,-1)] + 0.113446470384241*x0_B0[OPS_ACC0(0,2)] - 0.00412637789557492*x0_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 4) ? (
   -(-0.727822147724592*x0_B0[OPS_ACC0(0,1)] + 0.0451033223343881*x0_B0[OPS_ACC0(0,0)] - 0.082033432844602*x0_B0[OPS_ACC0(0,-2)] + 0.652141084861241*x0_B0[OPS_ACC0(0,-1)] - 0.00932597985049999*x0_B0[OPS_ACC0(0,3)] + 0.121937153224065*x0_B0[OPS_ACC0(0,2)])/Delta1block0
)
: (
   (rc2)*x0_B0[OPS_ACC0(0,1)]/Delta1block0 - rc2*x0_B0[OPS_ACC0(0,-1)]/Delta1block0 + (rc3)*x0_B0[OPS_ACC0(0,-2)]/Delta1block0 - rc3*x0_B0[OPS_ACC0(0,2)]/Delta1block0
)))))))));

//    wk1_B0[OPS_ACC1(0,0)] = inv_1*(-rc3*x0_B0[OPS_ACC0(0,2)] + (rc3)*x0_B0[OPS_ACC0(0,-2)] - rc2*x0_B0[OPS_ACC0(0,-1)] +
  //    (rc2)*x0_B0[OPS_ACC0(0,1)]);

}

void opensbliblock00Kernel021(const double *x1_B0,double *wk2_B0, const int *idx)
{

wk2_B0[OPS_ACC1(0,0)] = ((idx[0] == 0) ? (
   (3.00000000000002*x1_B0[OPS_ACC0(1,0)] - 1.83333333333334*x1_B0[OPS_ACC0(0,0)] + 1.06910884386911e-15*x1_B0[OPS_ACC0(5,0)] + 0.333333333333356*x1_B0[OPS_ACC0(3,0)] - 8.34617916606957e-15*x1_B0[OPS_ACC0(4,0)] - 1.50000000000003*x1_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 1) ? (
   (0.719443173328855*x1_B0[OPS_ACC0(1,0)] - 0.322484932882161*x1_B0[OPS_ACC0(0,0)] - 0.376283677513354*x1_B0[OPS_ACC0(-1,0)] - 0.0658051057710389*x1_B0[OPS_ACC0(3,0)] + 0.00571369039775442*x1_B0[OPS_ACC0(4,0)] + 0.0394168524399448*x1_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 2) ? (
   (0.521455851089587*x1_B0[OPS_ACC0(1,0)] + 0.197184333887745*x1_B0[OPS_ACC0(0,0)] - 0.791245592765872*x1_B0[OPS_ACC0(-1,0)] + 0.113446470384241*x1_B0[OPS_ACC0(-2,0)] - 0.00412637789557492*x1_B0[OPS_ACC0(3,0)] - 0.0367146847001262*x1_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 3) ? (
   (0.652141084861241*x1_B0[OPS_ACC0(1,0)] + 0.0451033223343881*x1_B0[OPS_ACC0(0,0)] - 0.727822147724592*x1_B0[OPS_ACC0(-1,0)] + 0.121937153224065*x1_B0[OPS_ACC0(-2,0)] - 0.082033432844602*x1_B0[OPS_ACC0(2,0)] - 0.00932597985049999*x1_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 1) ? (
   -(-1.83333333333334*x1_B0[OPS_ACC0(0,0)] - 1.50000000000003*x1_B0[OPS_ACC0(-2,0)] + 3.00000000000002*x1_B0[OPS_ACC0(-1,0)] + 1.06910884386911e-15*x1_B0[OPS_ACC0(-5,0)] - 8.34617916606957e-15*x1_B0[OPS_ACC0(-4,0)] + 0.333333333333356*x1_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 2) ? (
   -(-0.376283677513354*x1_B0[OPS_ACC0(1,0)] - 0.322484932882161*x1_B0[OPS_ACC0(0,0)] + 0.0394168524399448*x1_B0[OPS_ACC0(-2,0)] + 0.719443173328855*x1_B0[OPS_ACC0(-1,0)] + 0.00571369039775442*x1_B0[OPS_ACC0(-4,0)] - 0.0658051057710389*x1_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 3) ? (
   -(-0.791245592765872*x1_B0[OPS_ACC0(1,0)] + 0.197184333887745*x1_B0[OPS_ACC0(0,0)] - 0.0367146847001262*x1_B0[OPS_ACC0(-2,0)] + 0.521455851089587*x1_B0[OPS_ACC0(-1,0)] + 0.113446470384241*x1_B0[OPS_ACC0(2,0)] - 0.00412637789557492*x1_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 4) ? (
   -(-0.727822147724592*x1_B0[OPS_ACC0(1,0)] + 0.0451033223343881*x1_B0[OPS_ACC0(0,0)] - 0.082033432844602*x1_B0[OPS_ACC0(-2,0)] + 0.652141084861241*x1_B0[OPS_ACC0(-1,0)] - 0.00932597985049999*x1_B0[OPS_ACC0(3,0)] + 0.121937153224065*x1_B0[OPS_ACC0(2,0)])/Delta0block0
)
: (
   (rc2)*x1_B0[OPS_ACC0(1,0)]/Delta0block0 - rc2*x1_B0[OPS_ACC0(-1,0)]/Delta0block0 + (rc3)*x1_B0[OPS_ACC0(-2,0)]/Delta0block0 - rc3*x1_B0[OPS_ACC0(2,0)]/Delta0block0
)))))))));
//    wk2_B0[OPS_ACC1(0,0)] = inv_0*((rc3)*x1_B0[OPS_ACC0(-2,0)] + (rc2)*x1_B0[OPS_ACC0(1,0)] - rc2*x1_B0[OPS_ACC0(-1,0)]
  //    - rc3*x1_B0[OPS_ACC0(2,0)]);

}

void opensbliblock00Kernel022(const double *x0_B0,double *wk3_B0, const int *idx)
{

wk3_B0[OPS_ACC1(0,0)] = ((idx[0] == 0) ? (
   (3.00000000000002*x0_B0[OPS_ACC0(1,0)] - 1.83333333333334*x0_B0[OPS_ACC0(0,0)] + 1.06910884386911e-15*x0_B0[OPS_ACC0(5,0)] + 0.333333333333356*x0_B0[OPS_ACC0(3,0)] - 8.34617916606957e-15*x0_B0[OPS_ACC0(4,0)] - 1.50000000000003*x0_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 1) ? (
   (0.719443173328855*x0_B0[OPS_ACC0(1,0)] - 0.322484932882161*x0_B0[OPS_ACC0(0,0)] - 0.376283677513354*x0_B0[OPS_ACC0(-1,0)] - 0.0658051057710389*x0_B0[OPS_ACC0(3,0)] + 0.00571369039775442*x0_B0[OPS_ACC0(4,0)] + 0.0394168524399448*x0_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 2) ? (
   (0.521455851089587*x0_B0[OPS_ACC0(1,0)] + 0.197184333887745*x0_B0[OPS_ACC0(0,0)] - 0.791245592765872*x0_B0[OPS_ACC0(-1,0)] + 0.113446470384241*x0_B0[OPS_ACC0(-2,0)] - 0.00412637789557492*x0_B0[OPS_ACC0(3,0)] - 0.0367146847001262*x0_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 3) ? (
   (0.652141084861241*x0_B0[OPS_ACC0(1,0)] + 0.0451033223343881*x0_B0[OPS_ACC0(0,0)] - 0.727822147724592*x0_B0[OPS_ACC0(-1,0)] + 0.121937153224065*x0_B0[OPS_ACC0(-2,0)] - 0.082033432844602*x0_B0[OPS_ACC0(2,0)] - 0.00932597985049999*x0_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 1) ? (
   -(-1.83333333333334*x0_B0[OPS_ACC0(0,0)] - 1.50000000000003*x0_B0[OPS_ACC0(-2,0)] + 3.00000000000002*x0_B0[OPS_ACC0(-1,0)] + 1.06910884386911e-15*x0_B0[OPS_ACC0(-5,0)] - 8.34617916606957e-15*x0_B0[OPS_ACC0(-4,0)] + 0.333333333333356*x0_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 2) ? (
   -(-0.376283677513354*x0_B0[OPS_ACC0(1,0)] - 0.322484932882161*x0_B0[OPS_ACC0(0,0)] + 0.0394168524399448*x0_B0[OPS_ACC0(-2,0)] + 0.719443173328855*x0_B0[OPS_ACC0(-1,0)] + 0.00571369039775442*x0_B0[OPS_ACC0(-4,0)] - 0.0658051057710389*x0_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 3) ? (
   -(-0.791245592765872*x0_B0[OPS_ACC0(1,0)] + 0.197184333887745*x0_B0[OPS_ACC0(0,0)] - 0.0367146847001262*x0_B0[OPS_ACC0(-2,0)] + 0.521455851089587*x0_B0[OPS_ACC0(-1,0)] + 0.113446470384241*x0_B0[OPS_ACC0(2,0)] - 0.00412637789557492*x0_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 4) ? (
   -(-0.727822147724592*x0_B0[OPS_ACC0(1,0)] + 0.0451033223343881*x0_B0[OPS_ACC0(0,0)] - 0.082033432844602*x0_B0[OPS_ACC0(-2,0)] + 0.652141084861241*x0_B0[OPS_ACC0(-1,0)] - 0.00932597985049999*x0_B0[OPS_ACC0(3,0)] + 0.121937153224065*x0_B0[OPS_ACC0(2,0)])/Delta0block0
)
: (
   (rc2)*x0_B0[OPS_ACC0(1,0)]/Delta0block0 - rc2*x0_B0[OPS_ACC0(-1,0)]/Delta0block0 + (rc3)*x0_B0[OPS_ACC0(-2,0)]/Delta0block0 - rc3*x0_B0[OPS_ACC0(2,0)]/Delta0block0
)))))))));

//    wk3_B0[OPS_ACC1(0,0)] = inv_0*(-rc3*x0_B0[OPS_ACC0(2,0)] + (rc3)*x0_B0[OPS_ACC0(-2,0)] - rc2*x0_B0[OPS_ACC0(-1,0)] +
  //    (rc2)*x0_B0[OPS_ACC0(1,0)]);

}

 void opensbliblock00Kernel025(const double *wk3_B0,const double *wk0_B0,const double *wk1_B0,const double
*wk2_B0,double *D00_B0,double *D01_B0,double *D11_B0,double *D10_B0,double *detJ_B0)
{
   detJ_B0[OPS_ACC8(0,0)] = wk0_B0[OPS_ACC1(0,0)]*wk3_B0[OPS_ACC0(0,0)] - wk1_B0[OPS_ACC2(0,0)]*wk2_B0[OPS_ACC3(0,0)];

    D00_B0[OPS_ACC4(0,0)] = wk0_B0[OPS_ACC1(0,0)]/(wk0_B0[OPS_ACC1(0,0)]*wk3_B0[OPS_ACC0(0,0)] -
      wk1_B0[OPS_ACC2(0,0)]*wk2_B0[OPS_ACC3(0,0)]);

    D01_B0[OPS_ACC5(0,0)] = -wk1_B0[OPS_ACC2(0,0)]/(wk0_B0[OPS_ACC1(0,0)]*wk3_B0[OPS_ACC0(0,0)] -
      wk1_B0[OPS_ACC2(0,0)]*wk2_B0[OPS_ACC3(0,0)]);

    D10_B0[OPS_ACC7(0,0)] = -wk2_B0[OPS_ACC3(0,0)]/(wk0_B0[OPS_ACC1(0,0)]*wk3_B0[OPS_ACC0(0,0)] -
      wk1_B0[OPS_ACC2(0,0)]*wk2_B0[OPS_ACC3(0,0)]);

    D11_B0[OPS_ACC6(0,0)] = wk3_B0[OPS_ACC0(0,0)]/(wk0_B0[OPS_ACC1(0,0)]*wk3_B0[OPS_ACC0(0,0)] -
      wk1_B0[OPS_ACC2(0,0)]*wk2_B0[OPS_ACC3(0,0)]);

}

 void opensbliblock00Kernel030(const double *vorticity0_B0,const double *vorticity1_B0,double *vorticity1_old_B0,double
*vorticity0_old_B0)
{
   vorticity0_old_B0[OPS_ACC3(0,0)] = vorticity0_B0[OPS_ACC0(0,0)];

   vorticity1_old_B0[OPS_ACC2(0,0)] = vorticity1_B0[OPS_ACC1(0,0)];

}

void opensbliblock00Kernel005(double *Residual0_B0,double *Residual1_B0)
{
   Residual0_B0[OPS_ACC0(0,0)] = 0;

   Residual1_B0[OPS_ACC1(0,0)] = 0;

}

void opensbliblock00Kernel008(const double *velocity0_B0,double *wk2_B0, const int *idx)
{
wk2_B0[OPS_ACC1(0,0)] = ((idx[0] == 0) ? (
   (3.00000000000002*velocity0_B0[OPS_ACC0(1,0)] - 1.83333333333334*velocity0_B0[OPS_ACC0(0,0)] + 1.06910884386911e-15*velocity0_B0[OPS_ACC0(5,0)] + 0.333333333333356*velocity0_B0[OPS_ACC0(3,0)] - 8.34617916606957e-15*velocity0_B0[OPS_ACC0(4,0)] - 1.50000000000003*velocity0_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 1) ? (
   (0.719443173328855*velocity0_B0[OPS_ACC0(1,0)] - 0.322484932882161*velocity0_B0[OPS_ACC0(0,0)] - 0.376283677513354*velocity0_B0[OPS_ACC0(-1,0)] - 0.0658051057710389*velocity0_B0[OPS_ACC0(3,0)] + 0.00571369039775442*velocity0_B0[OPS_ACC0(4,0)] + 0.0394168524399448*velocity0_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 2) ? (
   (0.521455851089587*velocity0_B0[OPS_ACC0(1,0)] + 0.197184333887745*velocity0_B0[OPS_ACC0(0,0)] - 0.791245592765872*velocity0_B0[OPS_ACC0(-1,0)] + 0.113446470384241*velocity0_B0[OPS_ACC0(-2,0)] - 0.00412637789557492*velocity0_B0[OPS_ACC0(3,0)] - 0.0367146847001262*velocity0_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 3) ? (
   (0.652141084861241*velocity0_B0[OPS_ACC0(1,0)] + 0.0451033223343881*velocity0_B0[OPS_ACC0(0,0)] - 0.727822147724592*velocity0_B0[OPS_ACC0(-1,0)] + 0.121937153224065*velocity0_B0[OPS_ACC0(-2,0)] - 0.082033432844602*velocity0_B0[OPS_ACC0(2,0)] - 0.00932597985049999*velocity0_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 1) ? (
   -(-1.83333333333334*velocity0_B0[OPS_ACC0(0,0)] - 1.50000000000003*velocity0_B0[OPS_ACC0(-2,0)] + 3.00000000000002*velocity0_B0[OPS_ACC0(-1,0)] + 1.06910884386911e-15*velocity0_B0[OPS_ACC0(-5,0)] - 8.34617916606957e-15*velocity0_B0[OPS_ACC0(-4,0)] + 0.333333333333356*velocity0_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 2) ? (
   -(-0.376283677513354*velocity0_B0[OPS_ACC0(1,0)] - 0.322484932882161*velocity0_B0[OPS_ACC0(0,0)] + 0.0394168524399448*velocity0_B0[OPS_ACC0(-2,0)] + 0.719443173328855*velocity0_B0[OPS_ACC0(-1,0)] + 0.00571369039775442*velocity0_B0[OPS_ACC0(-4,0)] - 0.0658051057710389*velocity0_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 3) ? (
   -(-0.791245592765872*velocity0_B0[OPS_ACC0(1,0)] + 0.197184333887745*velocity0_B0[OPS_ACC0(0,0)] - 0.0367146847001262*velocity0_B0[OPS_ACC0(-2,0)] + 0.521455851089587*velocity0_B0[OPS_ACC0(-1,0)] + 0.113446470384241*velocity0_B0[OPS_ACC0(2,0)] - 0.00412637789557492*velocity0_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 4) ? (
   -(-0.727822147724592*velocity0_B0[OPS_ACC0(1,0)] + 0.0451033223343881*velocity0_B0[OPS_ACC0(0,0)] - 0.082033432844602*velocity0_B0[OPS_ACC0(-2,0)] + 0.652141084861241*velocity0_B0[OPS_ACC0(-1,0)] - 0.00932597985049999*velocity0_B0[OPS_ACC0(3,0)] + 0.121937153224065*velocity0_B0[OPS_ACC0(2,0)])/Delta0block0
)
: (
   (rc2)*velocity0_B0[OPS_ACC0(1,0)]/Delta0block0 - rc2*velocity0_B0[OPS_ACC0(-1,0)]/Delta0block0 + (rc3)*velocity0_B0[OPS_ACC0(-2,0)]/Delta0block0 - rc3*velocity0_B0[OPS_ACC0(2,0)]/Delta0block0
)))))))));

   // wk2_B0[OPS_ACC1(0,0)] = inv_0*((rc2)*velocity0_B0[OPS_ACC0(1,0)] - rc2*velocity0_B0[OPS_ACC0(-1,0)] +
   //   (rc3)*velocity0_B0[OPS_ACC0(-2,0)] - rc3*velocity0_B0[OPS_ACC0(2,0)]);

}

void opensbliblock00Kernel010(const double *velocity1_B0,double *wk3_B0, const int *idx)
{
wk3_B0[OPS_ACC1(0,0)] = ((idx[0] == 0) ? (
   (3.00000000000002*velocity1_B0[OPS_ACC0(1,0)] - 1.83333333333334*velocity1_B0[OPS_ACC0(0,0)] + 1.06910884386911e-15*velocity1_B0[OPS_ACC0(5,0)] + 0.333333333333356*velocity1_B0[OPS_ACC0(3,0)] - 8.34617916606957e-15*velocity1_B0[OPS_ACC0(4,0)] - 1.50000000000003*velocity1_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 1) ? (
   (0.719443173328855*velocity1_B0[OPS_ACC0(1,0)] - 0.322484932882161*velocity1_B0[OPS_ACC0(0,0)] - 0.376283677513354*velocity1_B0[OPS_ACC0(-1,0)] - 0.0658051057710389*velocity1_B0[OPS_ACC0(3,0)] + 0.00571369039775442*velocity1_B0[OPS_ACC0(4,0)] + 0.0394168524399448*velocity1_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 2) ? (
   (0.521455851089587*velocity1_B0[OPS_ACC0(1,0)] + 0.197184333887745*velocity1_B0[OPS_ACC0(0,0)] - 0.791245592765872*velocity1_B0[OPS_ACC0(-1,0)] + 0.113446470384241*velocity1_B0[OPS_ACC0(-2,0)] - 0.00412637789557492*velocity1_B0[OPS_ACC0(3,0)] - 0.0367146847001262*velocity1_B0[OPS_ACC0(2,0)])/Delta0block0
)
: ((idx[0] == 3) ? (
   (0.652141084861241*velocity1_B0[OPS_ACC0(1,0)] + 0.0451033223343881*velocity1_B0[OPS_ACC0(0,0)] - 0.727822147724592*velocity1_B0[OPS_ACC0(-1,0)] + 0.121937153224065*velocity1_B0[OPS_ACC0(-2,0)] - 0.082033432844602*velocity1_B0[OPS_ACC0(2,0)] - 0.00932597985049999*velocity1_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 1) ? (
   -(-1.83333333333334*velocity1_B0[OPS_ACC0(0,0)] - 1.50000000000003*velocity1_B0[OPS_ACC0(-2,0)] + 3.00000000000002*velocity1_B0[OPS_ACC0(-1,0)] + 1.06910884386911e-15*velocity1_B0[OPS_ACC0(-5,0)] - 8.34617916606957e-15*velocity1_B0[OPS_ACC0(-4,0)] + 0.333333333333356*velocity1_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 2) ? (
   -(-0.376283677513354*velocity1_B0[OPS_ACC0(1,0)] - 0.322484932882161*velocity1_B0[OPS_ACC0(0,0)] + 0.0394168524399448*velocity1_B0[OPS_ACC0(-2,0)] + 0.719443173328855*velocity1_B0[OPS_ACC0(-1,0)] + 0.00571369039775442*velocity1_B0[OPS_ACC0(-4,0)] - 0.0658051057710389*velocity1_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 3) ? (
   -(-0.791245592765872*velocity1_B0[OPS_ACC0(1,0)] + 0.197184333887745*velocity1_B0[OPS_ACC0(0,0)] - 0.0367146847001262*velocity1_B0[OPS_ACC0(-2,0)] + 0.521455851089587*velocity1_B0[OPS_ACC0(-1,0)] + 0.113446470384241*velocity1_B0[OPS_ACC0(2,0)] - 0.00412637789557492*velocity1_B0[OPS_ACC0(-3,0)])/Delta0block0
)
: ((idx[0] == block0np0 - 4) ? (
   -(-0.727822147724592*velocity1_B0[OPS_ACC0(1,0)] + 0.0451033223343881*velocity1_B0[OPS_ACC0(0,0)] - 0.082033432844602*velocity1_B0[OPS_ACC0(-2,0)] + 0.652141084861241*velocity1_B0[OPS_ACC0(-1,0)] - 0.00932597985049999*velocity1_B0[OPS_ACC0(3,0)] + 0.121937153224065*velocity1_B0[OPS_ACC0(2,0)])/Delta0block0
)
: (
   (rc2)*velocity1_B0[OPS_ACC0(1,0)]/Delta0block0 - rc2*velocity1_B0[OPS_ACC0(-1,0)]/Delta0block0 + (rc3)*velocity1_B0[OPS_ACC0(-2,0)]/Delta0block0 - rc3*velocity1_B0[OPS_ACC0(2,0)]/Delta0block0
)))))))));

    //wk3_B0[OPS_ACC1(0,0)] = inv_0*((rc2)*velocity1_B0[OPS_ACC0(1,0)] + (rc3)*velocity1_B0[OPS_ACC0(-2,0)] -
      //rc2*velocity1_B0[OPS_ACC0(-1,0)] - rc3*velocity1_B0[OPS_ACC0(2,0)]);

}

void opensbliblock00Kernel012(const double *velocity0_B0,double *wk4_B0, const int *idx)
{

wk4_B0[OPS_ACC1(0,0)] = ((idx[1] == 0) ? (
   (3.00000000000002*velocity0_B0[OPS_ACC0(0,1)] - 1.83333333333334*velocity0_B0[OPS_ACC0(0,0)] + 1.06910884386911e-15*velocity0_B0[OPS_ACC0(0,5)] + 0.333333333333356*velocity0_B0[OPS_ACC0(0,3)] - 8.34617916606957e-15*velocity0_B0[OPS_ACC0(0,4)] - 1.50000000000003*velocity0_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 1) ? (
   (0.719443173328855*velocity0_B0[OPS_ACC0(0,1)] - 0.322484932882161*velocity0_B0[OPS_ACC0(0,0)] - 0.376283677513354*velocity0_B0[OPS_ACC0(0,-1)] - 0.0658051057710389*velocity0_B0[OPS_ACC0(0,3)] + 0.00571369039775442*velocity0_B0[OPS_ACC0(0,4)] + 0.0394168524399448*velocity0_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 2) ? (
   (0.521455851089587*velocity0_B0[OPS_ACC0(0,1)] + 0.197184333887745*velocity0_B0[OPS_ACC0(0,0)] - 0.791245592765872*velocity0_B0[OPS_ACC0(0,-1)] + 0.113446470384241*velocity0_B0[OPS_ACC0(0,-2)] - 0.00412637789557492*velocity0_B0[OPS_ACC0(0,3)] - 0.0367146847001262*velocity0_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 3) ? (
   (0.652141084861241*velocity0_B0[OPS_ACC0(0,1)] + 0.0451033223343881*velocity0_B0[OPS_ACC0(0,0)] - 0.727822147724592*velocity0_B0[OPS_ACC0(0,-1)] + 0.121937153224065*velocity0_B0[OPS_ACC0(0,-2)] - 0.082033432844602*velocity0_B0[OPS_ACC0(0,2)] - 0.00932597985049999*velocity0_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 1) ? (
   -(-1.83333333333334*velocity0_B0[OPS_ACC0(0,0)] - 1.50000000000003*velocity0_B0[OPS_ACC0(0,-2)] + 3.00000000000002*velocity0_B0[OPS_ACC0(0,-1)] + 1.06910884386911e-15*velocity0_B0[OPS_ACC0(0,-5)] - 8.34617916606957e-15*velocity0_B0[OPS_ACC0(0,-4)] + 0.333333333333356*velocity0_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 2) ? (
   -(-0.376283677513354*velocity0_B0[OPS_ACC0(0,1)] - 0.322484932882161*velocity0_B0[OPS_ACC0(0,0)] + 0.0394168524399448*velocity0_B0[OPS_ACC0(0,-2)] + 0.719443173328855*velocity0_B0[OPS_ACC0(0,-1)] + 0.00571369039775442*velocity0_B0[OPS_ACC0(0,-4)] - 0.0658051057710389*velocity0_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 3) ? (
   -(-0.791245592765872*velocity0_B0[OPS_ACC0(0,1)] + 0.197184333887745*velocity0_B0[OPS_ACC0(0,0)] - 0.0367146847001262*velocity0_B0[OPS_ACC0(0,-2)] + 0.521455851089587*velocity0_B0[OPS_ACC0(0,-1)] + 0.113446470384241*velocity0_B0[OPS_ACC0(0,2)] - 0.00412637789557492*velocity0_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 4) ? (
   -(-0.727822147724592*velocity0_B0[OPS_ACC0(0,1)] + 0.0451033223343881*velocity0_B0[OPS_ACC0(0,0)] - 0.082033432844602*velocity0_B0[OPS_ACC0(0,-2)] + 0.652141084861241*velocity0_B0[OPS_ACC0(0,-1)] - 0.00932597985049999*velocity0_B0[OPS_ACC0(0,3)] + 0.121937153224065*velocity0_B0[OPS_ACC0(0,2)])/Delta1block0
)
: (
   (rc2)*velocity0_B0[OPS_ACC0(0,1)]/Delta1block0 - rc2*velocity0_B0[OPS_ACC0(0,-1)]/Delta1block0 + (rc3)*velocity0_B0[OPS_ACC0(0,-2)]/Delta1block0 - rc3*velocity0_B0[OPS_ACC0(0,2)]/Delta1block0
)))))))));

//    wk4_B0[OPS_ACC1(0,0)] = inv_1*((rc2)*velocity0_B0[OPS_ACC0(0,1)] - rc2*velocity0_B0[OPS_ACC0(0,-1)] +
  //    (rc3)*velocity0_B0[OPS_ACC0(0,-2)] - rc3*velocity0_B0[OPS_ACC0(0,2)]);

}

void opensbliblock00Kernel013(const double *velocity1_B0,double *wk5_B0, const int *idx)
{

wk5_B0[OPS_ACC1(0,0)] = ((idx[1] == 0) ? (
   (3.00000000000002*velocity1_B0[OPS_ACC0(0,1)] - 1.83333333333334*velocity1_B0[OPS_ACC0(0,0)] + 1.06910884386911e-15*velocity1_B0[OPS_ACC0(0,5)] + 0.333333333333356*velocity1_B0[OPS_ACC0(0,3)] - 8.34617916606957e-15*velocity1_B0[OPS_ACC0(0,4)] - 1.50000000000003*velocity1_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 1) ? (
   (0.719443173328855*velocity1_B0[OPS_ACC0(0,1)] - 0.322484932882161*velocity1_B0[OPS_ACC0(0,0)] - 0.376283677513354*velocity1_B0[OPS_ACC0(0,-1)] - 0.0658051057710389*velocity1_B0[OPS_ACC0(0,3)] + 0.00571369039775442*velocity1_B0[OPS_ACC0(0,4)] + 0.0394168524399448*velocity1_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 2) ? (
   (0.521455851089587*velocity1_B0[OPS_ACC0(0,1)] + 0.197184333887745*velocity1_B0[OPS_ACC0(0,0)] - 0.791245592765872*velocity1_B0[OPS_ACC0(0,-1)] + 0.113446470384241*velocity1_B0[OPS_ACC0(0,-2)] - 0.00412637789557492*velocity1_B0[OPS_ACC0(0,3)] - 0.0367146847001262*velocity1_B0[OPS_ACC0(0,2)])/Delta1block0
)
: ((idx[1] == 3) ? (
   (0.652141084861241*velocity1_B0[OPS_ACC0(0,1)] + 0.0451033223343881*velocity1_B0[OPS_ACC0(0,0)] - 0.727822147724592*velocity1_B0[OPS_ACC0(0,-1)] + 0.121937153224065*velocity1_B0[OPS_ACC0(0,-2)] - 0.082033432844602*velocity1_B0[OPS_ACC0(0,2)] - 0.00932597985049999*velocity1_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 1) ? (
   -(-1.83333333333334*velocity1_B0[OPS_ACC0(0,0)] - 1.50000000000003*velocity1_B0[OPS_ACC0(0,-2)] + 3.00000000000002*velocity1_B0[OPS_ACC0(0,-1)] + 1.06910884386911e-15*velocity1_B0[OPS_ACC0(0,-5)] - 8.34617916606957e-15*velocity1_B0[OPS_ACC0(0,-4)] + 0.333333333333356*velocity1_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 2) ? (
   -(-0.376283677513354*velocity1_B0[OPS_ACC0(0,1)] - 0.322484932882161*velocity1_B0[OPS_ACC0(0,0)] + 0.0394168524399448*velocity1_B0[OPS_ACC0(0,-2)] + 0.719443173328855*velocity1_B0[OPS_ACC0(0,-1)] + 0.00571369039775442*velocity1_B0[OPS_ACC0(0,-4)] - 0.0658051057710389*velocity1_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 3) ? (
   -(-0.791245592765872*velocity1_B0[OPS_ACC0(0,1)] + 0.197184333887745*velocity1_B0[OPS_ACC0(0,0)] - 0.0367146847001262*velocity1_B0[OPS_ACC0(0,-2)] + 0.521455851089587*velocity1_B0[OPS_ACC0(0,-1)] + 0.113446470384241*velocity1_B0[OPS_ACC0(0,2)] - 0.00412637789557492*velocity1_B0[OPS_ACC0(0,-3)])/Delta1block0
)
: ((idx[1] == block0np1 - 4) ? (
   -(-0.727822147724592*velocity1_B0[OPS_ACC0(0,1)] + 0.0451033223343881*velocity1_B0[OPS_ACC0(0,0)] - 0.082033432844602*velocity1_B0[OPS_ACC0(0,-2)] + 0.652141084861241*velocity1_B0[OPS_ACC0(0,-1)] - 0.00932597985049999*velocity1_B0[OPS_ACC0(0,3)] + 0.121937153224065*velocity1_B0[OPS_ACC0(0,2)])/Delta1block0
)
: (
   (rc2)*velocity1_B0[OPS_ACC0(0,1)]/Delta1block0 - rc2*velocity1_B0[OPS_ACC0(0,-1)]/Delta1block0 + (rc3)*velocity1_B0[OPS_ACC0(0,-2)]/Delta1block0 - rc3*velocity1_B0[OPS_ACC0(0,2)]/Delta1block0
)))))))));

    //wk5_B0[OPS_ACC1(0,0)] = inv_1*((rc2)*velocity1_B0[OPS_ACC0(0,1)] + (rc3)*velocity1_B0[OPS_ACC0(0,-2)] -
    //  rc2*velocity1_B0[OPS_ACC0(0,-1)] - rc3*velocity1_B0[OPS_ACC0(0,2)]);

}

 void opensbliblock00Kernel014(const double *D01_B0,const double *wk4_B0,const double *wk5_B0,const double *wk3_B0,const
double *D00_B0,const double *D11_B0,const double *D10_B0,const double *wk2_B0,double *Residual0_B0,double
*Residual1_B0)
{
    Residual0_B0[OPS_ACC8(0,0)] = D00_B0[OPS_ACC4(0,0)]*wk3_B0[OPS_ACC3(0,0)] -
      D01_B0[OPS_ACC0(0,0)]*wk2_B0[OPS_ACC7(0,0)] + D10_B0[OPS_ACC6(0,0)]*wk5_B0[OPS_ACC2(0,0)] -
      D11_B0[OPS_ACC5(0,0)]*wk4_B0[OPS_ACC1(0,0)] + Residual0_B0[OPS_ACC8(0,0)];

    Residual1_B0[OPS_ACC9(0,0)] = D00_B0[OPS_ACC4(0,0)]*wk3_B0[OPS_ACC3(0,0)] -
      D01_B0[OPS_ACC0(0,0)]*wk2_B0[OPS_ACC7(0,0)] + D10_B0[OPS_ACC6(0,0)]*wk5_B0[OPS_ACC2(0,0)] -
      D11_B0[OPS_ACC5(0,0)]*wk4_B0[OPS_ACC1(0,0)] + Residual1_B0[OPS_ACC9(0,0)];

}

 void opensbliblock00Kernel032(const double *vorticity1_old_B0,const double *vorticity0_old_B0,const double
*Residual0_B0,const double *Residual1_B0,double *vorticity0_B0,double *vorticity1_B0,const double *rknew)
{
   vorticity0_B0[OPS_ACC4(0,0)] = Residual0_B0[OPS_ACC2(0,0)];

   vorticity1_B0[OPS_ACC5(0,0)] = Residual1_B0[OPS_ACC3(0,0)];

}

 void opensbliblock00Kernel031(const double *Residual0_B0,const double *Residual1_B0,double *vorticity1_old_B0,double
*vorticity0_old_B0,const double *rkold)
{
   vorticity0_old_B0[OPS_ACC3(0,0)] = Residual0_B0[OPS_ACC0(0,0)] + vorticity0_old_B0[OPS_ACC3(0,0)];

   vorticity1_old_B0[OPS_ACC2(0,0)] = Residual1_B0[OPS_ACC1(0,0)] + vorticity1_old_B0[OPS_ACC2(0,0)];

}

void opensbliblock00KernelErrorIndicator(const double *vorticity, const double *x0, const double *x1, const double *hamming_window, const int *ldir, const int *start_idx, int *idx, double *N2sum, double *N4sum_r, double *N4sum_i, double *N8sum_r, double *N8sum_i, double *sbx0, double *sbx1)
{
    // Error indicator kernel.
    // For more details see: C. T. Jacobs, M. Zauner, N. De Tullio, Satya P. Jammy, David J. Lusher, N. D. Sandham (Submitted). "An error indicator for finite difference methods using spectral techniques with application to aerofoil simulation". Submitted to the ParCFD 2017 Special Issue of Computers & Fluids.

    double filtered;
    int index;

    const int nx_sb2 = 16;
    double hamming_window2[nx_sb2];
    
    double __complex__ z1 = 1.0I;
    double __complex__ z = 0.0 - 1.0I;
    double __complex__ zz = exp (__real__ (-z1*M_PI/4.0)) * (cos (__imag__ (-z1*M_PI/4.0)) + z1*sin(__imag__ (-z1*M_PI/4.0)));

    double __complex__ coefficient4 = z;
    double __complex__ coefficient8 = zz;

    //ops_printf("index = %d, %d, %f\n", idx[0], idx[1], vorticity[OPS_ACC0(0,0)]);
    
    for (int j=0; j<nx_sb2; j++){
       hamming_window2[j] = (0.54 - 0.46*cos(2*M_PI*j/(nx_sb2)))/0.54;
    }    

    if (*ldir == 0)
    {
        index = (idx[0] - *start_idx) % nx_sb2;
    }
    else
    {
        if (*ldir == 1)
        {
            index = (idx[1] - *start_idx) % nx_sb2;
        }
    }

    filtered = vorticity[OPS_ACC0(0,0)]*hamming_window2[index];
    //ops_printf("filtered[%d] = %lf\n", index, filtered);
    
    // Reconstruct the Fourier amplitudes of selected modes (N/2, N/4 and N/8) by simple summations.
    *N2sum = *N2sum + pow(-1.0, index)*filtered;
    //ops_printf("N2sum = %lf, %lf, %lf\n", pow(-1, index), filtered, *N2sum);
    if (index == 0)
    {
        *N4sum_r = *N4sum_r + filtered;
        *N8sum_r = *N8sum_r + filtered;
    }
    else
    {
        // Since we can't raise complex numbers to powers via a function, let's just do it ourselves here by multiplying by itself each time in the loop.
        for (int k=1; k<index; k++)
        {
            coefficient4 = coefficient4*z;
            coefficient8 = coefficient8*zz;
        }
        *N4sum_r = *N4sum_r + __real__ coefficient4*filtered;
        *N4sum_i = *N4sum_i + __imag__ coefficient4*filtered;
        *N8sum_r = *N8sum_r + __real__ coefficient8*filtered;
        *N8sum_i = *N8sum_i + __imag__ coefficient8*filtered;
    }

    *sbx0 += x0[OPS_ACC1(0,0)]/nx_sb2;
    *sbx1 += x1[OPS_ACC2(0,0)]/nx_sb2;
}

#endif
