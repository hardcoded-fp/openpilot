#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6288688874102308164) {
   out_6288688874102308164[0] = delta_x[0] + nom_x[0];
   out_6288688874102308164[1] = delta_x[1] + nom_x[1];
   out_6288688874102308164[2] = delta_x[2] + nom_x[2];
   out_6288688874102308164[3] = delta_x[3] + nom_x[3];
   out_6288688874102308164[4] = delta_x[4] + nom_x[4];
   out_6288688874102308164[5] = delta_x[5] + nom_x[5];
   out_6288688874102308164[6] = delta_x[6] + nom_x[6];
   out_6288688874102308164[7] = delta_x[7] + nom_x[7];
   out_6288688874102308164[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5232615150247248436) {
   out_5232615150247248436[0] = -nom_x[0] + true_x[0];
   out_5232615150247248436[1] = -nom_x[1] + true_x[1];
   out_5232615150247248436[2] = -nom_x[2] + true_x[2];
   out_5232615150247248436[3] = -nom_x[3] + true_x[3];
   out_5232615150247248436[4] = -nom_x[4] + true_x[4];
   out_5232615150247248436[5] = -nom_x[5] + true_x[5];
   out_5232615150247248436[6] = -nom_x[6] + true_x[6];
   out_5232615150247248436[7] = -nom_x[7] + true_x[7];
   out_5232615150247248436[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4159914697228368279) {
   out_4159914697228368279[0] = 1.0;
   out_4159914697228368279[1] = 0;
   out_4159914697228368279[2] = 0;
   out_4159914697228368279[3] = 0;
   out_4159914697228368279[4] = 0;
   out_4159914697228368279[5] = 0;
   out_4159914697228368279[6] = 0;
   out_4159914697228368279[7] = 0;
   out_4159914697228368279[8] = 0;
   out_4159914697228368279[9] = 0;
   out_4159914697228368279[10] = 1.0;
   out_4159914697228368279[11] = 0;
   out_4159914697228368279[12] = 0;
   out_4159914697228368279[13] = 0;
   out_4159914697228368279[14] = 0;
   out_4159914697228368279[15] = 0;
   out_4159914697228368279[16] = 0;
   out_4159914697228368279[17] = 0;
   out_4159914697228368279[18] = 0;
   out_4159914697228368279[19] = 0;
   out_4159914697228368279[20] = 1.0;
   out_4159914697228368279[21] = 0;
   out_4159914697228368279[22] = 0;
   out_4159914697228368279[23] = 0;
   out_4159914697228368279[24] = 0;
   out_4159914697228368279[25] = 0;
   out_4159914697228368279[26] = 0;
   out_4159914697228368279[27] = 0;
   out_4159914697228368279[28] = 0;
   out_4159914697228368279[29] = 0;
   out_4159914697228368279[30] = 1.0;
   out_4159914697228368279[31] = 0;
   out_4159914697228368279[32] = 0;
   out_4159914697228368279[33] = 0;
   out_4159914697228368279[34] = 0;
   out_4159914697228368279[35] = 0;
   out_4159914697228368279[36] = 0;
   out_4159914697228368279[37] = 0;
   out_4159914697228368279[38] = 0;
   out_4159914697228368279[39] = 0;
   out_4159914697228368279[40] = 1.0;
   out_4159914697228368279[41] = 0;
   out_4159914697228368279[42] = 0;
   out_4159914697228368279[43] = 0;
   out_4159914697228368279[44] = 0;
   out_4159914697228368279[45] = 0;
   out_4159914697228368279[46] = 0;
   out_4159914697228368279[47] = 0;
   out_4159914697228368279[48] = 0;
   out_4159914697228368279[49] = 0;
   out_4159914697228368279[50] = 1.0;
   out_4159914697228368279[51] = 0;
   out_4159914697228368279[52] = 0;
   out_4159914697228368279[53] = 0;
   out_4159914697228368279[54] = 0;
   out_4159914697228368279[55] = 0;
   out_4159914697228368279[56] = 0;
   out_4159914697228368279[57] = 0;
   out_4159914697228368279[58] = 0;
   out_4159914697228368279[59] = 0;
   out_4159914697228368279[60] = 1.0;
   out_4159914697228368279[61] = 0;
   out_4159914697228368279[62] = 0;
   out_4159914697228368279[63] = 0;
   out_4159914697228368279[64] = 0;
   out_4159914697228368279[65] = 0;
   out_4159914697228368279[66] = 0;
   out_4159914697228368279[67] = 0;
   out_4159914697228368279[68] = 0;
   out_4159914697228368279[69] = 0;
   out_4159914697228368279[70] = 1.0;
   out_4159914697228368279[71] = 0;
   out_4159914697228368279[72] = 0;
   out_4159914697228368279[73] = 0;
   out_4159914697228368279[74] = 0;
   out_4159914697228368279[75] = 0;
   out_4159914697228368279[76] = 0;
   out_4159914697228368279[77] = 0;
   out_4159914697228368279[78] = 0;
   out_4159914697228368279[79] = 0;
   out_4159914697228368279[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_4156207149520197820) {
   out_4156207149520197820[0] = state[0];
   out_4156207149520197820[1] = state[1];
   out_4156207149520197820[2] = state[2];
   out_4156207149520197820[3] = state[3];
   out_4156207149520197820[4] = state[4];
   out_4156207149520197820[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4156207149520197820[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4156207149520197820[7] = state[7];
   out_4156207149520197820[8] = state[8];
}
void F_fun(double *state, double dt, double *out_1014683453205324692) {
   out_1014683453205324692[0] = 1;
   out_1014683453205324692[1] = 0;
   out_1014683453205324692[2] = 0;
   out_1014683453205324692[3] = 0;
   out_1014683453205324692[4] = 0;
   out_1014683453205324692[5] = 0;
   out_1014683453205324692[6] = 0;
   out_1014683453205324692[7] = 0;
   out_1014683453205324692[8] = 0;
   out_1014683453205324692[9] = 0;
   out_1014683453205324692[10] = 1;
   out_1014683453205324692[11] = 0;
   out_1014683453205324692[12] = 0;
   out_1014683453205324692[13] = 0;
   out_1014683453205324692[14] = 0;
   out_1014683453205324692[15] = 0;
   out_1014683453205324692[16] = 0;
   out_1014683453205324692[17] = 0;
   out_1014683453205324692[18] = 0;
   out_1014683453205324692[19] = 0;
   out_1014683453205324692[20] = 1;
   out_1014683453205324692[21] = 0;
   out_1014683453205324692[22] = 0;
   out_1014683453205324692[23] = 0;
   out_1014683453205324692[24] = 0;
   out_1014683453205324692[25] = 0;
   out_1014683453205324692[26] = 0;
   out_1014683453205324692[27] = 0;
   out_1014683453205324692[28] = 0;
   out_1014683453205324692[29] = 0;
   out_1014683453205324692[30] = 1;
   out_1014683453205324692[31] = 0;
   out_1014683453205324692[32] = 0;
   out_1014683453205324692[33] = 0;
   out_1014683453205324692[34] = 0;
   out_1014683453205324692[35] = 0;
   out_1014683453205324692[36] = 0;
   out_1014683453205324692[37] = 0;
   out_1014683453205324692[38] = 0;
   out_1014683453205324692[39] = 0;
   out_1014683453205324692[40] = 1;
   out_1014683453205324692[41] = 0;
   out_1014683453205324692[42] = 0;
   out_1014683453205324692[43] = 0;
   out_1014683453205324692[44] = 0;
   out_1014683453205324692[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1014683453205324692[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1014683453205324692[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1014683453205324692[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1014683453205324692[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1014683453205324692[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1014683453205324692[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1014683453205324692[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1014683453205324692[53] = -9.8000000000000007*dt;
   out_1014683453205324692[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1014683453205324692[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1014683453205324692[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1014683453205324692[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1014683453205324692[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1014683453205324692[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1014683453205324692[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1014683453205324692[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1014683453205324692[62] = 0;
   out_1014683453205324692[63] = 0;
   out_1014683453205324692[64] = 0;
   out_1014683453205324692[65] = 0;
   out_1014683453205324692[66] = 0;
   out_1014683453205324692[67] = 0;
   out_1014683453205324692[68] = 0;
   out_1014683453205324692[69] = 0;
   out_1014683453205324692[70] = 1;
   out_1014683453205324692[71] = 0;
   out_1014683453205324692[72] = 0;
   out_1014683453205324692[73] = 0;
   out_1014683453205324692[74] = 0;
   out_1014683453205324692[75] = 0;
   out_1014683453205324692[76] = 0;
   out_1014683453205324692[77] = 0;
   out_1014683453205324692[78] = 0;
   out_1014683453205324692[79] = 0;
   out_1014683453205324692[80] = 1;
}
void h_25(double *state, double *unused, double *out_5671253352182223489) {
   out_5671253352182223489[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2927558261712726917) {
   out_2927558261712726917[0] = 0;
   out_2927558261712726917[1] = 0;
   out_2927558261712726917[2] = 0;
   out_2927558261712726917[3] = 0;
   out_2927558261712726917[4] = 0;
   out_2927558261712726917[5] = 0;
   out_2927558261712726917[6] = 1;
   out_2927558261712726917[7] = 0;
   out_2927558261712726917[8] = 0;
}
void h_24(double *state, double *unused, double *out_3382553279571021005) {
   out_3382553279571021005[0] = state[4];
   out_3382553279571021005[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7674067588679235951) {
   out_7674067588679235951[0] = 0;
   out_7674067588679235951[1] = 0;
   out_7674067588679235951[2] = 0;
   out_7674067588679235951[3] = 0;
   out_7674067588679235951[4] = 1;
   out_7674067588679235951[5] = 0;
   out_7674067588679235951[6] = 0;
   out_7674067588679235951[7] = 0;
   out_7674067588679235951[8] = 0;
   out_7674067588679235951[9] = 0;
   out_7674067588679235951[10] = 0;
   out_7674067588679235951[11] = 0;
   out_7674067588679235951[12] = 0;
   out_7674067588679235951[13] = 0;
   out_7674067588679235951[14] = 1;
   out_7674067588679235951[15] = 0;
   out_7674067588679235951[16] = 0;
   out_7674067588679235951[17] = 0;
}
void h_30(double *state, double *unused, double *out_6004655495176977191) {
   out_6004655495176977191[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3989132079778889838) {
   out_3989132079778889838[0] = 0;
   out_3989132079778889838[1] = 0;
   out_3989132079778889838[2] = 0;
   out_3989132079778889838[3] = 0;
   out_3989132079778889838[4] = 1;
   out_3989132079778889838[5] = 0;
   out_3989132079778889838[6] = 0;
   out_3989132079778889838[7] = 0;
   out_3989132079778889838[8] = 0;
}
void h_26(double *state, double *unused, double *out_1002982227516131546) {
   out_1002982227516131546[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6669061580586783141) {
   out_6669061580586783141[0] = 0;
   out_6669061580586783141[1] = 0;
   out_6669061580586783141[2] = 0;
   out_6669061580586783141[3] = 0;
   out_6669061580586783141[4] = 0;
   out_6669061580586783141[5] = 0;
   out_6669061580586783141[6] = 0;
   out_6669061580586783141[7] = 1;
   out_6669061580586783141[8] = 0;
}
void h_27(double *state, double *unused, double *out_7361701186917904449) {
   out_7361701186917904449[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1814368767978464927) {
   out_1814368767978464927[0] = 0;
   out_1814368767978464927[1] = 0;
   out_1814368767978464927[2] = 0;
   out_1814368767978464927[3] = 1;
   out_1814368767978464927[4] = 0;
   out_1814368767978464927[5] = 0;
   out_1814368767978464927[6] = 0;
   out_1814368767978464927[7] = 0;
   out_1814368767978464927[8] = 0;
}
void h_29(double *state, double *unused, double *out_2401001324467249462) {
   out_2401001324467249462[0] = state[1];
}
void H_29(double *state, double *unused, double *out_101006041108913894) {
   out_101006041108913894[0] = 0;
   out_101006041108913894[1] = 1;
   out_101006041108913894[2] = 0;
   out_101006041108913894[3] = 0;
   out_101006041108913894[4] = 0;
   out_101006041108913894[5] = 0;
   out_101006041108913894[6] = 0;
   out_101006041108913894[7] = 0;
   out_101006041108913894[8] = 0;
}
void h_28(double *state, double *unused, double *out_2181546259182111216) {
   out_2181546259182111216[0] = state[0];
}
void H_28(double *state, double *unused, double *out_4981392975960616680) {
   out_4981392975960616680[0] = 1;
   out_4981392975960616680[1] = 0;
   out_4981392975960616680[2] = 0;
   out_4981392975960616680[3] = 0;
   out_4981392975960616680[4] = 0;
   out_4981392975960616680[5] = 0;
   out_4981392975960616680[6] = 0;
   out_4981392975960616680[7] = 0;
   out_4981392975960616680[8] = 0;
}
void h_31(double *state, double *unused, double *out_425985762948188827) {
   out_425985762948188827[0] = state[8];
}
void H_31(double *state, double *unused, double *out_2896912299835766489) {
   out_2896912299835766489[0] = 0;
   out_2896912299835766489[1] = 0;
   out_2896912299835766489[2] = 0;
   out_2896912299835766489[3] = 0;
   out_2896912299835766489[4] = 0;
   out_2896912299835766489[5] = 0;
   out_2896912299835766489[6] = 0;
   out_2896912299835766489[7] = 0;
   out_2896912299835766489[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_6288688874102308164) {
  err_fun(nom_x, delta_x, out_6288688874102308164);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5232615150247248436) {
  inv_err_fun(nom_x, true_x, out_5232615150247248436);
}
void car_H_mod_fun(double *state, double *out_4159914697228368279) {
  H_mod_fun(state, out_4159914697228368279);
}
void car_f_fun(double *state, double dt, double *out_4156207149520197820) {
  f_fun(state,  dt, out_4156207149520197820);
}
void car_F_fun(double *state, double dt, double *out_1014683453205324692) {
  F_fun(state,  dt, out_1014683453205324692);
}
void car_h_25(double *state, double *unused, double *out_5671253352182223489) {
  h_25(state, unused, out_5671253352182223489);
}
void car_H_25(double *state, double *unused, double *out_2927558261712726917) {
  H_25(state, unused, out_2927558261712726917);
}
void car_h_24(double *state, double *unused, double *out_3382553279571021005) {
  h_24(state, unused, out_3382553279571021005);
}
void car_H_24(double *state, double *unused, double *out_7674067588679235951) {
  H_24(state, unused, out_7674067588679235951);
}
void car_h_30(double *state, double *unused, double *out_6004655495176977191) {
  h_30(state, unused, out_6004655495176977191);
}
void car_H_30(double *state, double *unused, double *out_3989132079778889838) {
  H_30(state, unused, out_3989132079778889838);
}
void car_h_26(double *state, double *unused, double *out_1002982227516131546) {
  h_26(state, unused, out_1002982227516131546);
}
void car_H_26(double *state, double *unused, double *out_6669061580586783141) {
  H_26(state, unused, out_6669061580586783141);
}
void car_h_27(double *state, double *unused, double *out_7361701186917904449) {
  h_27(state, unused, out_7361701186917904449);
}
void car_H_27(double *state, double *unused, double *out_1814368767978464927) {
  H_27(state, unused, out_1814368767978464927);
}
void car_h_29(double *state, double *unused, double *out_2401001324467249462) {
  h_29(state, unused, out_2401001324467249462);
}
void car_H_29(double *state, double *unused, double *out_101006041108913894) {
  H_29(state, unused, out_101006041108913894);
}
void car_h_28(double *state, double *unused, double *out_2181546259182111216) {
  h_28(state, unused, out_2181546259182111216);
}
void car_H_28(double *state, double *unused, double *out_4981392975960616680) {
  H_28(state, unused, out_4981392975960616680);
}
void car_h_31(double *state, double *unused, double *out_425985762948188827) {
  h_31(state, unused, out_425985762948188827);
}
void car_H_31(double *state, double *unused, double *out_2896912299835766489) {
  H_31(state, unused, out_2896912299835766489);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
