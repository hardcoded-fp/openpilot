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
void err_fun(double *nom_x, double *delta_x, double *out_3962033384618978572) {
   out_3962033384618978572[0] = delta_x[0] + nom_x[0];
   out_3962033384618978572[1] = delta_x[1] + nom_x[1];
   out_3962033384618978572[2] = delta_x[2] + nom_x[2];
   out_3962033384618978572[3] = delta_x[3] + nom_x[3];
   out_3962033384618978572[4] = delta_x[4] + nom_x[4];
   out_3962033384618978572[5] = delta_x[5] + nom_x[5];
   out_3962033384618978572[6] = delta_x[6] + nom_x[6];
   out_3962033384618978572[7] = delta_x[7] + nom_x[7];
   out_3962033384618978572[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8377453335461740860) {
   out_8377453335461740860[0] = -nom_x[0] + true_x[0];
   out_8377453335461740860[1] = -nom_x[1] + true_x[1];
   out_8377453335461740860[2] = -nom_x[2] + true_x[2];
   out_8377453335461740860[3] = -nom_x[3] + true_x[3];
   out_8377453335461740860[4] = -nom_x[4] + true_x[4];
   out_8377453335461740860[5] = -nom_x[5] + true_x[5];
   out_8377453335461740860[6] = -nom_x[6] + true_x[6];
   out_8377453335461740860[7] = -nom_x[7] + true_x[7];
   out_8377453335461740860[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_663431017939020115) {
   out_663431017939020115[0] = 1.0;
   out_663431017939020115[1] = 0;
   out_663431017939020115[2] = 0;
   out_663431017939020115[3] = 0;
   out_663431017939020115[4] = 0;
   out_663431017939020115[5] = 0;
   out_663431017939020115[6] = 0;
   out_663431017939020115[7] = 0;
   out_663431017939020115[8] = 0;
   out_663431017939020115[9] = 0;
   out_663431017939020115[10] = 1.0;
   out_663431017939020115[11] = 0;
   out_663431017939020115[12] = 0;
   out_663431017939020115[13] = 0;
   out_663431017939020115[14] = 0;
   out_663431017939020115[15] = 0;
   out_663431017939020115[16] = 0;
   out_663431017939020115[17] = 0;
   out_663431017939020115[18] = 0;
   out_663431017939020115[19] = 0;
   out_663431017939020115[20] = 1.0;
   out_663431017939020115[21] = 0;
   out_663431017939020115[22] = 0;
   out_663431017939020115[23] = 0;
   out_663431017939020115[24] = 0;
   out_663431017939020115[25] = 0;
   out_663431017939020115[26] = 0;
   out_663431017939020115[27] = 0;
   out_663431017939020115[28] = 0;
   out_663431017939020115[29] = 0;
   out_663431017939020115[30] = 1.0;
   out_663431017939020115[31] = 0;
   out_663431017939020115[32] = 0;
   out_663431017939020115[33] = 0;
   out_663431017939020115[34] = 0;
   out_663431017939020115[35] = 0;
   out_663431017939020115[36] = 0;
   out_663431017939020115[37] = 0;
   out_663431017939020115[38] = 0;
   out_663431017939020115[39] = 0;
   out_663431017939020115[40] = 1.0;
   out_663431017939020115[41] = 0;
   out_663431017939020115[42] = 0;
   out_663431017939020115[43] = 0;
   out_663431017939020115[44] = 0;
   out_663431017939020115[45] = 0;
   out_663431017939020115[46] = 0;
   out_663431017939020115[47] = 0;
   out_663431017939020115[48] = 0;
   out_663431017939020115[49] = 0;
   out_663431017939020115[50] = 1.0;
   out_663431017939020115[51] = 0;
   out_663431017939020115[52] = 0;
   out_663431017939020115[53] = 0;
   out_663431017939020115[54] = 0;
   out_663431017939020115[55] = 0;
   out_663431017939020115[56] = 0;
   out_663431017939020115[57] = 0;
   out_663431017939020115[58] = 0;
   out_663431017939020115[59] = 0;
   out_663431017939020115[60] = 1.0;
   out_663431017939020115[61] = 0;
   out_663431017939020115[62] = 0;
   out_663431017939020115[63] = 0;
   out_663431017939020115[64] = 0;
   out_663431017939020115[65] = 0;
   out_663431017939020115[66] = 0;
   out_663431017939020115[67] = 0;
   out_663431017939020115[68] = 0;
   out_663431017939020115[69] = 0;
   out_663431017939020115[70] = 1.0;
   out_663431017939020115[71] = 0;
   out_663431017939020115[72] = 0;
   out_663431017939020115[73] = 0;
   out_663431017939020115[74] = 0;
   out_663431017939020115[75] = 0;
   out_663431017939020115[76] = 0;
   out_663431017939020115[77] = 0;
   out_663431017939020115[78] = 0;
   out_663431017939020115[79] = 0;
   out_663431017939020115[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2189594786826403342) {
   out_2189594786826403342[0] = state[0];
   out_2189594786826403342[1] = state[1];
   out_2189594786826403342[2] = state[2];
   out_2189594786826403342[3] = state[3];
   out_2189594786826403342[4] = state[4];
   out_2189594786826403342[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2189594786826403342[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2189594786826403342[7] = state[7];
   out_2189594786826403342[8] = state[8];
}
void F_fun(double *state, double dt, double *out_2972862917221800866) {
   out_2972862917221800866[0] = 1;
   out_2972862917221800866[1] = 0;
   out_2972862917221800866[2] = 0;
   out_2972862917221800866[3] = 0;
   out_2972862917221800866[4] = 0;
   out_2972862917221800866[5] = 0;
   out_2972862917221800866[6] = 0;
   out_2972862917221800866[7] = 0;
   out_2972862917221800866[8] = 0;
   out_2972862917221800866[9] = 0;
   out_2972862917221800866[10] = 1;
   out_2972862917221800866[11] = 0;
   out_2972862917221800866[12] = 0;
   out_2972862917221800866[13] = 0;
   out_2972862917221800866[14] = 0;
   out_2972862917221800866[15] = 0;
   out_2972862917221800866[16] = 0;
   out_2972862917221800866[17] = 0;
   out_2972862917221800866[18] = 0;
   out_2972862917221800866[19] = 0;
   out_2972862917221800866[20] = 1;
   out_2972862917221800866[21] = 0;
   out_2972862917221800866[22] = 0;
   out_2972862917221800866[23] = 0;
   out_2972862917221800866[24] = 0;
   out_2972862917221800866[25] = 0;
   out_2972862917221800866[26] = 0;
   out_2972862917221800866[27] = 0;
   out_2972862917221800866[28] = 0;
   out_2972862917221800866[29] = 0;
   out_2972862917221800866[30] = 1;
   out_2972862917221800866[31] = 0;
   out_2972862917221800866[32] = 0;
   out_2972862917221800866[33] = 0;
   out_2972862917221800866[34] = 0;
   out_2972862917221800866[35] = 0;
   out_2972862917221800866[36] = 0;
   out_2972862917221800866[37] = 0;
   out_2972862917221800866[38] = 0;
   out_2972862917221800866[39] = 0;
   out_2972862917221800866[40] = 1;
   out_2972862917221800866[41] = 0;
   out_2972862917221800866[42] = 0;
   out_2972862917221800866[43] = 0;
   out_2972862917221800866[44] = 0;
   out_2972862917221800866[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2972862917221800866[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2972862917221800866[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2972862917221800866[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2972862917221800866[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2972862917221800866[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2972862917221800866[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2972862917221800866[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2972862917221800866[53] = -9.8000000000000007*dt;
   out_2972862917221800866[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2972862917221800866[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2972862917221800866[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2972862917221800866[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2972862917221800866[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2972862917221800866[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2972862917221800866[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2972862917221800866[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2972862917221800866[62] = 0;
   out_2972862917221800866[63] = 0;
   out_2972862917221800866[64] = 0;
   out_2972862917221800866[65] = 0;
   out_2972862917221800866[66] = 0;
   out_2972862917221800866[67] = 0;
   out_2972862917221800866[68] = 0;
   out_2972862917221800866[69] = 0;
   out_2972862917221800866[70] = 1;
   out_2972862917221800866[71] = 0;
   out_2972862917221800866[72] = 0;
   out_2972862917221800866[73] = 0;
   out_2972862917221800866[74] = 0;
   out_2972862917221800866[75] = 0;
   out_2972862917221800866[76] = 0;
   out_2972862917221800866[77] = 0;
   out_2972862917221800866[78] = 0;
   out_2972862917221800866[79] = 0;
   out_2972862917221800866[80] = 1;
}
void h_25(double *state, double *unused, double *out_5656616099551326110) {
   out_5656616099551326110[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6913354757324836494) {
   out_6913354757324836494[0] = 0;
   out_6913354757324836494[1] = 0;
   out_6913354757324836494[2] = 0;
   out_6913354757324836494[3] = 0;
   out_6913354757324836494[4] = 0;
   out_6913354757324836494[5] = 0;
   out_6913354757324836494[6] = 1;
   out_6913354757324836494[7] = 0;
   out_6913354757324836494[8] = 0;
}
void h_24(double *state, double *unused, double *out_4997502995021187120) {
   out_4997502995021187120[0] = state[4];
   out_4997502995021187120[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4297999924517583471) {
   out_4297999924517583471[0] = 0;
   out_4297999924517583471[1] = 0;
   out_4297999924517583471[2] = 0;
   out_4297999924517583471[3] = 0;
   out_4297999924517583471[4] = 1;
   out_4297999924517583471[5] = 0;
   out_4297999924517583471[6] = 0;
   out_4297999924517583471[7] = 0;
   out_4297999924517583471[8] = 0;
   out_4297999924517583471[9] = 0;
   out_4297999924517583471[10] = 0;
   out_4297999924517583471[11] = 0;
   out_4297999924517583471[12] = 0;
   out_4297999924517583471[13] = 0;
   out_4297999924517583471[14] = 1;
   out_4297999924517583471[15] = 0;
   out_4297999924517583471[16] = 0;
   out_4297999924517583471[17] = 0;
}
void h_30(double *state, double *unused, double *out_4141240455864945865) {
   out_4141240455864945865[0] = state[4];
}
void H_30(double *state, double *unused, double *out_2385658427197228296) {
   out_2385658427197228296[0] = 0;
   out_2385658427197228296[1] = 0;
   out_2385658427197228296[2] = 0;
   out_2385658427197228296[3] = 0;
   out_2385658427197228296[4] = 1;
   out_2385658427197228296[5] = 0;
   out_2385658427197228296[6] = 0;
   out_2385658427197228296[7] = 0;
   out_2385658427197228296[8] = 0;
}
void h_26(double *state, double *unused, double *out_860432811795899780) {
   out_860432811795899780[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3171851438450780270) {
   out_3171851438450780270[0] = 0;
   out_3171851438450780270[1] = 0;
   out_3171851438450780270[2] = 0;
   out_3171851438450780270[3] = 0;
   out_3171851438450780270[4] = 0;
   out_3171851438450780270[5] = 0;
   out_3171851438450780270[6] = 0;
   out_3171851438450780270[7] = 1;
   out_3171851438450780270[8] = 0;
}
void h_27(double *state, double *unused, double *out_7788454706497379861) {
   out_7788454706497379861[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4609252498381171513) {
   out_4609252498381171513[0] = 0;
   out_4609252498381171513[1] = 0;
   out_4609252498381171513[2] = 0;
   out_4609252498381171513[3] = 1;
   out_4609252498381171513[4] = 0;
   out_4609252498381171513[5] = 0;
   out_4609252498381171513[6] = 0;
   out_4609252498381171513[7] = 0;
   out_4609252498381171513[8] = 0;
}
void h_29(double *state, double *unused, double *out_4840728716141894900) {
   out_4840728716141894900[0] = state[1];
}
void H_29(double *state, double *unused, double *out_2895889771511620480) {
   out_2895889771511620480[0] = 0;
   out_2895889771511620480[1] = 1;
   out_2895889771511620480[2] = 0;
   out_2895889771511620480[3] = 0;
   out_2895889771511620480[4] = 0;
   out_2895889771511620480[5] = 0;
   out_2895889771511620480[6] = 0;
   out_2895889771511620480[7] = 0;
   out_2895889771511620480[8] = 0;
}
void h_28(double *state, double *unused, double *out_2786781438836534216) {
   out_2786781438836534216[0] = state[0];
}
void H_28(double *state, double *unused, double *out_4859520043076946731) {
   out_4859520043076946731[0] = 1;
   out_4859520043076946731[1] = 0;
   out_4859520043076946731[2] = 0;
   out_4859520043076946731[3] = 0;
   out_4859520043076946731[4] = 0;
   out_4859520043076946731[5] = 0;
   out_4859520043076946731[6] = 0;
   out_4859520043076946731[7] = 0;
   out_4859520043076946731[8] = 0;
}
void h_31(double *state, double *unused, double *out_1794932623122824458) {
   out_1794932623122824458[0] = state[8];
}
void H_31(double *state, double *unused, double *out_6944000719201796922) {
   out_6944000719201796922[0] = 0;
   out_6944000719201796922[1] = 0;
   out_6944000719201796922[2] = 0;
   out_6944000719201796922[3] = 0;
   out_6944000719201796922[4] = 0;
   out_6944000719201796922[5] = 0;
   out_6944000719201796922[6] = 0;
   out_6944000719201796922[7] = 0;
   out_6944000719201796922[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_3962033384618978572) {
  err_fun(nom_x, delta_x, out_3962033384618978572);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8377453335461740860) {
  inv_err_fun(nom_x, true_x, out_8377453335461740860);
}
void car_H_mod_fun(double *state, double *out_663431017939020115) {
  H_mod_fun(state, out_663431017939020115);
}
void car_f_fun(double *state, double dt, double *out_2189594786826403342) {
  f_fun(state,  dt, out_2189594786826403342);
}
void car_F_fun(double *state, double dt, double *out_2972862917221800866) {
  F_fun(state,  dt, out_2972862917221800866);
}
void car_h_25(double *state, double *unused, double *out_5656616099551326110) {
  h_25(state, unused, out_5656616099551326110);
}
void car_H_25(double *state, double *unused, double *out_6913354757324836494) {
  H_25(state, unused, out_6913354757324836494);
}
void car_h_24(double *state, double *unused, double *out_4997502995021187120) {
  h_24(state, unused, out_4997502995021187120);
}
void car_H_24(double *state, double *unused, double *out_4297999924517583471) {
  H_24(state, unused, out_4297999924517583471);
}
void car_h_30(double *state, double *unused, double *out_4141240455864945865) {
  h_30(state, unused, out_4141240455864945865);
}
void car_H_30(double *state, double *unused, double *out_2385658427197228296) {
  H_30(state, unused, out_2385658427197228296);
}
void car_h_26(double *state, double *unused, double *out_860432811795899780) {
  h_26(state, unused, out_860432811795899780);
}
void car_H_26(double *state, double *unused, double *out_3171851438450780270) {
  H_26(state, unused, out_3171851438450780270);
}
void car_h_27(double *state, double *unused, double *out_7788454706497379861) {
  h_27(state, unused, out_7788454706497379861);
}
void car_H_27(double *state, double *unused, double *out_4609252498381171513) {
  H_27(state, unused, out_4609252498381171513);
}
void car_h_29(double *state, double *unused, double *out_4840728716141894900) {
  h_29(state, unused, out_4840728716141894900);
}
void car_H_29(double *state, double *unused, double *out_2895889771511620480) {
  H_29(state, unused, out_2895889771511620480);
}
void car_h_28(double *state, double *unused, double *out_2786781438836534216) {
  h_28(state, unused, out_2786781438836534216);
}
void car_H_28(double *state, double *unused, double *out_4859520043076946731) {
  H_28(state, unused, out_4859520043076946731);
}
void car_h_31(double *state, double *unused, double *out_1794932623122824458) {
  h_31(state, unused, out_1794932623122824458);
}
void car_H_31(double *state, double *unused, double *out_6944000719201796922) {
  H_31(state, unused, out_6944000719201796922);
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
