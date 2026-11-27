#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

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
 *                      Code generated with SymPy 1.14.0                      *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4080457982728834630) {
   out_4080457982728834630[0] = delta_x[0] + nom_x[0];
   out_4080457982728834630[1] = delta_x[1] + nom_x[1];
   out_4080457982728834630[2] = delta_x[2] + nom_x[2];
   out_4080457982728834630[3] = delta_x[3] + nom_x[3];
   out_4080457982728834630[4] = delta_x[4] + nom_x[4];
   out_4080457982728834630[5] = delta_x[5] + nom_x[5];
   out_4080457982728834630[6] = delta_x[6] + nom_x[6];
   out_4080457982728834630[7] = delta_x[7] + nom_x[7];
   out_4080457982728834630[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2831246477908734191) {
   out_2831246477908734191[0] = -nom_x[0] + true_x[0];
   out_2831246477908734191[1] = -nom_x[1] + true_x[1];
   out_2831246477908734191[2] = -nom_x[2] + true_x[2];
   out_2831246477908734191[3] = -nom_x[3] + true_x[3];
   out_2831246477908734191[4] = -nom_x[4] + true_x[4];
   out_2831246477908734191[5] = -nom_x[5] + true_x[5];
   out_2831246477908734191[6] = -nom_x[6] + true_x[6];
   out_2831246477908734191[7] = -nom_x[7] + true_x[7];
   out_2831246477908734191[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_7119319816853775953) {
   out_7119319816853775953[0] = 1.0;
   out_7119319816853775953[1] = 0.0;
   out_7119319816853775953[2] = 0.0;
   out_7119319816853775953[3] = 0.0;
   out_7119319816853775953[4] = 0.0;
   out_7119319816853775953[5] = 0.0;
   out_7119319816853775953[6] = 0.0;
   out_7119319816853775953[7] = 0.0;
   out_7119319816853775953[8] = 0.0;
   out_7119319816853775953[9] = 0.0;
   out_7119319816853775953[10] = 1.0;
   out_7119319816853775953[11] = 0.0;
   out_7119319816853775953[12] = 0.0;
   out_7119319816853775953[13] = 0.0;
   out_7119319816853775953[14] = 0.0;
   out_7119319816853775953[15] = 0.0;
   out_7119319816853775953[16] = 0.0;
   out_7119319816853775953[17] = 0.0;
   out_7119319816853775953[18] = 0.0;
   out_7119319816853775953[19] = 0.0;
   out_7119319816853775953[20] = 1.0;
   out_7119319816853775953[21] = 0.0;
   out_7119319816853775953[22] = 0.0;
   out_7119319816853775953[23] = 0.0;
   out_7119319816853775953[24] = 0.0;
   out_7119319816853775953[25] = 0.0;
   out_7119319816853775953[26] = 0.0;
   out_7119319816853775953[27] = 0.0;
   out_7119319816853775953[28] = 0.0;
   out_7119319816853775953[29] = 0.0;
   out_7119319816853775953[30] = 1.0;
   out_7119319816853775953[31] = 0.0;
   out_7119319816853775953[32] = 0.0;
   out_7119319816853775953[33] = 0.0;
   out_7119319816853775953[34] = 0.0;
   out_7119319816853775953[35] = 0.0;
   out_7119319816853775953[36] = 0.0;
   out_7119319816853775953[37] = 0.0;
   out_7119319816853775953[38] = 0.0;
   out_7119319816853775953[39] = 0.0;
   out_7119319816853775953[40] = 1.0;
   out_7119319816853775953[41] = 0.0;
   out_7119319816853775953[42] = 0.0;
   out_7119319816853775953[43] = 0.0;
   out_7119319816853775953[44] = 0.0;
   out_7119319816853775953[45] = 0.0;
   out_7119319816853775953[46] = 0.0;
   out_7119319816853775953[47] = 0.0;
   out_7119319816853775953[48] = 0.0;
   out_7119319816853775953[49] = 0.0;
   out_7119319816853775953[50] = 1.0;
   out_7119319816853775953[51] = 0.0;
   out_7119319816853775953[52] = 0.0;
   out_7119319816853775953[53] = 0.0;
   out_7119319816853775953[54] = 0.0;
   out_7119319816853775953[55] = 0.0;
   out_7119319816853775953[56] = 0.0;
   out_7119319816853775953[57] = 0.0;
   out_7119319816853775953[58] = 0.0;
   out_7119319816853775953[59] = 0.0;
   out_7119319816853775953[60] = 1.0;
   out_7119319816853775953[61] = 0.0;
   out_7119319816853775953[62] = 0.0;
   out_7119319816853775953[63] = 0.0;
   out_7119319816853775953[64] = 0.0;
   out_7119319816853775953[65] = 0.0;
   out_7119319816853775953[66] = 0.0;
   out_7119319816853775953[67] = 0.0;
   out_7119319816853775953[68] = 0.0;
   out_7119319816853775953[69] = 0.0;
   out_7119319816853775953[70] = 1.0;
   out_7119319816853775953[71] = 0.0;
   out_7119319816853775953[72] = 0.0;
   out_7119319816853775953[73] = 0.0;
   out_7119319816853775953[74] = 0.0;
   out_7119319816853775953[75] = 0.0;
   out_7119319816853775953[76] = 0.0;
   out_7119319816853775953[77] = 0.0;
   out_7119319816853775953[78] = 0.0;
   out_7119319816853775953[79] = 0.0;
   out_7119319816853775953[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_6290842336792412547) {
   out_6290842336792412547[0] = state[0];
   out_6290842336792412547[1] = state[1];
   out_6290842336792412547[2] = state[2];
   out_6290842336792412547[3] = state[3];
   out_6290842336792412547[4] = state[4];
   out_6290842336792412547[5] = dt*(stiffness_front*(-state[2] - state[3] + state[7])*state[0]/state[1] + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/state[4] + ((-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/state[4] - state[4])*state[6] - 9.8100000000000005*state[8]) + state[5];
   out_6290842336792412547[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6290842336792412547[7] = state[7];
   out_6290842336792412547[8] = state[8];
}
void F_fun(double *state, double dt, double *out_5917174526903133239) {
   out_5917174526903133239[0] = 1;
   out_5917174526903133239[1] = 0;
   out_5917174526903133239[2] = 0;
   out_5917174526903133239[3] = 0;
   out_5917174526903133239[4] = 0;
   out_5917174526903133239[5] = 0;
   out_5917174526903133239[6] = 0;
   out_5917174526903133239[7] = 0;
   out_5917174526903133239[8] = 0;
   out_5917174526903133239[9] = 0;
   out_5917174526903133239[10] = 1;
   out_5917174526903133239[11] = 0;
   out_5917174526903133239[12] = 0;
   out_5917174526903133239[13] = 0;
   out_5917174526903133239[14] = 0;
   out_5917174526903133239[15] = 0;
   out_5917174526903133239[16] = 0;
   out_5917174526903133239[17] = 0;
   out_5917174526903133239[18] = 0;
   out_5917174526903133239[19] = 0;
   out_5917174526903133239[20] = 1;
   out_5917174526903133239[21] = 0;
   out_5917174526903133239[22] = 0;
   out_5917174526903133239[23] = 0;
   out_5917174526903133239[24] = 0;
   out_5917174526903133239[25] = 0;
   out_5917174526903133239[26] = 0;
   out_5917174526903133239[27] = 0;
   out_5917174526903133239[28] = 0;
   out_5917174526903133239[29] = 0;
   out_5917174526903133239[30] = 1;
   out_5917174526903133239[31] = 0;
   out_5917174526903133239[32] = 0;
   out_5917174526903133239[33] = 0;
   out_5917174526903133239[34] = 0;
   out_5917174526903133239[35] = 0;
   out_5917174526903133239[36] = 0;
   out_5917174526903133239[37] = 0;
   out_5917174526903133239[38] = 0;
   out_5917174526903133239[39] = 0;
   out_5917174526903133239[40] = 1;
   out_5917174526903133239[41] = 0;
   out_5917174526903133239[42] = 0;
   out_5917174526903133239[43] = 0;
   out_5917174526903133239[44] = 0;
   out_5917174526903133239[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/state[1] + (-stiffness_front - stiffness_rear)*state[5]/state[4] + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/state[4]);
   out_5917174526903133239[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/pow(state[1], 2);
   out_5917174526903133239[47] = -dt*stiffness_front*state[0]/state[1];
   out_5917174526903133239[48] = -dt*stiffness_front*state[0]/state[1];
   out_5917174526903133239[49] = dt*(-(-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/pow(state[4], 2) + (-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/pow(state[4], 2) - 1)*state[6]);
   out_5917174526903133239[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/state[4] + 1;
   out_5917174526903133239[51] = dt*((-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/state[4] - state[4]);
   out_5917174526903133239[52] = dt*stiffness_front*state[0]/state[1];
   out_5917174526903133239[53] = -9.8100000000000005*dt;
   out_5917174526903133239[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5917174526903133239[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5917174526903133239[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5917174526903133239[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5917174526903133239[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5917174526903133239[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5917174526903133239[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5917174526903133239[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5917174526903133239[62] = 0;
   out_5917174526903133239[63] = 0;
   out_5917174526903133239[64] = 0;
   out_5917174526903133239[65] = 0;
   out_5917174526903133239[66] = 0;
   out_5917174526903133239[67] = 0;
   out_5917174526903133239[68] = 0;
   out_5917174526903133239[69] = 0;
   out_5917174526903133239[70] = 1;
   out_5917174526903133239[71] = 0;
   out_5917174526903133239[72] = 0;
   out_5917174526903133239[73] = 0;
   out_5917174526903133239[74] = 0;
   out_5917174526903133239[75] = 0;
   out_5917174526903133239[76] = 0;
   out_5917174526903133239[77] = 0;
   out_5917174526903133239[78] = 0;
   out_5917174526903133239[79] = 0;
   out_5917174526903133239[80] = 1;
}
void h_25(double *state, double *unused, double *out_6854990307410438286) {
   out_6854990307410438286[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6423732431671110823) {
   out_6423732431671110823[0] = 0;
   out_6423732431671110823[1] = 0;
   out_6423732431671110823[2] = 0;
   out_6423732431671110823[3] = 0;
   out_6423732431671110823[4] = 0;
   out_6423732431671110823[5] = 0;
   out_6423732431671110823[6] = 1;
   out_6423732431671110823[7] = 0;
   out_6423732431671110823[8] = 0;
}
void h_24(double *state, double *unused, double *out_851954178914977484) {
   out_851954178914977484[0] = state[4];
   out_851954178914977484[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8040688490640672895) {
   out_8040688490640672895[0] = 0;
   out_8040688490640672895[1] = 0;
   out_8040688490640672895[2] = 0;
   out_8040688490640672895[3] = 0;
   out_8040688490640672895[4] = 1;
   out_8040688490640672895[5] = 0;
   out_8040688490640672895[6] = 0;
   out_8040688490640672895[7] = 0;
   out_8040688490640672895[8] = 0;
   out_8040688490640672895[9] = 0;
   out_8040688490640672895[10] = 0;
   out_8040688490640672895[11] = 0;
   out_8040688490640672895[12] = 0;
   out_8040688490640672895[13] = 0;
   out_8040688490640672895[14] = 1;
   out_8040688490640672895[15] = 0;
   out_8040688490640672895[16] = 0;
   out_8040688490640672895[17] = 0;
}
void h_30(double *state, double *unused, double *out_3698814271414589971) {
   out_3698814271414589971[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3905399473163862196) {
   out_3905399473163862196[0] = 0;
   out_3905399473163862196[1] = 0;
   out_3905399473163862196[2] = 0;
   out_3905399473163862196[3] = 0;
   out_3905399473163862196[4] = 1;
   out_3905399473163862196[5] = 0;
   out_3905399473163862196[6] = 0;
   out_3905399473163862196[7] = 0;
   out_3905399473163862196[8] = 0;
}
void h_26(double *state, double *unused, double *out_6735600322858922656) {
   out_6735600322858922656[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3119206461910310222) {
   out_3119206461910310222[0] = 0;
   out_3119206461910310222[1] = 0;
   out_3119206461910310222[2] = 0;
   out_3119206461910310222[3] = 0;
   out_3119206461910310222[4] = 0;
   out_3119206461910310222[5] = 0;
   out_3119206461910310222[6] = 0;
   out_3119206461910310222[7] = 1;
   out_3119206461910310222[8] = 0;
}
void h_27(double *state, double *unused, double *out_608539249950264818) {
   out_608539249950264818[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1681805401979918979) {
   out_1681805401979918979[0] = 0;
   out_1681805401979918979[1] = 0;
   out_1681805401979918979[2] = 0;
   out_1681805401979918979[3] = 1;
   out_1681805401979918979[4] = 0;
   out_1681805401979918979[5] = 0;
   out_1681805401979918979[6] = 0;
   out_1681805401979918979[7] = 0;
   out_1681805401979918979[8] = 0;
}
void h_29(double *state, double *unused, double *out_2008366436678909011) {
   out_2008366436678909011[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3395168128849470012) {
   out_3395168128849470012[0] = 0;
   out_3395168128849470012[1] = 1;
   out_3395168128849470012[2] = 0;
   out_3395168128849470012[3] = 0;
   out_3395168128849470012[4] = 0;
   out_3395168128849470012[5] = 0;
   out_3395168128849470012[6] = 0;
   out_3395168128849470012[7] = 0;
   out_3395168128849470012[8] = 0;
}
void h_28(double *state, double *unused, double *out_7611457861953118210) {
   out_7611457861953118210[0] = state[0];
}
void H_28(double *state, double *unused, double *out_8477567145919000586) {
   out_8477567145919000586[0] = 1;
   out_8477567145919000586[1] = 0;
   out_8477567145919000586[2] = 0;
   out_8477567145919000586[3] = 0;
   out_8477567145919000586[4] = 0;
   out_8477567145919000586[5] = 0;
   out_8477567145919000586[6] = 0;
   out_8477567145919000586[7] = 0;
   out_8477567145919000586[8] = 0;
}
void h_31(double *state, double *unused, double *out_7130184369694944175) {
   out_7130184369694944175[0] = state[8];
}
void H_31(double *state, double *unused, double *out_652942818840706430) {
   out_652942818840706430[0] = 0;
   out_652942818840706430[1] = 0;
   out_652942818840706430[2] = 0;
   out_652942818840706430[3] = 0;
   out_652942818840706430[4] = 0;
   out_652942818840706430[5] = 0;
   out_652942818840706430[6] = 0;
   out_652942818840706430[7] = 0;
   out_652942818840706430[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_4080457982728834630) {
  err_fun(nom_x, delta_x, out_4080457982728834630);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2831246477908734191) {
  inv_err_fun(nom_x, true_x, out_2831246477908734191);
}
void car_H_mod_fun(double *state, double *out_7119319816853775953) {
  H_mod_fun(state, out_7119319816853775953);
}
void car_f_fun(double *state, double dt, double *out_6290842336792412547) {
  f_fun(state,  dt, out_6290842336792412547);
}
void car_F_fun(double *state, double dt, double *out_5917174526903133239) {
  F_fun(state,  dt, out_5917174526903133239);
}
void car_h_25(double *state, double *unused, double *out_6854990307410438286) {
  h_25(state, unused, out_6854990307410438286);
}
void car_H_25(double *state, double *unused, double *out_6423732431671110823) {
  H_25(state, unused, out_6423732431671110823);
}
void car_h_24(double *state, double *unused, double *out_851954178914977484) {
  h_24(state, unused, out_851954178914977484);
}
void car_H_24(double *state, double *unused, double *out_8040688490640672895) {
  H_24(state, unused, out_8040688490640672895);
}
void car_h_30(double *state, double *unused, double *out_3698814271414589971) {
  h_30(state, unused, out_3698814271414589971);
}
void car_H_30(double *state, double *unused, double *out_3905399473163862196) {
  H_30(state, unused, out_3905399473163862196);
}
void car_h_26(double *state, double *unused, double *out_6735600322858922656) {
  h_26(state, unused, out_6735600322858922656);
}
void car_H_26(double *state, double *unused, double *out_3119206461910310222) {
  H_26(state, unused, out_3119206461910310222);
}
void car_h_27(double *state, double *unused, double *out_608539249950264818) {
  h_27(state, unused, out_608539249950264818);
}
void car_H_27(double *state, double *unused, double *out_1681805401979918979) {
  H_27(state, unused, out_1681805401979918979);
}
void car_h_29(double *state, double *unused, double *out_2008366436678909011) {
  h_29(state, unused, out_2008366436678909011);
}
void car_H_29(double *state, double *unused, double *out_3395168128849470012) {
  H_29(state, unused, out_3395168128849470012);
}
void car_h_28(double *state, double *unused, double *out_7611457861953118210) {
  h_28(state, unused, out_7611457861953118210);
}
void car_H_28(double *state, double *unused, double *out_8477567145919000586) {
  H_28(state, unused, out_8477567145919000586);
}
void car_h_31(double *state, double *unused, double *out_7130184369694944175) {
  h_31(state, unused, out_7130184369694944175);
}
void car_H_31(double *state, double *unused, double *out_652942818840706430) {
  H_31(state, unused, out_652942818840706430);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
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
