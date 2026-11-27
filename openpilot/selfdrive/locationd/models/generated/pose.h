#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void pose_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void pose_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void pose_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void pose_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void pose_err_fun(double *nom_x, double *delta_x, double *out_4750261472376034529);
void pose_inv_err_fun(double *nom_x, double *true_x, double *out_8277444623959904026);
void pose_H_mod_fun(double *state, double *out_7806277314590803126);
void pose_f_fun(double *state, double dt, double *out_4548647388261971274);
void pose_F_fun(double *state, double dt, double *out_1059548545756754699);
void pose_h_4(double *state, double *unused, double *out_1272267543225290969);
void pose_H_4(double *state, double *unused, double *out_4532055235799853873);
void pose_h_10(double *state, double *unused, double *out_4823025017146951895);
void pose_H_10(double *state, double *unused, double *out_7328110893683068042);
void pose_h_13(double *state, double *unused, double *out_3553628753437543153);
void pose_H_13(double *state, double *unused, double *out_7744329061132186674);
void pose_h_14(double *state, double *unused, double *out_2142424780536854880);
void pose_H_14(double *state, double *unused, double *out_8495296092139338402);
void pose_predict(double *in_x, double *in_P, double *in_Q, double dt);
}