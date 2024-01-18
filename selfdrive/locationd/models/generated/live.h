#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_60673759599383242);
void live_err_fun(double *nom_x, double *delta_x, double *out_115236777840615556);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_2198698251275647311);
void live_H_mod_fun(double *state, double *out_8103795288019997415);
void live_f_fun(double *state, double dt, double *out_3618810602285908789);
void live_F_fun(double *state, double dt, double *out_812707499149624584);
void live_h_4(double *state, double *unused, double *out_519059569381175517);
void live_H_4(double *state, double *unused, double *out_6993697842471443277);
void live_h_9(double *state, double *unused, double *out_183230378108869121);
void live_H_9(double *state, double *unused, double *out_4165827295973660869);
void live_h_10(double *state, double *unused, double *out_7103314064859147806);
void live_H_10(double *state, double *unused, double *out_181541369712314370);
void live_h_12(double *state, double *unused, double *out_3105114867520165171);
void live_H_12(double *state, double *unused, double *out_3785917917555657847);
void live_h_35(double *state, double *unused, double *out_2361744896182109929);
void live_H_35(double *state, double *unused, double *out_1040354885230644138);
void live_h_32(double *state, double *unused, double *out_9034153276957235957);
void live_H_32(double *state, double *unused, double *out_7955219388924760494);
void live_h_13(double *state, double *unused, double *out_1821120093556423515);
void live_H_13(double *state, double *unused, double *out_580227250654406660);
void live_h_14(double *state, double *unused, double *out_183230378108869121);
void live_H_14(double *state, double *unused, double *out_4165827295973660869);
void live_h_33(double *state, double *unused, double *out_739960791988425354);
void live_H_33(double *state, double *unused, double *out_2110202119408213466);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}