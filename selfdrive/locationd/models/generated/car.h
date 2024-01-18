#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3962033384618978572);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8377453335461740860);
void car_H_mod_fun(double *state, double *out_663431017939020115);
void car_f_fun(double *state, double dt, double *out_2189594786826403342);
void car_F_fun(double *state, double dt, double *out_2972862917221800866);
void car_h_25(double *state, double *unused, double *out_5656616099551326110);
void car_H_25(double *state, double *unused, double *out_6913354757324836494);
void car_h_24(double *state, double *unused, double *out_4997502995021187120);
void car_H_24(double *state, double *unused, double *out_4297999924517583471);
void car_h_30(double *state, double *unused, double *out_4141240455864945865);
void car_H_30(double *state, double *unused, double *out_2385658427197228296);
void car_h_26(double *state, double *unused, double *out_860432811795899780);
void car_H_26(double *state, double *unused, double *out_3171851438450780270);
void car_h_27(double *state, double *unused, double *out_7788454706497379861);
void car_H_27(double *state, double *unused, double *out_4609252498381171513);
void car_h_29(double *state, double *unused, double *out_4840728716141894900);
void car_H_29(double *state, double *unused, double *out_2895889771511620480);
void car_h_28(double *state, double *unused, double *out_2786781438836534216);
void car_H_28(double *state, double *unused, double *out_4859520043076946731);
void car_h_31(double *state, double *unused, double *out_1794932623122824458);
void car_H_31(double *state, double *unused, double *out_6944000719201796922);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}