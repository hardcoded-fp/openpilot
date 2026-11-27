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
void car_err_fun(double *nom_x, double *delta_x, double *out_4080457982728834630);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_2831246477908734191);
void car_H_mod_fun(double *state, double *out_7119319816853775953);
void car_f_fun(double *state, double dt, double *out_6290842336792412547);
void car_F_fun(double *state, double dt, double *out_5917174526903133239);
void car_h_25(double *state, double *unused, double *out_6854990307410438286);
void car_H_25(double *state, double *unused, double *out_6423732431671110823);
void car_h_24(double *state, double *unused, double *out_851954178914977484);
void car_H_24(double *state, double *unused, double *out_8040688490640672895);
void car_h_30(double *state, double *unused, double *out_3698814271414589971);
void car_H_30(double *state, double *unused, double *out_3905399473163862196);
void car_h_26(double *state, double *unused, double *out_6735600322858922656);
void car_H_26(double *state, double *unused, double *out_3119206461910310222);
void car_h_27(double *state, double *unused, double *out_608539249950264818);
void car_H_27(double *state, double *unused, double *out_1681805401979918979);
void car_h_29(double *state, double *unused, double *out_2008366436678909011);
void car_H_29(double *state, double *unused, double *out_3395168128849470012);
void car_h_28(double *state, double *unused, double *out_7611457861953118210);
void car_H_28(double *state, double *unused, double *out_8477567145919000586);
void car_h_31(double *state, double *unused, double *out_7130184369694944175);
void car_H_31(double *state, double *unused, double *out_652942818840706430);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}