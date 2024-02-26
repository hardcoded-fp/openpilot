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
void car_err_fun(double *nom_x, double *delta_x, double *out_6288688874102308164);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5232615150247248436);
void car_H_mod_fun(double *state, double *out_4159914697228368279);
void car_f_fun(double *state, double dt, double *out_4156207149520197820);
void car_F_fun(double *state, double dt, double *out_1014683453205324692);
void car_h_25(double *state, double *unused, double *out_5671253352182223489);
void car_H_25(double *state, double *unused, double *out_2927558261712726917);
void car_h_24(double *state, double *unused, double *out_3382553279571021005);
void car_H_24(double *state, double *unused, double *out_7674067588679235951);
void car_h_30(double *state, double *unused, double *out_6004655495176977191);
void car_H_30(double *state, double *unused, double *out_3989132079778889838);
void car_h_26(double *state, double *unused, double *out_1002982227516131546);
void car_H_26(double *state, double *unused, double *out_6669061580586783141);
void car_h_27(double *state, double *unused, double *out_7361701186917904449);
void car_H_27(double *state, double *unused, double *out_1814368767978464927);
void car_h_29(double *state, double *unused, double *out_2401001324467249462);
void car_H_29(double *state, double *unused, double *out_101006041108913894);
void car_h_28(double *state, double *unused, double *out_2181546259182111216);
void car_H_28(double *state, double *unused, double *out_4981392975960616680);
void car_h_31(double *state, double *unused, double *out_425985762948188827);
void car_H_31(double *state, double *unused, double *out_2896912299835766489);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}