#include "pose.h"

namespace {
#define DIM 18
#define EDIM 18
#define MEDIM 18
typedef void (*Hfun)(double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814727903251177;
const static double MAHA_THRESH_10 = 7.814727903251177;
const static double MAHA_THRESH_13 = 7.814727903251177;
const static double MAHA_THRESH_14 = 7.814727903251177;

/******************************************************************************
 *                      Code generated with SymPy 1.14.0                      *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4750261472376034529) {
   out_4750261472376034529[0] = delta_x[0] + nom_x[0];
   out_4750261472376034529[1] = delta_x[1] + nom_x[1];
   out_4750261472376034529[2] = delta_x[2] + nom_x[2];
   out_4750261472376034529[3] = delta_x[3] + nom_x[3];
   out_4750261472376034529[4] = delta_x[4] + nom_x[4];
   out_4750261472376034529[5] = delta_x[5] + nom_x[5];
   out_4750261472376034529[6] = delta_x[6] + nom_x[6];
   out_4750261472376034529[7] = delta_x[7] + nom_x[7];
   out_4750261472376034529[8] = delta_x[8] + nom_x[8];
   out_4750261472376034529[9] = delta_x[9] + nom_x[9];
   out_4750261472376034529[10] = delta_x[10] + nom_x[10];
   out_4750261472376034529[11] = delta_x[11] + nom_x[11];
   out_4750261472376034529[12] = delta_x[12] + nom_x[12];
   out_4750261472376034529[13] = delta_x[13] + nom_x[13];
   out_4750261472376034529[14] = delta_x[14] + nom_x[14];
   out_4750261472376034529[15] = delta_x[15] + nom_x[15];
   out_4750261472376034529[16] = delta_x[16] + nom_x[16];
   out_4750261472376034529[17] = delta_x[17] + nom_x[17];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8277444623959904026) {
   out_8277444623959904026[0] = -nom_x[0] + true_x[0];
   out_8277444623959904026[1] = -nom_x[1] + true_x[1];
   out_8277444623959904026[2] = -nom_x[2] + true_x[2];
   out_8277444623959904026[3] = -nom_x[3] + true_x[3];
   out_8277444623959904026[4] = -nom_x[4] + true_x[4];
   out_8277444623959904026[5] = -nom_x[5] + true_x[5];
   out_8277444623959904026[6] = -nom_x[6] + true_x[6];
   out_8277444623959904026[7] = -nom_x[7] + true_x[7];
   out_8277444623959904026[8] = -nom_x[8] + true_x[8];
   out_8277444623959904026[9] = -nom_x[9] + true_x[9];
   out_8277444623959904026[10] = -nom_x[10] + true_x[10];
   out_8277444623959904026[11] = -nom_x[11] + true_x[11];
   out_8277444623959904026[12] = -nom_x[12] + true_x[12];
   out_8277444623959904026[13] = -nom_x[13] + true_x[13];
   out_8277444623959904026[14] = -nom_x[14] + true_x[14];
   out_8277444623959904026[15] = -nom_x[15] + true_x[15];
   out_8277444623959904026[16] = -nom_x[16] + true_x[16];
   out_8277444623959904026[17] = -nom_x[17] + true_x[17];
}
void H_mod_fun(double *state, double *out_7806277314590803126) {
   out_7806277314590803126[0] = 1.0;
   out_7806277314590803126[1] = 0.0;
   out_7806277314590803126[2] = 0.0;
   out_7806277314590803126[3] = 0.0;
   out_7806277314590803126[4] = 0.0;
   out_7806277314590803126[5] = 0.0;
   out_7806277314590803126[6] = 0.0;
   out_7806277314590803126[7] = 0.0;
   out_7806277314590803126[8] = 0.0;
   out_7806277314590803126[9] = 0.0;
   out_7806277314590803126[10] = 0.0;
   out_7806277314590803126[11] = 0.0;
   out_7806277314590803126[12] = 0.0;
   out_7806277314590803126[13] = 0.0;
   out_7806277314590803126[14] = 0.0;
   out_7806277314590803126[15] = 0.0;
   out_7806277314590803126[16] = 0.0;
   out_7806277314590803126[17] = 0.0;
   out_7806277314590803126[18] = 0.0;
   out_7806277314590803126[19] = 1.0;
   out_7806277314590803126[20] = 0.0;
   out_7806277314590803126[21] = 0.0;
   out_7806277314590803126[22] = 0.0;
   out_7806277314590803126[23] = 0.0;
   out_7806277314590803126[24] = 0.0;
   out_7806277314590803126[25] = 0.0;
   out_7806277314590803126[26] = 0.0;
   out_7806277314590803126[27] = 0.0;
   out_7806277314590803126[28] = 0.0;
   out_7806277314590803126[29] = 0.0;
   out_7806277314590803126[30] = 0.0;
   out_7806277314590803126[31] = 0.0;
   out_7806277314590803126[32] = 0.0;
   out_7806277314590803126[33] = 0.0;
   out_7806277314590803126[34] = 0.0;
   out_7806277314590803126[35] = 0.0;
   out_7806277314590803126[36] = 0.0;
   out_7806277314590803126[37] = 0.0;
   out_7806277314590803126[38] = 1.0;
   out_7806277314590803126[39] = 0.0;
   out_7806277314590803126[40] = 0.0;
   out_7806277314590803126[41] = 0.0;
   out_7806277314590803126[42] = 0.0;
   out_7806277314590803126[43] = 0.0;
   out_7806277314590803126[44] = 0.0;
   out_7806277314590803126[45] = 0.0;
   out_7806277314590803126[46] = 0.0;
   out_7806277314590803126[47] = 0.0;
   out_7806277314590803126[48] = 0.0;
   out_7806277314590803126[49] = 0.0;
   out_7806277314590803126[50] = 0.0;
   out_7806277314590803126[51] = 0.0;
   out_7806277314590803126[52] = 0.0;
   out_7806277314590803126[53] = 0.0;
   out_7806277314590803126[54] = 0.0;
   out_7806277314590803126[55] = 0.0;
   out_7806277314590803126[56] = 0.0;
   out_7806277314590803126[57] = 1.0;
   out_7806277314590803126[58] = 0.0;
   out_7806277314590803126[59] = 0.0;
   out_7806277314590803126[60] = 0.0;
   out_7806277314590803126[61] = 0.0;
   out_7806277314590803126[62] = 0.0;
   out_7806277314590803126[63] = 0.0;
   out_7806277314590803126[64] = 0.0;
   out_7806277314590803126[65] = 0.0;
   out_7806277314590803126[66] = 0.0;
   out_7806277314590803126[67] = 0.0;
   out_7806277314590803126[68] = 0.0;
   out_7806277314590803126[69] = 0.0;
   out_7806277314590803126[70] = 0.0;
   out_7806277314590803126[71] = 0.0;
   out_7806277314590803126[72] = 0.0;
   out_7806277314590803126[73] = 0.0;
   out_7806277314590803126[74] = 0.0;
   out_7806277314590803126[75] = 0.0;
   out_7806277314590803126[76] = 1.0;
   out_7806277314590803126[77] = 0.0;
   out_7806277314590803126[78] = 0.0;
   out_7806277314590803126[79] = 0.0;
   out_7806277314590803126[80] = 0.0;
   out_7806277314590803126[81] = 0.0;
   out_7806277314590803126[82] = 0.0;
   out_7806277314590803126[83] = 0.0;
   out_7806277314590803126[84] = 0.0;
   out_7806277314590803126[85] = 0.0;
   out_7806277314590803126[86] = 0.0;
   out_7806277314590803126[87] = 0.0;
   out_7806277314590803126[88] = 0.0;
   out_7806277314590803126[89] = 0.0;
   out_7806277314590803126[90] = 0.0;
   out_7806277314590803126[91] = 0.0;
   out_7806277314590803126[92] = 0.0;
   out_7806277314590803126[93] = 0.0;
   out_7806277314590803126[94] = 0.0;
   out_7806277314590803126[95] = 1.0;
   out_7806277314590803126[96] = 0.0;
   out_7806277314590803126[97] = 0.0;
   out_7806277314590803126[98] = 0.0;
   out_7806277314590803126[99] = 0.0;
   out_7806277314590803126[100] = 0.0;
   out_7806277314590803126[101] = 0.0;
   out_7806277314590803126[102] = 0.0;
   out_7806277314590803126[103] = 0.0;
   out_7806277314590803126[104] = 0.0;
   out_7806277314590803126[105] = 0.0;
   out_7806277314590803126[106] = 0.0;
   out_7806277314590803126[107] = 0.0;
   out_7806277314590803126[108] = 0.0;
   out_7806277314590803126[109] = 0.0;
   out_7806277314590803126[110] = 0.0;
   out_7806277314590803126[111] = 0.0;
   out_7806277314590803126[112] = 0.0;
   out_7806277314590803126[113] = 0.0;
   out_7806277314590803126[114] = 1.0;
   out_7806277314590803126[115] = 0.0;
   out_7806277314590803126[116] = 0.0;
   out_7806277314590803126[117] = 0.0;
   out_7806277314590803126[118] = 0.0;
   out_7806277314590803126[119] = 0.0;
   out_7806277314590803126[120] = 0.0;
   out_7806277314590803126[121] = 0.0;
   out_7806277314590803126[122] = 0.0;
   out_7806277314590803126[123] = 0.0;
   out_7806277314590803126[124] = 0.0;
   out_7806277314590803126[125] = 0.0;
   out_7806277314590803126[126] = 0.0;
   out_7806277314590803126[127] = 0.0;
   out_7806277314590803126[128] = 0.0;
   out_7806277314590803126[129] = 0.0;
   out_7806277314590803126[130] = 0.0;
   out_7806277314590803126[131] = 0.0;
   out_7806277314590803126[132] = 0.0;
   out_7806277314590803126[133] = 1.0;
   out_7806277314590803126[134] = 0.0;
   out_7806277314590803126[135] = 0.0;
   out_7806277314590803126[136] = 0.0;
   out_7806277314590803126[137] = 0.0;
   out_7806277314590803126[138] = 0.0;
   out_7806277314590803126[139] = 0.0;
   out_7806277314590803126[140] = 0.0;
   out_7806277314590803126[141] = 0.0;
   out_7806277314590803126[142] = 0.0;
   out_7806277314590803126[143] = 0.0;
   out_7806277314590803126[144] = 0.0;
   out_7806277314590803126[145] = 0.0;
   out_7806277314590803126[146] = 0.0;
   out_7806277314590803126[147] = 0.0;
   out_7806277314590803126[148] = 0.0;
   out_7806277314590803126[149] = 0.0;
   out_7806277314590803126[150] = 0.0;
   out_7806277314590803126[151] = 0.0;
   out_7806277314590803126[152] = 1.0;
   out_7806277314590803126[153] = 0.0;
   out_7806277314590803126[154] = 0.0;
   out_7806277314590803126[155] = 0.0;
   out_7806277314590803126[156] = 0.0;
   out_7806277314590803126[157] = 0.0;
   out_7806277314590803126[158] = 0.0;
   out_7806277314590803126[159] = 0.0;
   out_7806277314590803126[160] = 0.0;
   out_7806277314590803126[161] = 0.0;
   out_7806277314590803126[162] = 0.0;
   out_7806277314590803126[163] = 0.0;
   out_7806277314590803126[164] = 0.0;
   out_7806277314590803126[165] = 0.0;
   out_7806277314590803126[166] = 0.0;
   out_7806277314590803126[167] = 0.0;
   out_7806277314590803126[168] = 0.0;
   out_7806277314590803126[169] = 0.0;
   out_7806277314590803126[170] = 0.0;
   out_7806277314590803126[171] = 1.0;
   out_7806277314590803126[172] = 0.0;
   out_7806277314590803126[173] = 0.0;
   out_7806277314590803126[174] = 0.0;
   out_7806277314590803126[175] = 0.0;
   out_7806277314590803126[176] = 0.0;
   out_7806277314590803126[177] = 0.0;
   out_7806277314590803126[178] = 0.0;
   out_7806277314590803126[179] = 0.0;
   out_7806277314590803126[180] = 0.0;
   out_7806277314590803126[181] = 0.0;
   out_7806277314590803126[182] = 0.0;
   out_7806277314590803126[183] = 0.0;
   out_7806277314590803126[184] = 0.0;
   out_7806277314590803126[185] = 0.0;
   out_7806277314590803126[186] = 0.0;
   out_7806277314590803126[187] = 0.0;
   out_7806277314590803126[188] = 0.0;
   out_7806277314590803126[189] = 0.0;
   out_7806277314590803126[190] = 1.0;
   out_7806277314590803126[191] = 0.0;
   out_7806277314590803126[192] = 0.0;
   out_7806277314590803126[193] = 0.0;
   out_7806277314590803126[194] = 0.0;
   out_7806277314590803126[195] = 0.0;
   out_7806277314590803126[196] = 0.0;
   out_7806277314590803126[197] = 0.0;
   out_7806277314590803126[198] = 0.0;
   out_7806277314590803126[199] = 0.0;
   out_7806277314590803126[200] = 0.0;
   out_7806277314590803126[201] = 0.0;
   out_7806277314590803126[202] = 0.0;
   out_7806277314590803126[203] = 0.0;
   out_7806277314590803126[204] = 0.0;
   out_7806277314590803126[205] = 0.0;
   out_7806277314590803126[206] = 0.0;
   out_7806277314590803126[207] = 0.0;
   out_7806277314590803126[208] = 0.0;
   out_7806277314590803126[209] = 1.0;
   out_7806277314590803126[210] = 0.0;
   out_7806277314590803126[211] = 0.0;
   out_7806277314590803126[212] = 0.0;
   out_7806277314590803126[213] = 0.0;
   out_7806277314590803126[214] = 0.0;
   out_7806277314590803126[215] = 0.0;
   out_7806277314590803126[216] = 0.0;
   out_7806277314590803126[217] = 0.0;
   out_7806277314590803126[218] = 0.0;
   out_7806277314590803126[219] = 0.0;
   out_7806277314590803126[220] = 0.0;
   out_7806277314590803126[221] = 0.0;
   out_7806277314590803126[222] = 0.0;
   out_7806277314590803126[223] = 0.0;
   out_7806277314590803126[224] = 0.0;
   out_7806277314590803126[225] = 0.0;
   out_7806277314590803126[226] = 0.0;
   out_7806277314590803126[227] = 0.0;
   out_7806277314590803126[228] = 1.0;
   out_7806277314590803126[229] = 0.0;
   out_7806277314590803126[230] = 0.0;
   out_7806277314590803126[231] = 0.0;
   out_7806277314590803126[232] = 0.0;
   out_7806277314590803126[233] = 0.0;
   out_7806277314590803126[234] = 0.0;
   out_7806277314590803126[235] = 0.0;
   out_7806277314590803126[236] = 0.0;
   out_7806277314590803126[237] = 0.0;
   out_7806277314590803126[238] = 0.0;
   out_7806277314590803126[239] = 0.0;
   out_7806277314590803126[240] = 0.0;
   out_7806277314590803126[241] = 0.0;
   out_7806277314590803126[242] = 0.0;
   out_7806277314590803126[243] = 0.0;
   out_7806277314590803126[244] = 0.0;
   out_7806277314590803126[245] = 0.0;
   out_7806277314590803126[246] = 0.0;
   out_7806277314590803126[247] = 1.0;
   out_7806277314590803126[248] = 0.0;
   out_7806277314590803126[249] = 0.0;
   out_7806277314590803126[250] = 0.0;
   out_7806277314590803126[251] = 0.0;
   out_7806277314590803126[252] = 0.0;
   out_7806277314590803126[253] = 0.0;
   out_7806277314590803126[254] = 0.0;
   out_7806277314590803126[255] = 0.0;
   out_7806277314590803126[256] = 0.0;
   out_7806277314590803126[257] = 0.0;
   out_7806277314590803126[258] = 0.0;
   out_7806277314590803126[259] = 0.0;
   out_7806277314590803126[260] = 0.0;
   out_7806277314590803126[261] = 0.0;
   out_7806277314590803126[262] = 0.0;
   out_7806277314590803126[263] = 0.0;
   out_7806277314590803126[264] = 0.0;
   out_7806277314590803126[265] = 0.0;
   out_7806277314590803126[266] = 1.0;
   out_7806277314590803126[267] = 0.0;
   out_7806277314590803126[268] = 0.0;
   out_7806277314590803126[269] = 0.0;
   out_7806277314590803126[270] = 0.0;
   out_7806277314590803126[271] = 0.0;
   out_7806277314590803126[272] = 0.0;
   out_7806277314590803126[273] = 0.0;
   out_7806277314590803126[274] = 0.0;
   out_7806277314590803126[275] = 0.0;
   out_7806277314590803126[276] = 0.0;
   out_7806277314590803126[277] = 0.0;
   out_7806277314590803126[278] = 0.0;
   out_7806277314590803126[279] = 0.0;
   out_7806277314590803126[280] = 0.0;
   out_7806277314590803126[281] = 0.0;
   out_7806277314590803126[282] = 0.0;
   out_7806277314590803126[283] = 0.0;
   out_7806277314590803126[284] = 0.0;
   out_7806277314590803126[285] = 1.0;
   out_7806277314590803126[286] = 0.0;
   out_7806277314590803126[287] = 0.0;
   out_7806277314590803126[288] = 0.0;
   out_7806277314590803126[289] = 0.0;
   out_7806277314590803126[290] = 0.0;
   out_7806277314590803126[291] = 0.0;
   out_7806277314590803126[292] = 0.0;
   out_7806277314590803126[293] = 0.0;
   out_7806277314590803126[294] = 0.0;
   out_7806277314590803126[295] = 0.0;
   out_7806277314590803126[296] = 0.0;
   out_7806277314590803126[297] = 0.0;
   out_7806277314590803126[298] = 0.0;
   out_7806277314590803126[299] = 0.0;
   out_7806277314590803126[300] = 0.0;
   out_7806277314590803126[301] = 0.0;
   out_7806277314590803126[302] = 0.0;
   out_7806277314590803126[303] = 0.0;
   out_7806277314590803126[304] = 1.0;
   out_7806277314590803126[305] = 0.0;
   out_7806277314590803126[306] = 0.0;
   out_7806277314590803126[307] = 0.0;
   out_7806277314590803126[308] = 0.0;
   out_7806277314590803126[309] = 0.0;
   out_7806277314590803126[310] = 0.0;
   out_7806277314590803126[311] = 0.0;
   out_7806277314590803126[312] = 0.0;
   out_7806277314590803126[313] = 0.0;
   out_7806277314590803126[314] = 0.0;
   out_7806277314590803126[315] = 0.0;
   out_7806277314590803126[316] = 0.0;
   out_7806277314590803126[317] = 0.0;
   out_7806277314590803126[318] = 0.0;
   out_7806277314590803126[319] = 0.0;
   out_7806277314590803126[320] = 0.0;
   out_7806277314590803126[321] = 0.0;
   out_7806277314590803126[322] = 0.0;
   out_7806277314590803126[323] = 1.0;
}
void f_fun(double *state, double dt, double *out_4548647388261971274) {
   out_4548647388261971274[0] = atan2((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), -(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]));
   out_4548647388261971274[1] = asin(sin(dt*state[7])*cos(state[0])*cos(state[1]) - sin(dt*state[8])*sin(state[0])*cos(dt*state[7])*cos(state[1]) + sin(state[1])*cos(dt*state[7])*cos(dt*state[8]));
   out_4548647388261971274[2] = atan2(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), -(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]));
   out_4548647388261971274[3] = dt*state[12] + state[3];
   out_4548647388261971274[4] = dt*state[13] + state[4];
   out_4548647388261971274[5] = dt*state[14] + state[5];
   out_4548647388261971274[6] = state[6];
   out_4548647388261971274[7] = state[7];
   out_4548647388261971274[8] = state[8];
   out_4548647388261971274[9] = state[9];
   out_4548647388261971274[10] = state[10];
   out_4548647388261971274[11] = state[11];
   out_4548647388261971274[12] = state[12];
   out_4548647388261971274[13] = state[13];
   out_4548647388261971274[14] = state[14];
   out_4548647388261971274[15] = state[15];
   out_4548647388261971274[16] = state[16];
   out_4548647388261971274[17] = state[17];
}
void F_fun(double *state, double dt, double *out_1059548545756754699) {
   out_1059548545756754699[0] = ((-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*cos(state[0])*cos(state[1]) - sin(state[0])*cos(dt*state[6])*cos(dt*state[7])*cos(state[1]))*(-(sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) + (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) - sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2)) + ((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*cos(state[0])*cos(state[1]) - sin(dt*state[6])*sin(state[0])*cos(dt*state[7])*cos(state[1]))*(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2));
   out_1059548545756754699[1] = ((-sin(dt*state[6])*sin(dt*state[8]) - sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*cos(state[1]) - (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*sin(state[1]) - sin(state[1])*cos(dt*state[6])*cos(dt*state[7])*cos(state[0]))*(-(sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) + (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) - sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2)) + (-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))*(-(sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*sin(state[1]) + (-sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) + sin(dt*state[8])*cos(dt*state[6]))*cos(state[1]) - sin(dt*state[6])*sin(state[1])*cos(dt*state[7])*cos(state[0]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2));
   out_1059548545756754699[2] = 0;
   out_1059548545756754699[3] = 0;
   out_1059548545756754699[4] = 0;
   out_1059548545756754699[5] = 0;
   out_1059548545756754699[6] = (-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))*(dt*cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]) + (-dt*sin(dt*state[6])*sin(dt*state[8]) - dt*sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-dt*sin(dt*state[6])*cos(dt*state[8]) + dt*sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2)) + (-(sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) + (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) - sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))*(-dt*sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]) + (-dt*sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) - dt*cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) + (dt*sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - dt*sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2));
   out_1059548545756754699[7] = (-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))*(-dt*sin(dt*state[6])*sin(dt*state[7])*cos(state[0])*cos(state[1]) + dt*sin(dt*state[6])*sin(dt*state[8])*sin(state[0])*cos(dt*state[7])*cos(state[1]) - dt*sin(dt*state[6])*sin(state[1])*cos(dt*state[7])*cos(dt*state[8]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2)) + (-(sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) + (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) - sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))*(-dt*sin(dt*state[7])*cos(dt*state[6])*cos(state[0])*cos(state[1]) + dt*sin(dt*state[8])*sin(state[0])*cos(dt*state[6])*cos(dt*state[7])*cos(state[1]) - dt*sin(state[1])*cos(dt*state[6])*cos(dt*state[7])*cos(dt*state[8]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2));
   out_1059548545756754699[8] = ((dt*sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + dt*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (dt*sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - dt*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]))*(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2)) + ((dt*sin(dt*state[6])*sin(dt*state[8]) + dt*sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) + (-dt*sin(dt*state[6])*cos(dt*state[8]) + dt*sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]))*(-(sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) + (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) - sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]))/(pow(-(sin(dt*state[6])*sin(dt*state[8]) + sin(dt*state[7])*cos(dt*state[6])*cos(dt*state[8]))*sin(state[1]) + (-sin(dt*state[6])*cos(dt*state[8]) + sin(dt*state[7])*sin(dt*state[8])*cos(dt*state[6]))*sin(state[0])*cos(state[1]) + cos(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2) + pow((sin(dt*state[6])*sin(dt*state[7])*sin(dt*state[8]) + cos(dt*state[6])*cos(dt*state[8]))*sin(state[0])*cos(state[1]) - (sin(dt*state[6])*sin(dt*state[7])*cos(dt*state[8]) - sin(dt*state[8])*cos(dt*state[6]))*sin(state[1]) + sin(dt*state[6])*cos(dt*state[7])*cos(state[0])*cos(state[1]), 2));
   out_1059548545756754699[9] = 0;
   out_1059548545756754699[10] = 0;
   out_1059548545756754699[11] = 0;
   out_1059548545756754699[12] = 0;
   out_1059548545756754699[13] = 0;
   out_1059548545756754699[14] = 0;
   out_1059548545756754699[15] = 0;
   out_1059548545756754699[16] = 0;
   out_1059548545756754699[17] = 0;
   out_1059548545756754699[18] = (-sin(dt*state[7])*sin(state[0])*cos(state[1]) - sin(dt*state[8])*cos(dt*state[7])*cos(state[0])*cos(state[1]))/sqrt(1 - pow(sin(dt*state[7])*cos(state[0])*cos(state[1]) - sin(dt*state[8])*sin(state[0])*cos(dt*state[7])*cos(state[1]) + sin(state[1])*cos(dt*state[7])*cos(dt*state[8]), 2));
   out_1059548545756754699[19] = (-sin(dt*state[7])*sin(state[1])*cos(state[0]) + sin(dt*state[8])*sin(state[0])*sin(state[1])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1]))/sqrt(1 - pow(sin(dt*state[7])*cos(state[0])*cos(state[1]) - sin(dt*state[8])*sin(state[0])*cos(dt*state[7])*cos(state[1]) + sin(state[1])*cos(dt*state[7])*cos(dt*state[8]), 2));
   out_1059548545756754699[20] = 0;
   out_1059548545756754699[21] = 0;
   out_1059548545756754699[22] = 0;
   out_1059548545756754699[23] = 0;
   out_1059548545756754699[24] = 0;
   out_1059548545756754699[25] = (dt*sin(dt*state[7])*sin(dt*state[8])*sin(state[0])*cos(state[1]) - dt*sin(dt*state[7])*sin(state[1])*cos(dt*state[8]) + dt*cos(dt*state[7])*cos(state[0])*cos(state[1]))/sqrt(1 - pow(sin(dt*state[7])*cos(state[0])*cos(state[1]) - sin(dt*state[8])*sin(state[0])*cos(dt*state[7])*cos(state[1]) + sin(state[1])*cos(dt*state[7])*cos(dt*state[8]), 2));
   out_1059548545756754699[26] = (-dt*sin(dt*state[8])*sin(state[1])*cos(dt*state[7]) - dt*sin(state[0])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]))/sqrt(1 - pow(sin(dt*state[7])*cos(state[0])*cos(state[1]) - sin(dt*state[8])*sin(state[0])*cos(dt*state[7])*cos(state[1]) + sin(state[1])*cos(dt*state[7])*cos(dt*state[8]), 2));
   out_1059548545756754699[27] = 0;
   out_1059548545756754699[28] = 0;
   out_1059548545756754699[29] = 0;
   out_1059548545756754699[30] = 0;
   out_1059548545756754699[31] = 0;
   out_1059548545756754699[32] = 0;
   out_1059548545756754699[33] = 0;
   out_1059548545756754699[34] = 0;
   out_1059548545756754699[35] = 0;
   out_1059548545756754699[36] = ((sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[7]))*((-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) - (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) - sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2)) + ((-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[7]))*(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2));
   out_1059548545756754699[37] = (-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]))*(-sin(dt*state[7])*sin(state[2])*cos(state[0])*cos(state[1]) + sin(dt*state[8])*sin(state[0])*sin(state[2])*cos(dt*state[7])*cos(state[1]) - sin(state[1])*sin(state[2])*cos(dt*state[7])*cos(dt*state[8]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2)) + ((-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) - (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) - sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]))*(-sin(dt*state[7])*cos(state[0])*cos(state[1])*cos(state[2]) + sin(dt*state[8])*sin(state[0])*cos(dt*state[7])*cos(state[1])*cos(state[2]) - sin(state[1])*cos(dt*state[7])*cos(dt*state[8])*cos(state[2]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2));
   out_1059548545756754699[38] = ((-sin(state[0])*sin(state[2]) - sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]))*(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2)) + ((-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (-sin(state[0])*sin(state[1])*sin(state[2]) - cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) - sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]))*((-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) - (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) - sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2));
   out_1059548545756754699[39] = 0;
   out_1059548545756754699[40] = 0;
   out_1059548545756754699[41] = 0;
   out_1059548545756754699[42] = 0;
   out_1059548545756754699[43] = (-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]))*(dt*(sin(state[0])*cos(state[2]) - sin(state[1])*sin(state[2])*cos(state[0]))*cos(dt*state[7]) - dt*(sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[7])*sin(dt*state[8]) - dt*sin(dt*state[7])*sin(state[2])*cos(dt*state[8])*cos(state[1]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2)) + ((-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) - (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) - sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]))*(dt*(-sin(state[0])*sin(state[2]) - sin(state[1])*cos(state[0])*cos(state[2]))*cos(dt*state[7]) - dt*(sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[7])*sin(dt*state[8]) - dt*sin(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2));
   out_1059548545756754699[44] = (dt*(sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*cos(dt*state[7])*cos(dt*state[8]) - dt*sin(dt*state[8])*sin(state[2])*cos(dt*state[7])*cos(state[1]))*(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2)) + (dt*(sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*cos(dt*state[7])*cos(dt*state[8]) - dt*sin(dt*state[8])*cos(dt*state[7])*cos(state[1])*cos(state[2]))*((-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) - (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) - sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]))/(pow(-(sin(state[0])*sin(state[2]) + sin(state[1])*cos(state[0])*cos(state[2]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*cos(state[2]) - sin(state[2])*cos(state[0]))*sin(dt*state[8])*cos(dt*state[7]) + cos(dt*state[7])*cos(dt*state[8])*cos(state[1])*cos(state[2]), 2) + pow(-(-sin(state[0])*cos(state[2]) + sin(state[1])*sin(state[2])*cos(state[0]))*sin(dt*state[7]) + (sin(state[0])*sin(state[1])*sin(state[2]) + cos(state[0])*cos(state[2]))*sin(dt*state[8])*cos(dt*state[7]) + sin(state[2])*cos(dt*state[7])*cos(dt*state[8])*cos(state[1]), 2));
   out_1059548545756754699[45] = 0;
   out_1059548545756754699[46] = 0;
   out_1059548545756754699[47] = 0;
   out_1059548545756754699[48] = 0;
   out_1059548545756754699[49] = 0;
   out_1059548545756754699[50] = 0;
   out_1059548545756754699[51] = 0;
   out_1059548545756754699[52] = 0;
   out_1059548545756754699[53] = 0;
   out_1059548545756754699[54] = 0;
   out_1059548545756754699[55] = 0;
   out_1059548545756754699[56] = 0;
   out_1059548545756754699[57] = 1;
   out_1059548545756754699[58] = 0;
   out_1059548545756754699[59] = 0;
   out_1059548545756754699[60] = 0;
   out_1059548545756754699[61] = 0;
   out_1059548545756754699[62] = 0;
   out_1059548545756754699[63] = 0;
   out_1059548545756754699[64] = 0;
   out_1059548545756754699[65] = 0;
   out_1059548545756754699[66] = dt;
   out_1059548545756754699[67] = 0;
   out_1059548545756754699[68] = 0;
   out_1059548545756754699[69] = 0;
   out_1059548545756754699[70] = 0;
   out_1059548545756754699[71] = 0;
   out_1059548545756754699[72] = 0;
   out_1059548545756754699[73] = 0;
   out_1059548545756754699[74] = 0;
   out_1059548545756754699[75] = 0;
   out_1059548545756754699[76] = 1;
   out_1059548545756754699[77] = 0;
   out_1059548545756754699[78] = 0;
   out_1059548545756754699[79] = 0;
   out_1059548545756754699[80] = 0;
   out_1059548545756754699[81] = 0;
   out_1059548545756754699[82] = 0;
   out_1059548545756754699[83] = 0;
   out_1059548545756754699[84] = 0;
   out_1059548545756754699[85] = dt;
   out_1059548545756754699[86] = 0;
   out_1059548545756754699[87] = 0;
   out_1059548545756754699[88] = 0;
   out_1059548545756754699[89] = 0;
   out_1059548545756754699[90] = 0;
   out_1059548545756754699[91] = 0;
   out_1059548545756754699[92] = 0;
   out_1059548545756754699[93] = 0;
   out_1059548545756754699[94] = 0;
   out_1059548545756754699[95] = 1;
   out_1059548545756754699[96] = 0;
   out_1059548545756754699[97] = 0;
   out_1059548545756754699[98] = 0;
   out_1059548545756754699[99] = 0;
   out_1059548545756754699[100] = 0;
   out_1059548545756754699[101] = 0;
   out_1059548545756754699[102] = 0;
   out_1059548545756754699[103] = 0;
   out_1059548545756754699[104] = dt;
   out_1059548545756754699[105] = 0;
   out_1059548545756754699[106] = 0;
   out_1059548545756754699[107] = 0;
   out_1059548545756754699[108] = 0;
   out_1059548545756754699[109] = 0;
   out_1059548545756754699[110] = 0;
   out_1059548545756754699[111] = 0;
   out_1059548545756754699[112] = 0;
   out_1059548545756754699[113] = 0;
   out_1059548545756754699[114] = 1;
   out_1059548545756754699[115] = 0;
   out_1059548545756754699[116] = 0;
   out_1059548545756754699[117] = 0;
   out_1059548545756754699[118] = 0;
   out_1059548545756754699[119] = 0;
   out_1059548545756754699[120] = 0;
   out_1059548545756754699[121] = 0;
   out_1059548545756754699[122] = 0;
   out_1059548545756754699[123] = 0;
   out_1059548545756754699[124] = 0;
   out_1059548545756754699[125] = 0;
   out_1059548545756754699[126] = 0;
   out_1059548545756754699[127] = 0;
   out_1059548545756754699[128] = 0;
   out_1059548545756754699[129] = 0;
   out_1059548545756754699[130] = 0;
   out_1059548545756754699[131] = 0;
   out_1059548545756754699[132] = 0;
   out_1059548545756754699[133] = 1;
   out_1059548545756754699[134] = 0;
   out_1059548545756754699[135] = 0;
   out_1059548545756754699[136] = 0;
   out_1059548545756754699[137] = 0;
   out_1059548545756754699[138] = 0;
   out_1059548545756754699[139] = 0;
   out_1059548545756754699[140] = 0;
   out_1059548545756754699[141] = 0;
   out_1059548545756754699[142] = 0;
   out_1059548545756754699[143] = 0;
   out_1059548545756754699[144] = 0;
   out_1059548545756754699[145] = 0;
   out_1059548545756754699[146] = 0;
   out_1059548545756754699[147] = 0;
   out_1059548545756754699[148] = 0;
   out_1059548545756754699[149] = 0;
   out_1059548545756754699[150] = 0;
   out_1059548545756754699[151] = 0;
   out_1059548545756754699[152] = 1;
   out_1059548545756754699[153] = 0;
   out_1059548545756754699[154] = 0;
   out_1059548545756754699[155] = 0;
   out_1059548545756754699[156] = 0;
   out_1059548545756754699[157] = 0;
   out_1059548545756754699[158] = 0;
   out_1059548545756754699[159] = 0;
   out_1059548545756754699[160] = 0;
   out_1059548545756754699[161] = 0;
   out_1059548545756754699[162] = 0;
   out_1059548545756754699[163] = 0;
   out_1059548545756754699[164] = 0;
   out_1059548545756754699[165] = 0;
   out_1059548545756754699[166] = 0;
   out_1059548545756754699[167] = 0;
   out_1059548545756754699[168] = 0;
   out_1059548545756754699[169] = 0;
   out_1059548545756754699[170] = 0;
   out_1059548545756754699[171] = 1;
   out_1059548545756754699[172] = 0;
   out_1059548545756754699[173] = 0;
   out_1059548545756754699[174] = 0;
   out_1059548545756754699[175] = 0;
   out_1059548545756754699[176] = 0;
   out_1059548545756754699[177] = 0;
   out_1059548545756754699[178] = 0;
   out_1059548545756754699[179] = 0;
   out_1059548545756754699[180] = 0;
   out_1059548545756754699[181] = 0;
   out_1059548545756754699[182] = 0;
   out_1059548545756754699[183] = 0;
   out_1059548545756754699[184] = 0;
   out_1059548545756754699[185] = 0;
   out_1059548545756754699[186] = 0;
   out_1059548545756754699[187] = 0;
   out_1059548545756754699[188] = 0;
   out_1059548545756754699[189] = 0;
   out_1059548545756754699[190] = 1;
   out_1059548545756754699[191] = 0;
   out_1059548545756754699[192] = 0;
   out_1059548545756754699[193] = 0;
   out_1059548545756754699[194] = 0;
   out_1059548545756754699[195] = 0;
   out_1059548545756754699[196] = 0;
   out_1059548545756754699[197] = 0;
   out_1059548545756754699[198] = 0;
   out_1059548545756754699[199] = 0;
   out_1059548545756754699[200] = 0;
   out_1059548545756754699[201] = 0;
   out_1059548545756754699[202] = 0;
   out_1059548545756754699[203] = 0;
   out_1059548545756754699[204] = 0;
   out_1059548545756754699[205] = 0;
   out_1059548545756754699[206] = 0;
   out_1059548545756754699[207] = 0;
   out_1059548545756754699[208] = 0;
   out_1059548545756754699[209] = 1;
   out_1059548545756754699[210] = 0;
   out_1059548545756754699[211] = 0;
   out_1059548545756754699[212] = 0;
   out_1059548545756754699[213] = 0;
   out_1059548545756754699[214] = 0;
   out_1059548545756754699[215] = 0;
   out_1059548545756754699[216] = 0;
   out_1059548545756754699[217] = 0;
   out_1059548545756754699[218] = 0;
   out_1059548545756754699[219] = 0;
   out_1059548545756754699[220] = 0;
   out_1059548545756754699[221] = 0;
   out_1059548545756754699[222] = 0;
   out_1059548545756754699[223] = 0;
   out_1059548545756754699[224] = 0;
   out_1059548545756754699[225] = 0;
   out_1059548545756754699[226] = 0;
   out_1059548545756754699[227] = 0;
   out_1059548545756754699[228] = 1;
   out_1059548545756754699[229] = 0;
   out_1059548545756754699[230] = 0;
   out_1059548545756754699[231] = 0;
   out_1059548545756754699[232] = 0;
   out_1059548545756754699[233] = 0;
   out_1059548545756754699[234] = 0;
   out_1059548545756754699[235] = 0;
   out_1059548545756754699[236] = 0;
   out_1059548545756754699[237] = 0;
   out_1059548545756754699[238] = 0;
   out_1059548545756754699[239] = 0;
   out_1059548545756754699[240] = 0;
   out_1059548545756754699[241] = 0;
   out_1059548545756754699[242] = 0;
   out_1059548545756754699[243] = 0;
   out_1059548545756754699[244] = 0;
   out_1059548545756754699[245] = 0;
   out_1059548545756754699[246] = 0;
   out_1059548545756754699[247] = 1;
   out_1059548545756754699[248] = 0;
   out_1059548545756754699[249] = 0;
   out_1059548545756754699[250] = 0;
   out_1059548545756754699[251] = 0;
   out_1059548545756754699[252] = 0;
   out_1059548545756754699[253] = 0;
   out_1059548545756754699[254] = 0;
   out_1059548545756754699[255] = 0;
   out_1059548545756754699[256] = 0;
   out_1059548545756754699[257] = 0;
   out_1059548545756754699[258] = 0;
   out_1059548545756754699[259] = 0;
   out_1059548545756754699[260] = 0;
   out_1059548545756754699[261] = 0;
   out_1059548545756754699[262] = 0;
   out_1059548545756754699[263] = 0;
   out_1059548545756754699[264] = 0;
   out_1059548545756754699[265] = 0;
   out_1059548545756754699[266] = 1;
   out_1059548545756754699[267] = 0;
   out_1059548545756754699[268] = 0;
   out_1059548545756754699[269] = 0;
   out_1059548545756754699[270] = 0;
   out_1059548545756754699[271] = 0;
   out_1059548545756754699[272] = 0;
   out_1059548545756754699[273] = 0;
   out_1059548545756754699[274] = 0;
   out_1059548545756754699[275] = 0;
   out_1059548545756754699[276] = 0;
   out_1059548545756754699[277] = 0;
   out_1059548545756754699[278] = 0;
   out_1059548545756754699[279] = 0;
   out_1059548545756754699[280] = 0;
   out_1059548545756754699[281] = 0;
   out_1059548545756754699[282] = 0;
   out_1059548545756754699[283] = 0;
   out_1059548545756754699[284] = 0;
   out_1059548545756754699[285] = 1;
   out_1059548545756754699[286] = 0;
   out_1059548545756754699[287] = 0;
   out_1059548545756754699[288] = 0;
   out_1059548545756754699[289] = 0;
   out_1059548545756754699[290] = 0;
   out_1059548545756754699[291] = 0;
   out_1059548545756754699[292] = 0;
   out_1059548545756754699[293] = 0;
   out_1059548545756754699[294] = 0;
   out_1059548545756754699[295] = 0;
   out_1059548545756754699[296] = 0;
   out_1059548545756754699[297] = 0;
   out_1059548545756754699[298] = 0;
   out_1059548545756754699[299] = 0;
   out_1059548545756754699[300] = 0;
   out_1059548545756754699[301] = 0;
   out_1059548545756754699[302] = 0;
   out_1059548545756754699[303] = 0;
   out_1059548545756754699[304] = 1;
   out_1059548545756754699[305] = 0;
   out_1059548545756754699[306] = 0;
   out_1059548545756754699[307] = 0;
   out_1059548545756754699[308] = 0;
   out_1059548545756754699[309] = 0;
   out_1059548545756754699[310] = 0;
   out_1059548545756754699[311] = 0;
   out_1059548545756754699[312] = 0;
   out_1059548545756754699[313] = 0;
   out_1059548545756754699[314] = 0;
   out_1059548545756754699[315] = 0;
   out_1059548545756754699[316] = 0;
   out_1059548545756754699[317] = 0;
   out_1059548545756754699[318] = 0;
   out_1059548545756754699[319] = 0;
   out_1059548545756754699[320] = 0;
   out_1059548545756754699[321] = 0;
   out_1059548545756754699[322] = 0;
   out_1059548545756754699[323] = 1;
}
void h_4(double *state, double *unused, double *out_1272267543225290969) {
   out_1272267543225290969[0] = state[6] + state[9];
   out_1272267543225290969[1] = state[7] + state[10];
   out_1272267543225290969[2] = state[8] + state[11];
}
void H_4(double *state, double *unused, double *out_4532055235799853873) {
   out_4532055235799853873[0] = 0;
   out_4532055235799853873[1] = 0;
   out_4532055235799853873[2] = 0;
   out_4532055235799853873[3] = 0;
   out_4532055235799853873[4] = 0;
   out_4532055235799853873[5] = 0;
   out_4532055235799853873[6] = 1;
   out_4532055235799853873[7] = 0;
   out_4532055235799853873[8] = 0;
   out_4532055235799853873[9] = 1;
   out_4532055235799853873[10] = 0;
   out_4532055235799853873[11] = 0;
   out_4532055235799853873[12] = 0;
   out_4532055235799853873[13] = 0;
   out_4532055235799853873[14] = 0;
   out_4532055235799853873[15] = 0;
   out_4532055235799853873[16] = 0;
   out_4532055235799853873[17] = 0;
   out_4532055235799853873[18] = 0;
   out_4532055235799853873[19] = 0;
   out_4532055235799853873[20] = 0;
   out_4532055235799853873[21] = 0;
   out_4532055235799853873[22] = 0;
   out_4532055235799853873[23] = 0;
   out_4532055235799853873[24] = 0;
   out_4532055235799853873[25] = 1;
   out_4532055235799853873[26] = 0;
   out_4532055235799853873[27] = 0;
   out_4532055235799853873[28] = 1;
   out_4532055235799853873[29] = 0;
   out_4532055235799853873[30] = 0;
   out_4532055235799853873[31] = 0;
   out_4532055235799853873[32] = 0;
   out_4532055235799853873[33] = 0;
   out_4532055235799853873[34] = 0;
   out_4532055235799853873[35] = 0;
   out_4532055235799853873[36] = 0;
   out_4532055235799853873[37] = 0;
   out_4532055235799853873[38] = 0;
   out_4532055235799853873[39] = 0;
   out_4532055235799853873[40] = 0;
   out_4532055235799853873[41] = 0;
   out_4532055235799853873[42] = 0;
   out_4532055235799853873[43] = 0;
   out_4532055235799853873[44] = 1;
   out_4532055235799853873[45] = 0;
   out_4532055235799853873[46] = 0;
   out_4532055235799853873[47] = 1;
   out_4532055235799853873[48] = 0;
   out_4532055235799853873[49] = 0;
   out_4532055235799853873[50] = 0;
   out_4532055235799853873[51] = 0;
   out_4532055235799853873[52] = 0;
   out_4532055235799853873[53] = 0;
}
void h_10(double *state, double *unused, double *out_4823025017146951895) {
   out_4823025017146951895[0] = 9.8100000000000005*sin(state[1]) - state[4]*state[8] + state[5]*state[7] + state[12] + state[15];
   out_4823025017146951895[1] = -9.8100000000000005*sin(state[0])*cos(state[1]) + state[3]*state[8] - state[5]*state[6] + state[13] + state[16];
   out_4823025017146951895[2] = -9.8100000000000005*cos(state[0])*cos(state[1]) - state[3]*state[7] + state[4]*state[6] + state[14] + state[17];
}
void H_10(double *state, double *unused, double *out_7328110893683068042) {
   out_7328110893683068042[0] = 0;
   out_7328110893683068042[1] = 9.8100000000000005*cos(state[1]);
   out_7328110893683068042[2] = 0;
   out_7328110893683068042[3] = 0;
   out_7328110893683068042[4] = -state[8];
   out_7328110893683068042[5] = state[7];
   out_7328110893683068042[6] = 0;
   out_7328110893683068042[7] = state[5];
   out_7328110893683068042[8] = -state[4];
   out_7328110893683068042[9] = 0;
   out_7328110893683068042[10] = 0;
   out_7328110893683068042[11] = 0;
   out_7328110893683068042[12] = 1;
   out_7328110893683068042[13] = 0;
   out_7328110893683068042[14] = 0;
   out_7328110893683068042[15] = 1;
   out_7328110893683068042[16] = 0;
   out_7328110893683068042[17] = 0;
   out_7328110893683068042[18] = -9.8100000000000005*cos(state[0])*cos(state[1]);
   out_7328110893683068042[19] = 9.8100000000000005*sin(state[0])*sin(state[1]);
   out_7328110893683068042[20] = 0;
   out_7328110893683068042[21] = state[8];
   out_7328110893683068042[22] = 0;
   out_7328110893683068042[23] = -state[6];
   out_7328110893683068042[24] = -state[5];
   out_7328110893683068042[25] = 0;
   out_7328110893683068042[26] = state[3];
   out_7328110893683068042[27] = 0;
   out_7328110893683068042[28] = 0;
   out_7328110893683068042[29] = 0;
   out_7328110893683068042[30] = 0;
   out_7328110893683068042[31] = 1;
   out_7328110893683068042[32] = 0;
   out_7328110893683068042[33] = 0;
   out_7328110893683068042[34] = 1;
   out_7328110893683068042[35] = 0;
   out_7328110893683068042[36] = 9.8100000000000005*sin(state[0])*cos(state[1]);
   out_7328110893683068042[37] = 9.8100000000000005*sin(state[1])*cos(state[0]);
   out_7328110893683068042[38] = 0;
   out_7328110893683068042[39] = -state[7];
   out_7328110893683068042[40] = state[6];
   out_7328110893683068042[41] = 0;
   out_7328110893683068042[42] = state[4];
   out_7328110893683068042[43] = -state[3];
   out_7328110893683068042[44] = 0;
   out_7328110893683068042[45] = 0;
   out_7328110893683068042[46] = 0;
   out_7328110893683068042[47] = 0;
   out_7328110893683068042[48] = 0;
   out_7328110893683068042[49] = 0;
   out_7328110893683068042[50] = 1;
   out_7328110893683068042[51] = 0;
   out_7328110893683068042[52] = 0;
   out_7328110893683068042[53] = 1;
}
void h_13(double *state, double *unused, double *out_3553628753437543153) {
   out_3553628753437543153[0] = state[3];
   out_3553628753437543153[1] = state[4];
   out_3553628753437543153[2] = state[5];
}
void H_13(double *state, double *unused, double *out_7744329061132186674) {
   out_7744329061132186674[0] = 0;
   out_7744329061132186674[1] = 0;
   out_7744329061132186674[2] = 0;
   out_7744329061132186674[3] = 1;
   out_7744329061132186674[4] = 0;
   out_7744329061132186674[5] = 0;
   out_7744329061132186674[6] = 0;
   out_7744329061132186674[7] = 0;
   out_7744329061132186674[8] = 0;
   out_7744329061132186674[9] = 0;
   out_7744329061132186674[10] = 0;
   out_7744329061132186674[11] = 0;
   out_7744329061132186674[12] = 0;
   out_7744329061132186674[13] = 0;
   out_7744329061132186674[14] = 0;
   out_7744329061132186674[15] = 0;
   out_7744329061132186674[16] = 0;
   out_7744329061132186674[17] = 0;
   out_7744329061132186674[18] = 0;
   out_7744329061132186674[19] = 0;
   out_7744329061132186674[20] = 0;
   out_7744329061132186674[21] = 0;
   out_7744329061132186674[22] = 1;
   out_7744329061132186674[23] = 0;
   out_7744329061132186674[24] = 0;
   out_7744329061132186674[25] = 0;
   out_7744329061132186674[26] = 0;
   out_7744329061132186674[27] = 0;
   out_7744329061132186674[28] = 0;
   out_7744329061132186674[29] = 0;
   out_7744329061132186674[30] = 0;
   out_7744329061132186674[31] = 0;
   out_7744329061132186674[32] = 0;
   out_7744329061132186674[33] = 0;
   out_7744329061132186674[34] = 0;
   out_7744329061132186674[35] = 0;
   out_7744329061132186674[36] = 0;
   out_7744329061132186674[37] = 0;
   out_7744329061132186674[38] = 0;
   out_7744329061132186674[39] = 0;
   out_7744329061132186674[40] = 0;
   out_7744329061132186674[41] = 1;
   out_7744329061132186674[42] = 0;
   out_7744329061132186674[43] = 0;
   out_7744329061132186674[44] = 0;
   out_7744329061132186674[45] = 0;
   out_7744329061132186674[46] = 0;
   out_7744329061132186674[47] = 0;
   out_7744329061132186674[48] = 0;
   out_7744329061132186674[49] = 0;
   out_7744329061132186674[50] = 0;
   out_7744329061132186674[51] = 0;
   out_7744329061132186674[52] = 0;
   out_7744329061132186674[53] = 0;
}
void h_14(double *state, double *unused, double *out_2142424780536854880) {
   out_2142424780536854880[0] = state[6];
   out_2142424780536854880[1] = state[7];
   out_2142424780536854880[2] = state[8];
}
void H_14(double *state, double *unused, double *out_8495296092139338402) {
   out_8495296092139338402[0] = 0;
   out_8495296092139338402[1] = 0;
   out_8495296092139338402[2] = 0;
   out_8495296092139338402[3] = 0;
   out_8495296092139338402[4] = 0;
   out_8495296092139338402[5] = 0;
   out_8495296092139338402[6] = 1;
   out_8495296092139338402[7] = 0;
   out_8495296092139338402[8] = 0;
   out_8495296092139338402[9] = 0;
   out_8495296092139338402[10] = 0;
   out_8495296092139338402[11] = 0;
   out_8495296092139338402[12] = 0;
   out_8495296092139338402[13] = 0;
   out_8495296092139338402[14] = 0;
   out_8495296092139338402[15] = 0;
   out_8495296092139338402[16] = 0;
   out_8495296092139338402[17] = 0;
   out_8495296092139338402[18] = 0;
   out_8495296092139338402[19] = 0;
   out_8495296092139338402[20] = 0;
   out_8495296092139338402[21] = 0;
   out_8495296092139338402[22] = 0;
   out_8495296092139338402[23] = 0;
   out_8495296092139338402[24] = 0;
   out_8495296092139338402[25] = 1;
   out_8495296092139338402[26] = 0;
   out_8495296092139338402[27] = 0;
   out_8495296092139338402[28] = 0;
   out_8495296092139338402[29] = 0;
   out_8495296092139338402[30] = 0;
   out_8495296092139338402[31] = 0;
   out_8495296092139338402[32] = 0;
   out_8495296092139338402[33] = 0;
   out_8495296092139338402[34] = 0;
   out_8495296092139338402[35] = 0;
   out_8495296092139338402[36] = 0;
   out_8495296092139338402[37] = 0;
   out_8495296092139338402[38] = 0;
   out_8495296092139338402[39] = 0;
   out_8495296092139338402[40] = 0;
   out_8495296092139338402[41] = 0;
   out_8495296092139338402[42] = 0;
   out_8495296092139338402[43] = 0;
   out_8495296092139338402[44] = 1;
   out_8495296092139338402[45] = 0;
   out_8495296092139338402[46] = 0;
   out_8495296092139338402[47] = 0;
   out_8495296092139338402[48] = 0;
   out_8495296092139338402[49] = 0;
   out_8495296092139338402[50] = 0;
   out_8495296092139338402[51] = 0;
   out_8495296092139338402[52] = 0;
   out_8495296092139338402[53] = 0;
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

void pose_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<3, 3, 0>(in_x, in_P, h_4, H_4, NULL, in_z, in_R, in_ea, MAHA_THRESH_4);
}
void pose_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<3, 3, 0>(in_x, in_P, h_10, H_10, NULL, in_z, in_R, in_ea, MAHA_THRESH_10);
}
void pose_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<3, 3, 0>(in_x, in_P, h_13, H_13, NULL, in_z, in_R, in_ea, MAHA_THRESH_13);
}
void pose_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<3, 3, 0>(in_x, in_P, h_14, H_14, NULL, in_z, in_R, in_ea, MAHA_THRESH_14);
}
void pose_err_fun(double *nom_x, double *delta_x, double *out_4750261472376034529) {
  err_fun(nom_x, delta_x, out_4750261472376034529);
}
void pose_inv_err_fun(double *nom_x, double *true_x, double *out_8277444623959904026) {
  inv_err_fun(nom_x, true_x, out_8277444623959904026);
}
void pose_H_mod_fun(double *state, double *out_7806277314590803126) {
  H_mod_fun(state, out_7806277314590803126);
}
void pose_f_fun(double *state, double dt, double *out_4548647388261971274) {
  f_fun(state,  dt, out_4548647388261971274);
}
void pose_F_fun(double *state, double dt, double *out_1059548545756754699) {
  F_fun(state,  dt, out_1059548545756754699);
}
void pose_h_4(double *state, double *unused, double *out_1272267543225290969) {
  h_4(state, unused, out_1272267543225290969);
}
void pose_H_4(double *state, double *unused, double *out_4532055235799853873) {
  H_4(state, unused, out_4532055235799853873);
}
void pose_h_10(double *state, double *unused, double *out_4823025017146951895) {
  h_10(state, unused, out_4823025017146951895);
}
void pose_H_10(double *state, double *unused, double *out_7328110893683068042) {
  H_10(state, unused, out_7328110893683068042);
}
void pose_h_13(double *state, double *unused, double *out_3553628753437543153) {
  h_13(state, unused, out_3553628753437543153);
}
void pose_H_13(double *state, double *unused, double *out_7744329061132186674) {
  H_13(state, unused, out_7744329061132186674);
}
void pose_h_14(double *state, double *unused, double *out_2142424780536854880) {
  h_14(state, unused, out_2142424780536854880);
}
void pose_H_14(double *state, double *unused, double *out_8495296092139338402) {
  H_14(state, unused, out_8495296092139338402);
}
void pose_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
}

const EKF pose = {
  .name = "pose",
  .kinds = { 4, 10, 13, 14 },
  .feature_kinds = {  },
  .f_fun = pose_f_fun,
  .F_fun = pose_F_fun,
  .err_fun = pose_err_fun,
  .inv_err_fun = pose_inv_err_fun,
  .H_mod_fun = pose_H_mod_fun,
  .predict = pose_predict,
  .hs = {
    { 4, pose_h_4 },
    { 10, pose_h_10 },
    { 13, pose_h_13 },
    { 14, pose_h_14 },
  },
  .Hs = {
    { 4, pose_H_4 },
    { 10, pose_H_10 },
    { 13, pose_H_13 },
    { 14, pose_H_14 },
  },
  .updates = {
    { 4, pose_update_4 },
    { 10, pose_update_10 },
    { 13, pose_update_13 },
    { 14, pose_update_14 },
  },
  .Hes = {
  },
  .sets = {
  },
  .extra_routines = {
  },
};

ekf_lib_init(pose)
