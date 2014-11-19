/* Include files */

#include <stddef.h>
#include "blas.h"
#include "hw8_two_tanks_sfun.h"
#include "c3_hw8_two_tanks.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "hw8_two_tanks_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)
#define c3_IN_NO_ACTIVE_CHILD          ((uint8_T)0U)
#define c3_IN_Balancing                ((uint8_T)1U)
#define c3_IN_From_tank_1              ((uint8_T)2U)
#define c3_IN_From_tank_2              ((uint8_T)3U)
#define c3_IN_Separated                ((uint8_T)4U)
#define c3_const_h                     (0.3)
#define c3_const_A                     (0.0154)
#define c3_const_rho                   (1000.0)
#define c3_const_g                     (9.81)
#define c3_const_Kin                   (0.06)
#define c3_const_Ka                    (0.001)
#define c3_const_Kout                  (0.001)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c3_debug_family_names[2] = { "nargin", "nargout" };

static const char * c3_b_debug_family_names[3] = { "nargin", "nargout", "q" };

static const char * c3_c_debug_family_names[6] = { "delta", "nargin", "nargout",
  "ha", "hb", "qa" };

static const char * c3_d_debug_family_names[3] = { "nargin", "nargout", "q" };

static const char * c3_e_debug_family_names[2] = { "nargin", "nargout" };

static const char * c3_f_debug_family_names[2] = { "nargin", "nargout" };

static const char * c3_g_debug_family_names[6] = { "delta", "nargin", "nargout",
  "ha", "hb", "q" };

static const char * c3_h_debug_family_names[2] = { "nargin", "nargout" };

static const char * c3_i_debug_family_names[2] = { "nargin", "nargout" };

static const char * c3_j_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static const char * c3_k_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static const char * c3_l_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static const char * c3_m_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static const char * c3_n_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static const char * c3_o_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static const char * c3_p_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static const char * c3_q_debug_family_names[3] = { "nargin", "nargout",
  "sf_internal_predicateOutput" };

static boolean_T c3_dataWrittenToVector[12];

/* Function Declarations */
static void initialize_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void initialize_params_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct *
  chartInstance);
static void enable_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void disable_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void c3_update_debugger_state_c3_hw8_two_tanks
  (SFc3_hw8_two_tanksInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c3_hw8_two_tanks
  (SFc3_hw8_two_tanksInstanceStruct *chartInstance);
static void set_sim_state_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_st);
static void c3_set_sim_state_side_effects_c3_hw8_two_tanks
  (SFc3_hw8_two_tanksInstanceStruct *chartInstance);
static void finalize_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void sf_gateway_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void zeroCrossings_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void derivatives_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void outputs_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void initSimStructsc3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void c3_eml_ini_fcn_to_be_inlined_23(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void c3_eml_term_fcn_to_be_inlined_23(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static real_T c3_sign(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                      c3_x);
static void c3_isBuiltInNumeric(SFc3_hw8_two_tanksInstanceStruct *chartInstance);
static real_T c3_eml_scalar_sign(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T c3_x);
static real_T c3_abs(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                     c3_x);
static real_T c3_eml_scalar_abs(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T c3_x);
static real_T c3_sqrt(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                      c3_x);
static void c3_eml_error(SFc3_hw8_two_tanksInstanceStruct *chartInstance);
static real_T c3_eml_scalar_sqrt(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T c3_x);
static real_T c3_mrdivide(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T c3_b_A, real_T c3_B);
static void c3_assert(SFc3_hw8_two_tanksInstanceStruct *chartInstance);
static real_T c3_rdivide(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
  c3_x, real_T c3_y);
static void c3_eml_scalexp_compatible(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static real_T c3_eml_div(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
  c3_x, real_T c3_y);
static real_T c3_div(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                     c3_x, real_T c3_y);
static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber);
static const mxArray *c3_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const char_T c3_u[30]);
static const mxArray *c3_b_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const char_T c3_u[4]);
static const mxArray *c3_c_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const real_T c3_u);
static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData);
static real_T c3_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_nargout, const char_T *c3_identifier);
static real_T c3_b_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_d_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const boolean_T c3_u);
static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static boolean_T c3_c_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_sf_internal_predicateOutput, const char_T
  *c3_identifier);
static boolean_T c3_d_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static void c3_info_helper(const mxArray **c3_info);
static const mxArray *c3_e_emlrt_marshallOut(const char * c3_u);
static const mxArray *c3_f_emlrt_marshallOut(const uint32_T c3_u);
static const mxArray *c3_g_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const int32_T c3_u);
static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static int32_T c3_e_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_b_sfEvent, const char_T *c3_identifier);
static int32_T c3_f_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_h_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const uint8_T c3_u);
static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData);
static uint8_T c3_g_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_b_tp_From_tank_1, const char_T
  *c3_identifier);
static uint8_T c3_h_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData);
static const mxArray *c3_i_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);
static void c3_i_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u);
static const mxArray *c3_j_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_b_setSimStateSideEffectsInfo, const char_T
  *c3_identifier);
static const mxArray *c3_k_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId);
static void c3_updateDataWrittenToVector(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, uint32_T c3_vectorIndex);
static void c3_errorIfDataNotWrittenToFcn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, uint32_T c3_vectorIndex, uint32_T c3_dataNumber);
static void c3_b_sign(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                      *c3_x);
static void c3_b_eml_scalar_sign(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T *c3_x);
static void c3_b_sqrt(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                      *c3_x);
static void c3_b_eml_scalar_sqrt(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T *c3_x);
static void init_dsm_address_info(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  chartInstance->c3_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c3_doSetSimStateSideEffects = 0U;
  chartInstance->c3_setSimStateSideEffectsInfo = NULL;
  chartInstance->c3_tp_Balancing = 0U;
  chartInstance->c3_tp_From_tank_1 = 0U;
  chartInstance->c3_tp_From_tank_2 = 0U;
  chartInstance->c3_tp_Separated = 0U;
  chartInstance->c3_is_active_c3_hw8_two_tanks = 0U;
  chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_NO_ACTIVE_CHILD;
  chartInstance->c3_h = 0.3;
  chartInstance->c3_A = 0.0154;
  chartInstance->c3_rho = 1000.0;
  chartInstance->c3_g = 9.81;
  chartInstance->c3_Kin = 0.06;
  chartInstance->c3_Ka = 0.001;
  chartInstance->c3_Kout = 0.001;
}

static void initialize_params_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct *
  chartInstance)
{
  (void)chartInstance;
}

static void enable_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c3_update_debugger_state_c3_hw8_two_tanks
  (SFc3_hw8_two_tanksInstanceStruct *chartInstance)
{
  uint32_T c3_prevAniVal;
  c3_prevAniVal = _SFD_GET_ANIMATION();
  _SFD_SET_ANIMATION(0U);
  _SFD_SET_HONOR_BREAKPOINTS(0U);
  if (chartInstance->c3_is_active_c3_hw8_two_tanks == 1U) {
    _SFD_CC_CALL(CHART_ACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
  }

  if (chartInstance->c3_is_c3_hw8_two_tanks == c3_IN_From_tank_1) {
    _SFD_CS_CALL(STATE_ACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
  } else {
    _SFD_CS_CALL(STATE_INACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
  }

  if (chartInstance->c3_is_c3_hw8_two_tanks == c3_IN_Separated) {
    _SFD_CS_CALL(STATE_ACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
  } else {
    _SFD_CS_CALL(STATE_INACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
  }

  if (chartInstance->c3_is_c3_hw8_two_tanks == c3_IN_Balancing) {
    _SFD_CS_CALL(STATE_ACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
  } else {
    _SFD_CS_CALL(STATE_INACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
  }

  if (chartInstance->c3_is_c3_hw8_two_tanks == c3_IN_From_tank_2) {
    _SFD_CS_CALL(STATE_ACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
  } else {
    _SFD_CS_CALL(STATE_INACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
  }

  _SFD_SET_ANIMATION(c3_prevAniVal);
  _SFD_SET_HONOR_BREAKPOINTS(1U);
  _SFD_ANIMATE();
}

static const mxArray *get_sim_state_c3_hw8_two_tanks
  (SFc3_hw8_two_tanksInstanceStruct *chartInstance)
{
  const mxArray *c3_st = NULL;
  c3_st = NULL;
  sf_mex_assign(&c3_st, c3_i_emlrt_marshallOut(chartInstance), false);
  return c3_st;
}

static void set_sim_state_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_st)
{
  c3_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_st));
  chartInstance->c3_doSetSimStateSideEffects = 1U;
  c3_update_debugger_state_c3_hw8_two_tanks(chartInstance);
  sf_mex_destroy(&c3_st);
}

static void c3_set_sim_state_side_effects_c3_hw8_two_tanks
  (SFc3_hw8_two_tanksInstanceStruct *chartInstance)
{
  if (chartInstance->c3_doSetSimStateSideEffects != 0) {
    if (chartInstance->c3_is_c3_hw8_two_tanks == c3_IN_Balancing) {
      chartInstance->c3_tp_Balancing = 1U;
    } else {
      chartInstance->c3_tp_Balancing = 0U;
    }

    if (chartInstance->c3_is_c3_hw8_two_tanks == c3_IN_From_tank_1) {
      chartInstance->c3_tp_From_tank_1 = 1U;
    } else {
      chartInstance->c3_tp_From_tank_1 = 0U;
    }

    if (chartInstance->c3_is_c3_hw8_two_tanks == c3_IN_From_tank_2) {
      chartInstance->c3_tp_From_tank_2 = 1U;
    } else {
      chartInstance->c3_tp_From_tank_2 = 0U;
    }

    if (chartInstance->c3_is_c3_hw8_two_tanks == c3_IN_Separated) {
      chartInstance->c3_tp_Separated = 1U;
    } else {
      chartInstance->c3_tp_Separated = 0U;
    }

    chartInstance->c3_doSetSimStateSideEffects = 0U;
  }
}

static void finalize_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  sf_mex_destroy(&chartInstance->c3_setSimStateSideEffectsInfo);
}

static void sf_gateway_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  uint32_T c3_debug_family_var_map[2];
  real_T c3_nargin = 0.0;
  real_T c3_nargout = 0.0;
  uint32_T c3_b_debug_family_var_map[3];
  real_T c3_b_nargin = 0.0;
  real_T c3_b_nargout = 1.0;
  boolean_T c3_out;
  real_T c3_c_nargin = 0.0;
  real_T c3_c_nargout = 1.0;
  boolean_T c3_b_out;
  real_T c3_d_nargin = 0.0;
  real_T c3_d_nargout = 1.0;
  boolean_T c3_c_out;
  real_T c3_e_nargin = 0.0;
  real_T c3_e_nargout = 1.0;
  boolean_T c3_d_out;
  real_T c3_f_nargin = 0.0;
  real_T c3_f_nargout = 1.0;
  boolean_T c3_e_out;
  real_T c3_g_nargin = 0.0;
  real_T c3_g_nargout = 1.0;
  boolean_T c3_f_out;
  real_T c3_h_nargin = 0.0;
  real_T c3_h_nargout = 1.0;
  boolean_T c3_g_out;
  real_T c3_i_nargin = 0.0;
  real_T c3_i_nargout = 1.0;
  boolean_T c3_h_out;
  real_T c3_j_nargin = 0.0;
  real_T c3_j_nargout = 1.0;
  boolean_T c3_i_out;
  real_T c3_k_nargin = 0.0;
  real_T c3_k_nargout = 1.0;
  boolean_T c3_j_out;
  real_T c3_l_nargin = 0.0;
  real_T c3_l_nargout = 1.0;
  boolean_T c3_k_out;
  real_T c3_m_nargin = 0.0;
  real_T c3_m_nargout = 1.0;
  boolean_T c3_l_out;
  real_T c3_n_nargin = 0.0;
  real_T c3_n_nargout = 0.0;
  real_T c3_o_nargin = 0.0;
  real_T c3_o_nargout = 1.0;
  real_T c3_q;
  real_T c3_hoistedGlobal;
  real_T c3_b_hoistedGlobal;
  real_T c3_ha;
  real_T c3_hb;
  uint32_T c3_c_debug_family_var_map[6];
  real_T c3_delta;
  real_T c3_p_nargin = 2.0;
  real_T c3_p_nargout = 1.0;
  real_T c3_b_q;
  real_T c3_d0;
  real_T c3_d1;
  real_T c3_c_hoistedGlobal;
  real_T c3_d_hoistedGlobal;
  real_T c3_b_ha;
  real_T c3_b_hb;
  real_T c3_b_delta;
  real_T c3_q_nargin = 2.0;
  real_T c3_q_nargout = 1.0;
  real_T c3_c_q;
  real_T c3_d2;
  real_T c3_d3;
  real_T c3_r_nargin = 0.0;
  real_T c3_r_nargout = 1.0;
  real_T c3_d_q;
  real_T c3_d4;
  real_T c3_d5;
  real_T c3_s_nargin = 0.0;
  real_T c3_s_nargout = 1.0;
  real_T c3_e_q;
  real_T c3_e_hoistedGlobal;
  real_T c3_f_hoistedGlobal;
  real_T c3_c_ha;
  real_T c3_c_hb;
  real_T c3_c_delta;
  real_T c3_t_nargin = 2.0;
  real_T c3_t_nargout = 1.0;
  real_T c3_f_q;
  real_T c3_d6;
  real_T c3_d7;
  real_T c3_u_nargin = 0.0;
  real_T c3_u_nargout = 1.0;
  real_T c3_g_q;
  real_T c3_d8;
  real_T c3_d9;
  real_T c3_v_nargin = 0.0;
  real_T c3_v_nargout = 0.0;
  real_T c3_w_nargin = 0.0;
  real_T c3_w_nargout = 1.0;
  real_T c3_h_q;
  real_T c3_g_hoistedGlobal;
  real_T c3_d_ha;
  real_T c3_d_hb;
  real_T c3_d_delta;
  real_T c3_x_nargin = 2.0;
  real_T c3_x_nargout = 1.0;
  real_T c3_qa;
  real_T c3_d10;
  real_T c3_d11;
  real_T c3_h_hoistedGlobal;
  real_T c3_e_ha;
  real_T c3_e_hb;
  real_T c3_e_delta;
  real_T c3_y_nargin = 2.0;
  real_T c3_y_nargout = 1.0;
  real_T c3_b_qa;
  real_T c3_d12;
  real_T c3_d13;
  real_T c3_ab_nargin = 0.0;
  real_T c3_ab_nargout = 1.0;
  real_T c3_i_q;
  real_T c3_d14;
  real_T c3_d15;
  real_T c3_bb_nargin = 0.0;
  real_T c3_bb_nargout = 1.0;
  real_T c3_j_q;
  real_T c3_i_hoistedGlobal;
  real_T c3_f_ha;
  real_T c3_f_hb;
  real_T c3_f_delta;
  real_T c3_cb_nargin = 2.0;
  real_T c3_cb_nargout = 1.0;
  real_T c3_c_qa;
  real_T c3_d16;
  real_T c3_d17;
  real_T c3_db_nargin = 0.0;
  real_T c3_db_nargout = 1.0;
  real_T c3_k_q;
  real_T c3_d18;
  real_T c3_d19;
  real_T c3_eb_nargin = 0.0;
  real_T c3_eb_nargout = 0.0;
  real_T c3_fb_nargin = 0.0;
  real_T c3_fb_nargout = 1.0;
  real_T c3_l_q;
  real_T c3_j_hoistedGlobal;
  real_T c3_g_ha;
  real_T c3_g_hb;
  real_T c3_g_delta;
  real_T c3_gb_nargin = 2.0;
  real_T c3_gb_nargout = 1.0;
  real_T c3_d_qa;
  real_T c3_d20;
  real_T c3_d21;
  real_T c3_k_hoistedGlobal;
  real_T c3_h_ha;
  real_T c3_h_hb;
  real_T c3_h_delta;
  real_T c3_hb_nargin = 2.0;
  real_T c3_hb_nargout = 1.0;
  real_T c3_e_qa;
  real_T c3_d22;
  real_T c3_d23;
  real_T c3_ib_nargin = 0.0;
  real_T c3_ib_nargout = 1.0;
  real_T c3_m_q;
  real_T c3_d24;
  real_T c3_d25;
  real_T c3_jb_nargin = 0.0;
  real_T c3_jb_nargout = 1.0;
  real_T c3_n_q;
  real_T c3_l_hoistedGlobal;
  real_T c3_i_ha;
  real_T c3_i_hb;
  real_T c3_i_delta;
  real_T c3_kb_nargin = 2.0;
  real_T c3_kb_nargout = 1.0;
  real_T c3_f_qa;
  real_T c3_d26;
  real_T c3_d27;
  real_T c3_lb_nargin = 0.0;
  real_T c3_lb_nargout = 1.0;
  real_T c3_o_q;
  real_T c3_d28;
  real_T c3_d29;
  real_T c3_mb_nargin = 0.0;
  real_T c3_mb_nargout = 0.0;
  real_T c3_nb_nargin = 0.0;
  real_T c3_nb_nargout = 1.0;
  real_T c3_p_q;
  real_T c3_ob_nargin = 0.0;
  real_T c3_ob_nargout = 1.0;
  real_T c3_q_q;
  real_T c3_d30;
  real_T c3_d31;
  real_T c3_pb_nargin = 0.0;
  real_T c3_pb_nargout = 1.0;
  real_T c3_r_q;
  real_T c3_qb_nargin = 0.0;
  real_T c3_qb_nargout = 1.0;
  real_T c3_s_q;
  real_T c3_d32;
  real_T c3_d33;
  real_T *c3_mode_out;
  real_T *c3_h1;
  real_T *c3_h2;
  real_T *c3_Qin;
  real_T *c3_h1_out;
  real_T *c3_Qa;
  real_T *c3_h2_out;
  real_T *c3_Qout;
  real_T *c3_x1;
  real_T *c3_x2;
  real_T *c3_flow;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T guard5 = false;
  boolean_T guard6 = false;
  boolean_T guard7 = false;
  boolean_T guard8 = false;
  boolean_T guard9 = false;
  boolean_T guard10 = false;
  boolean_T guard11 = false;
  boolean_T guard12 = false;
  c3_flow = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c3_x2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c3_x1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c3_Qout = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c3_h2_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c3_Qa = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c3_h1_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c3_Qin = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c3_h2 = (real_T *)(ssGetContStates(chartInstance->S) + 1);
  c3_h1 = (real_T *)(ssGetContStates(chartInstance->S) + 0);
  c3_mode_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c3_set_sim_state_side_effects_c3_hw8_two_tanks(chartInstance);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  if (ssIsMajorTimeStep(chartInstance->S) != 0) {
    chartInstance->c3_lastMajorTime = _sfTime_;
    chartInstance->c3_stateChanged = false;
    _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_h1, 1U);
    _SFD_DATA_RANGE_CHECK(chartInstance->c3_h, 2U);
    _SFD_DATA_RANGE_CHECK(*c3_h2, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_DATA_RANGE_CHECK(chartInstance->c3_A, 8U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    _SFD_DATA_RANGE_CHECK(*c3_x1, 10U);
    _SFD_DATA_RANGE_CHECK(*c3_x2, 11U);
    _SFD_DATA_RANGE_CHECK(chartInstance->c3_rho, 12U);
    _SFD_DATA_RANGE_CHECK(chartInstance->c3_g, 13U);
    _SFD_DATA_RANGE_CHECK(chartInstance->c3_Kin, 14U);
    _SFD_DATA_RANGE_CHECK(chartInstance->c3_Ka, 15U);
    _SFD_DATA_RANGE_CHECK(chartInstance->c3_Kout, 16U);
    _SFD_DATA_RANGE_CHECK(*c3_flow, 17U);
    chartInstance->c3_sfEvent = CALL_EVENT;
    _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
    if (chartInstance->c3_is_active_c3_hw8_two_tanks == 0U) {
      _SFD_CC_CALL(CHART_ENTER_ENTRY_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
      chartInstance->c3_stateChanged = true;
      chartInstance->c3_is_active_c3_hw8_two_tanks = 1U;
      _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
      _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 16U, chartInstance->c3_sfEvent);
      _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_i_debug_family_names,
        c3_debug_family_var_map);
      _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 0U, c3_sf_marshallOut,
        c3_sf_marshallIn);
      _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 1U, c3_sf_marshallOut,
        c3_sf_marshallIn);
      *c3_h1 = *c3_x1;
      c3_updateDataWrittenToVector(chartInstance, 1U);
      _SFD_DATA_RANGE_CHECK(*c3_h1, 1U);
      c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
      sf_mex_printf("%s =\\n", "h1");
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "disp", 0U, 1U, 14,
                        c3_c_emlrt_marshallOut(chartInstance, *c3_h1));
      *c3_h2 = *c3_x2;
      c3_updateDataWrittenToVector(chartInstance, 2U);
      _SFD_DATA_RANGE_CHECK(*c3_h2, 3U);
      c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
      sf_mex_printf("%s =\\n", "h2");
      sf_mex_call_debug(sfGlobalDebugInstanceStruct, "disp", 0U, 1U, 14,
                        c3_c_emlrt_marshallOut(chartInstance, *c3_h2));
      _SFD_SYMBOL_SCOPE_POP();
      chartInstance->c3_stateChanged = true;
      chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_Separated;
      _SFD_CS_CALL(STATE_ACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
      chartInstance->c3_tp_Separated = 1U;
    } else {
      switch (chartInstance->c3_is_c3_hw8_two_tanks) {
       case c3_IN_Balancing:
        CV_CHART_EVAL(0, 0, 1);
        _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 0U,
                     chartInstance->c3_sfEvent);
        _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 11U, chartInstance->c3_sfEvent);
        _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 0U,
                     chartInstance->c3_sfEvent);
        _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_j_debug_family_names,
          c3_b_debug_family_var_map);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargin, 0U, c3_sf_marshallOut,
          c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargout, 1U,
          c3_sf_marshallOut, c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_out, 2U, c3_b_sf_marshallOut,
          c3_b_sf_marshallIn);
        c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
        guard12 = false;
        if (CV_EML_COND(0, 0, 0, *c3_h1 < c3_const_h)) {
          if (CV_EML_COND(0, 0, 1, *c3_h2 < c3_const_h)) {
            CV_EML_MCDC(0, 0, 0, true);
            CV_EML_IF(0, 0, 0, true);
            c3_out = true;
          } else {
            guard12 = true;
          }
        } else {
          c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
          guard12 = true;
        }

        if (guard12 == true) {
          CV_EML_MCDC(0, 0, 0, false);
          CV_EML_IF(0, 0, 0, false);
          c3_out = false;
        }

        _SFD_SYMBOL_SCOPE_POP();
        if (c3_out) {
          _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
          chartInstance->c3_tp_Balancing = 0U;
          _SFD_CS_CALL(STATE_INACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
          chartInstance->c3_stateChanged = true;
          chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_Separated;
          _SFD_CS_CALL(STATE_ACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
          chartInstance->c3_tp_Separated = 1U;
        } else {
          _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 10U, chartInstance->c3_sfEvent);
          _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 1U,
                       chartInstance->c3_sfEvent);
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_k_debug_family_names,
            c3_b_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_nargin, 0U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_nargout, 1U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_out, 2U,
            c3_b_sf_marshallOut, c3_b_sf_marshallIn);
          c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
          guard11 = false;
          if (CV_EML_COND(1, 0, 0, *c3_h1 > c3_const_h)) {
            if (CV_EML_COND(1, 0, 1, *c3_h2 < c3_const_h)) {
              CV_EML_MCDC(1, 0, 0, true);
              CV_EML_IF(1, 0, 0, true);
              c3_b_out = true;
            } else {
              guard11 = true;
            }
          } else {
            c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
            guard11 = true;
          }

          if (guard11 == true) {
            CV_EML_MCDC(1, 0, 0, false);
            CV_EML_IF(1, 0, 0, false);
            c3_b_out = false;
          }

          _SFD_SYMBOL_SCOPE_POP();
          if (c3_b_out) {
            _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
            chartInstance->c3_tp_Balancing = 0U;
            _SFD_CS_CALL(STATE_INACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
            chartInstance->c3_stateChanged = true;
            chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_From_tank_1;
            _SFD_CS_CALL(STATE_ACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
            chartInstance->c3_tp_From_tank_1 = 1U;
          } else {
            _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 12U, chartInstance->c3_sfEvent);
            _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 3U,
                         chartInstance->c3_sfEvent);
            _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_m_debug_family_names,
              c3_b_debug_family_var_map);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_nargin, 0U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_nargout, 1U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_out, 2U,
              c3_b_sf_marshallOut, c3_b_sf_marshallIn);
            c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
            guard10 = false;
            if (CV_EML_COND(3, 0, 0, *c3_h1 < c3_const_h)) {
              if (CV_EML_COND(3, 0, 1, *c3_h2 > c3_const_h)) {
                CV_EML_MCDC(3, 0, 0, true);
                CV_EML_IF(3, 0, 0, true);
                c3_c_out = true;
              } else {
                guard10 = true;
              }
            } else {
              c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
              guard10 = true;
            }

            if (guard10 == true) {
              CV_EML_MCDC(3, 0, 0, false);
              CV_EML_IF(3, 0, 0, false);
              c3_c_out = false;
            }

            _SFD_SYMBOL_SCOPE_POP();
            if (c3_c_out) {
              _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
              chartInstance->c3_tp_Balancing = 0U;
              _SFD_CS_CALL(STATE_INACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
              chartInstance->c3_stateChanged = true;
              chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_From_tank_2;
              _SFD_CS_CALL(STATE_ACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
              chartInstance->c3_tp_From_tank_2 = 1U;
            } else {
              _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U,
                           chartInstance->c3_sfEvent);
            }
          }
        }
        break;

       case c3_IN_From_tank_1:
        CV_CHART_EVAL(0, 0, 2);
        _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 1U,
                     chartInstance->c3_sfEvent);
        _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 7U, chartInstance->c3_sfEvent);
        _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 0U,
                     chartInstance->c3_sfEvent);
        _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_j_debug_family_names,
          c3_b_debug_family_var_map);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_nargin, 0U, c3_sf_marshallOut,
          c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_nargout, 1U,
          c3_sf_marshallOut, c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_out, 2U, c3_b_sf_marshallOut,
          c3_b_sf_marshallIn);
        c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
        guard9 = false;
        if (CV_EML_COND(0, 0, 0, *c3_h1 < c3_const_h)) {
          if (CV_EML_COND(0, 0, 1, *c3_h2 < c3_const_h)) {
            CV_EML_MCDC(0, 0, 0, true);
            CV_EML_IF(0, 0, 0, true);
            c3_d_out = true;
          } else {
            guard9 = true;
          }
        } else {
          c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
          guard9 = true;
        }

        if (guard9 == true) {
          CV_EML_MCDC(0, 0, 0, false);
          CV_EML_IF(0, 0, 0, false);
          c3_d_out = false;
        }

        _SFD_SYMBOL_SCOPE_POP();
        if (c3_d_out) {
          _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
          chartInstance->c3_tp_From_tank_1 = 0U;
          _SFD_CS_CALL(STATE_INACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
          chartInstance->c3_stateChanged = true;
          chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_Separated;
          _SFD_CS_CALL(STATE_ACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
          chartInstance->c3_tp_Separated = 1U;
        } else {
          _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 8U, chartInstance->c3_sfEvent);
          _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 3U,
                       chartInstance->c3_sfEvent);
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_m_debug_family_names,
            c3_b_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_nargin, 0U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_nargout, 1U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_out, 2U,
            c3_b_sf_marshallOut, c3_b_sf_marshallIn);
          c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
          guard8 = false;
          if (CV_EML_COND(3, 0, 0, *c3_h1 < c3_const_h)) {
            if (CV_EML_COND(3, 0, 1, *c3_h2 > c3_const_h)) {
              CV_EML_MCDC(3, 0, 0, true);
              CV_EML_IF(3, 0, 0, true);
              c3_e_out = true;
            } else {
              guard8 = true;
            }
          } else {
            c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
            guard8 = true;
          }

          if (guard8 == true) {
            CV_EML_MCDC(3, 0, 0, false);
            CV_EML_IF(3, 0, 0, false);
            c3_e_out = false;
          }

          _SFD_SYMBOL_SCOPE_POP();
          if (c3_e_out) {
            _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
            chartInstance->c3_tp_From_tank_1 = 0U;
            _SFD_CS_CALL(STATE_INACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
            chartInstance->c3_stateChanged = true;
            chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_From_tank_2;
            _SFD_CS_CALL(STATE_ACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
            chartInstance->c3_tp_From_tank_2 = 1U;
          } else {
            _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 9U, chartInstance->c3_sfEvent);
            _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 2U,
                         chartInstance->c3_sfEvent);
            _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_l_debug_family_names,
              c3_b_debug_family_var_map);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_nargin, 0U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_nargout, 1U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_out, 2U,
              c3_b_sf_marshallOut, c3_b_sf_marshallIn);
            c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
            guard7 = false;
            if (CV_EML_COND(2, 0, 0, *c3_h1 > c3_const_h)) {
              if (CV_EML_COND(2, 0, 1, *c3_h2 > c3_const_h)) {
                CV_EML_MCDC(2, 0, 0, true);
                CV_EML_IF(2, 0, 0, true);
                c3_f_out = true;
              } else {
                guard7 = true;
              }
            } else {
              c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
              guard7 = true;
            }

            if (guard7 == true) {
              CV_EML_MCDC(2, 0, 0, false);
              CV_EML_IF(2, 0, 0, false);
              c3_f_out = false;
            }

            _SFD_SYMBOL_SCOPE_POP();
            if (c3_f_out) {
              _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
              chartInstance->c3_tp_From_tank_1 = 0U;
              _SFD_CS_CALL(STATE_INACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
              chartInstance->c3_stateChanged = true;
              chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_Balancing;
              _SFD_CS_CALL(STATE_ACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
              chartInstance->c3_tp_Balancing = 1U;
            } else {
              _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U,
                           chartInstance->c3_sfEvent);
            }
          }
        }
        break;

       case c3_IN_From_tank_2:
        CV_CHART_EVAL(0, 0, 3);
        _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 2U,
                     chartInstance->c3_sfEvent);
        _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 13U, chartInstance->c3_sfEvent);
        _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 0U,
                     chartInstance->c3_sfEvent);
        _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_j_debug_family_names,
          c3_b_debug_family_var_map);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_nargin, 0U, c3_sf_marshallOut,
          c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_nargout, 1U,
          c3_sf_marshallOut, c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_out, 2U, c3_b_sf_marshallOut,
          c3_b_sf_marshallIn);
        c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
        guard6 = false;
        if (CV_EML_COND(0, 0, 0, *c3_h1 < c3_const_h)) {
          if (CV_EML_COND(0, 0, 1, *c3_h2 < c3_const_h)) {
            CV_EML_MCDC(0, 0, 0, true);
            CV_EML_IF(0, 0, 0, true);
            c3_g_out = true;
          } else {
            guard6 = true;
          }
        } else {
          c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
          guard6 = true;
        }

        if (guard6 == true) {
          CV_EML_MCDC(0, 0, 0, false);
          CV_EML_IF(0, 0, 0, false);
          c3_g_out = false;
        }

        _SFD_SYMBOL_SCOPE_POP();
        if (c3_g_out) {
          _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
          chartInstance->c3_tp_From_tank_2 = 0U;
          _SFD_CS_CALL(STATE_INACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
          chartInstance->c3_stateChanged = true;
          chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_Separated;
          _SFD_CS_CALL(STATE_ACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
          chartInstance->c3_tp_Separated = 1U;
        } else {
          _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 14U, chartInstance->c3_sfEvent);
          _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 1U,
                       chartInstance->c3_sfEvent);
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_k_debug_family_names,
            c3_b_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_nargin, 0U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_nargout, 1U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_out, 2U,
            c3_b_sf_marshallOut, c3_b_sf_marshallIn);
          c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
          guard5 = false;
          if (CV_EML_COND(1, 0, 0, *c3_h1 > c3_const_h)) {
            if (CV_EML_COND(1, 0, 1, *c3_h2 < c3_const_h)) {
              CV_EML_MCDC(1, 0, 0, true);
              CV_EML_IF(1, 0, 0, true);
              c3_h_out = true;
            } else {
              guard5 = true;
            }
          } else {
            c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
            guard5 = true;
          }

          if (guard5 == true) {
            CV_EML_MCDC(1, 0, 0, false);
            CV_EML_IF(1, 0, 0, false);
            c3_h_out = false;
          }

          _SFD_SYMBOL_SCOPE_POP();
          if (c3_h_out) {
            _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
            chartInstance->c3_tp_From_tank_2 = 0U;
            _SFD_CS_CALL(STATE_INACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
            chartInstance->c3_stateChanged = true;
            chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_From_tank_1;
            _SFD_CS_CALL(STATE_ACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
            chartInstance->c3_tp_From_tank_1 = 1U;
          } else {
            _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 15U, chartInstance->c3_sfEvent);
            _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 2U,
                         chartInstance->c3_sfEvent);
            _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_l_debug_family_names,
              c3_b_debug_family_var_map);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_nargin, 0U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_nargout, 1U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_out, 2U,
              c3_b_sf_marshallOut, c3_b_sf_marshallIn);
            c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
            guard4 = false;
            if (CV_EML_COND(2, 0, 0, *c3_h1 > c3_const_h)) {
              if (CV_EML_COND(2, 0, 1, *c3_h2 > c3_const_h)) {
                CV_EML_MCDC(2, 0, 0, true);
                CV_EML_IF(2, 0, 0, true);
                c3_i_out = true;
              } else {
                guard4 = true;
              }
            } else {
              c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
              guard4 = true;
            }

            if (guard4 == true) {
              CV_EML_MCDC(2, 0, 0, false);
              CV_EML_IF(2, 0, 0, false);
              c3_i_out = false;
            }

            _SFD_SYMBOL_SCOPE_POP();
            if (c3_i_out) {
              _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
              chartInstance->c3_tp_From_tank_2 = 0U;
              _SFD_CS_CALL(STATE_INACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
              chartInstance->c3_stateChanged = true;
              chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_Balancing;
              _SFD_CS_CALL(STATE_ACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
              chartInstance->c3_tp_Balancing = 1U;
            } else {
              _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 2U,
                           chartInstance->c3_sfEvent);
            }
          }
        }
        break;

       case c3_IN_Separated:
        CV_CHART_EVAL(0, 0, 4);
        _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 3U,
                     chartInstance->c3_sfEvent);
        _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 4U, chartInstance->c3_sfEvent);
        _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 1U,
                     chartInstance->c3_sfEvent);
        _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_k_debug_family_names,
          c3_b_debug_family_var_map);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_nargin, 0U, c3_sf_marshallOut,
          c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_nargout, 1U,
          c3_sf_marshallOut, c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_out, 2U, c3_b_sf_marshallOut,
          c3_b_sf_marshallIn);
        c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
        guard3 = false;
        if (CV_EML_COND(1, 0, 0, *c3_h1 > c3_const_h)) {
          if (CV_EML_COND(1, 0, 1, *c3_h2 < c3_const_h)) {
            CV_EML_MCDC(1, 0, 0, true);
            CV_EML_IF(1, 0, 0, true);
            c3_j_out = true;
          } else {
            guard3 = true;
          }
        } else {
          c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
          guard3 = true;
        }

        if (guard3 == true) {
          CV_EML_MCDC(1, 0, 0, false);
          CV_EML_IF(1, 0, 0, false);
          c3_j_out = false;
        }

        _SFD_SYMBOL_SCOPE_POP();
        if (c3_j_out) {
          _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
          chartInstance->c3_tp_Separated = 0U;
          _SFD_CS_CALL(STATE_INACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
          chartInstance->c3_stateChanged = true;
          chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_From_tank_1;
          _SFD_CS_CALL(STATE_ACTIVE_TAG, 1U, chartInstance->c3_sfEvent);
          chartInstance->c3_tp_From_tank_1 = 1U;
        } else {
          _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 5U, chartInstance->c3_sfEvent);
          _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 2U,
                       chartInstance->c3_sfEvent);
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_l_debug_family_names,
            c3_b_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_nargin, 0U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_nargout, 1U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_out, 2U,
            c3_b_sf_marshallOut, c3_b_sf_marshallIn);
          c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
          guard2 = false;
          if (CV_EML_COND(2, 0, 0, *c3_h1 > c3_const_h)) {
            if (CV_EML_COND(2, 0, 1, *c3_h2 > c3_const_h)) {
              CV_EML_MCDC(2, 0, 0, true);
              CV_EML_IF(2, 0, 0, true);
              c3_k_out = true;
            } else {
              guard2 = true;
            }
          } else {
            c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
            guard2 = true;
          }

          if (guard2 == true) {
            CV_EML_MCDC(2, 0, 0, false);
            CV_EML_IF(2, 0, 0, false);
            c3_k_out = false;
          }

          _SFD_SYMBOL_SCOPE_POP();
          if (c3_k_out) {
            _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
            chartInstance->c3_tp_Separated = 0U;
            _SFD_CS_CALL(STATE_INACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
            chartInstance->c3_stateChanged = true;
            chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_Balancing;
            _SFD_CS_CALL(STATE_ACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
            chartInstance->c3_tp_Balancing = 1U;
          } else {
            _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 6U, chartInstance->c3_sfEvent);
            _SFD_CT_CALL(TRANSITION_BEFORE_PROCESSING_TAG, 3U,
                         chartInstance->c3_sfEvent);
            _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_m_debug_family_names,
              c3_b_debug_family_var_map);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_m_nargin, 0U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_m_nargout, 1U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_out, 2U,
              c3_b_sf_marshallOut, c3_b_sf_marshallIn);
            c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
            guard1 = false;
            if (CV_EML_COND(3, 0, 0, *c3_h1 < c3_const_h)) {
              if (CV_EML_COND(3, 0, 1, *c3_h2 > c3_const_h)) {
                CV_EML_MCDC(3, 0, 0, true);
                CV_EML_IF(3, 0, 0, true);
                c3_l_out = true;
              } else {
                guard1 = true;
              }
            } else {
              c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
              guard1 = true;
            }

            if (guard1 == true) {
              CV_EML_MCDC(3, 0, 0, false);
              CV_EML_IF(3, 0, 0, false);
              c3_l_out = false;
            }

            _SFD_SYMBOL_SCOPE_POP();
            if (c3_l_out) {
              _SFD_CT_CALL(TRANSITION_ACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
              chartInstance->c3_tp_Separated = 0U;
              _SFD_CS_CALL(STATE_INACTIVE_TAG, 3U, chartInstance->c3_sfEvent);
              chartInstance->c3_stateChanged = true;
              chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_From_tank_2;
              _SFD_CS_CALL(STATE_ACTIVE_TAG, 2U, chartInstance->c3_sfEvent);
              chartInstance->c3_tp_From_tank_2 = 1U;
            } else {
              _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 3U,
                           chartInstance->c3_sfEvent);
            }
          }
        }
        break;

       default:
        CV_CHART_EVAL(0, 0, 0);
        chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_NO_ACTIVE_CHILD;
        _SFD_CS_CALL(STATE_INACTIVE_TAG, 0U, chartInstance->c3_sfEvent);
        break;
      }
    }

    _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
    if (chartInstance->c3_stateChanged) {
      ssSetSolverNeedsReset(chartInstance->S);
    }
  }

  _sfTime_ = sf_get_time(chartInstance->S);
  switch (chartInstance->c3_is_c3_hw8_two_tanks) {
   case c3_IN_Balancing:
    CV_CHART_EVAL(0, 0, 1);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_f_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_n_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_n_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_o_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_o_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_hoistedGlobal = *c3_h1;
    c3_b_hoistedGlobal = *c3_h2;
    c3_ha = c3_hoistedGlobal;
    c3_hb = c3_b_hoistedGlobal;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_g_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_q, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(4, 0);
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 3);
    c3_delta = c3_ha - c3_hb;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 4);
    c3_d0 = c3_delta;
    c3_b_sign(chartInstance, &c3_d0);
    c3_d1 = 9810.0 * c3_abs(chartInstance, c3_delta);
    c3_b_sqrt(chartInstance, &c3_d1);
    c3_b_q = c3_d0 * c3_const_Ka * c3_d1;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, -4);
    _SFD_SYMBOL_SCOPE_POP();
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_c_hoistedGlobal = *c3_h1;
    c3_d_hoistedGlobal = *c3_h2;
    c3_b_ha = c3_c_hoistedGlobal;
    c3_b_hb = c3_d_hoistedGlobal;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_g_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_q, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(4, 0);
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 3);
    c3_b_delta = c3_b_ha - c3_b_hb;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 4);
    c3_d2 = c3_b_delta;
    c3_b_sign(chartInstance, &c3_d2);
    c3_d3 = 9810.0 * c3_abs(chartInstance, c3_b_delta);
    c3_b_sqrt(chartInstance, &c3_d3);
    c3_c_q = c3_d2 * c3_const_Ka * c3_d3;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, -4);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d4 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d4);
    c3_d5 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d5);
    c3_d_q = c3_d4 * c3_const_Kout * c3_d5;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_mode_out = 2.0;
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_s_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_s_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_e_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qin = c3_e_q;
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_e_hoistedGlobal = *c3_h1;
    c3_f_hoistedGlobal = *c3_h2;
    c3_c_ha = c3_e_hoistedGlobal;
    c3_c_hb = c3_f_hoistedGlobal;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_g_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_q, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(4, 0);
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 3);
    c3_c_delta = c3_c_ha - c3_c_hb;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 4);
    c3_d6 = c3_c_delta;
    c3_b_sign(chartInstance, &c3_d6);
    c3_d7 = 9810.0 * c3_abs(chartInstance, c3_c_delta);
    c3_b_sqrt(chartInstance, &c3_d7);
    c3_f_q = c3_d6 * c3_const_Ka * c3_d7;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, -4);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qa = c3_f_q;
    c3_updateDataWrittenToVector(chartInstance, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_u_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_u_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d8 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d8);
    c3_d9 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d9);
    c3_g_q = c3_d8 * c3_const_Kout * c3_d9;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qout = c3_g_q;
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    *c3_h1_out = *c3_h1;
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    *c3_h2_out = *c3_h2;
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
    break;

   case c3_IN_From_tank_1:
    CV_CHART_EVAL(0, 0, 2);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_v_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_v_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_w_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_w_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_h_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_g_hoistedGlobal = *c3_h1;
    c3_d_ha = c3_g_hoistedGlobal;
    c3_d_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_x_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_x_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_d_delta = c3_d_ha - c3_d_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d10 = c3_d_delta;
    c3_b_sign(chartInstance, &c3_d10);
    c3_d11 = 9810.0 * c3_abs(chartInstance, c3_d_delta);
    c3_b_sqrt(chartInstance, &c3_d11);
    c3_qa = c3_d10 * c3_const_Ka * c3_d11;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_h_hoistedGlobal = *c3_h1;
    c3_e_ha = c3_h_hoistedGlobal;
    c3_e_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_y_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_y_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_e_delta = c3_e_ha - c3_e_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d12 = c3_e_delta;
    c3_b_sign(chartInstance, &c3_d12);
    c3_d13 = 9810.0 * c3_abs(chartInstance, c3_e_delta);
    c3_b_sqrt(chartInstance, &c3_d13);
    c3_b_qa = c3_d12 * c3_const_Ka * c3_d13;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ab_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ab_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d14 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d14);
    c3_d15 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d15);
    c3_i_q = c3_d14 * c3_const_Kout * c3_d15;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_mode_out = 1.0;
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_bb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_bb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_j_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qin = c3_j_q;
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_i_hoistedGlobal = *c3_h1;
    c3_f_ha = c3_i_hoistedGlobal;
    c3_f_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_cb_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_cb_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_f_delta = c3_f_ha - c3_f_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d16 = c3_f_delta;
    c3_b_sign(chartInstance, &c3_d16);
    c3_d17 = 9810.0 * c3_abs(chartInstance, c3_f_delta);
    c3_b_sqrt(chartInstance, &c3_d17);
    c3_c_qa = c3_d16 * c3_const_Ka * c3_d17;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qa = c3_c_qa;
    c3_updateDataWrittenToVector(chartInstance, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_db_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_db_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d18 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d18);
    c3_d19 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d19);
    c3_k_q = c3_d18 * c3_const_Kout * c3_d19;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qout = c3_k_q;
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    *c3_h1_out = *c3_h1;
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    *c3_h2_out = *c3_h2;
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c3_sfEvent);
    break;

   case c3_IN_From_tank_2:
    CV_CHART_EVAL(0, 0, 3);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_h_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_eb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_eb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_fb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_fb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_l_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_j_hoistedGlobal = *c3_h2;
    c3_g_ha = c3_j_hoistedGlobal;
    c3_g_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_gb_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_gb_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_g_delta = c3_g_ha - c3_g_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d20 = c3_g_delta;
    c3_b_sign(chartInstance, &c3_d20);
    c3_d21 = 9810.0 * c3_abs(chartInstance, c3_g_delta);
    c3_b_sqrt(chartInstance, &c3_d21);
    c3_d_qa = c3_d20 * c3_const_Ka * c3_d21;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_k_hoistedGlobal = *c3_h2;
    c3_h_ha = c3_k_hoistedGlobal;
    c3_h_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_hb_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_hb_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_h_delta = c3_h_ha - c3_h_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d22 = c3_h_delta;
    c3_b_sign(chartInstance, &c3_d22);
    c3_d23 = 9810.0 * c3_abs(chartInstance, c3_h_delta);
    c3_b_sqrt(chartInstance, &c3_d23);
    c3_e_qa = c3_d22 * c3_const_Ka * c3_d23;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ib_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ib_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_m_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d24 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d24);
    c3_d25 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d25);
    c3_m_q = c3_d24 * c3_const_Kout * c3_d25;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_mode_out = -1.0;
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_jb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_jb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_n_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_n_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qin = c3_n_q;
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_l_hoistedGlobal = *c3_h2;
    c3_i_ha = c3_l_hoistedGlobal;
    c3_i_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_kb_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_kb_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_i_delta = c3_i_ha - c3_i_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d26 = c3_i_delta;
    c3_b_sign(chartInstance, &c3_d26);
    c3_d27 = 9810.0 * c3_abs(chartInstance, c3_i_delta);
    c3_b_sqrt(chartInstance, &c3_d27);
    c3_f_qa = c3_d26 * c3_const_Ka * c3_d27;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qa = c3_f_qa;
    c3_updateDataWrittenToVector(chartInstance, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_lb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_lb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_o_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d28 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d28);
    c3_d29 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d29);
    c3_o_q = c3_d28 * c3_const_Kout * c3_d29;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qout = c3_o_q;
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    *c3_h1_out = *c3_h1;
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    *c3_h2_out = *c3_h2;
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
    break;

   case c3_IN_Separated:
    CV_CHART_EVAL(0, 0, 4);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 3U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_e_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_mb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_mb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_p_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ob_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ob_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d30 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d30);
    c3_d31 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d31);
    c3_q_q = c3_d30 * c3_const_Kout * c3_d31;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_mode_out = 0.0;
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_pb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_pb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_r_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qin = c3_r_q;
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_qb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_qb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_s_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d32 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d32);
    c3_d33 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d33);
    c3_s_q = c3_d32 * c3_const_Kout * c3_d33;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qout = c3_s_q;
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    *c3_h1_out = *c3_h1;
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    *c3_h2_out = *c3_h2;
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 3U, chartInstance->c3_sfEvent);
    break;

   default:
    CV_CHART_EVAL(0, 0, 0);
    break;
  }

  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_hw8_two_tanksMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void zeroCrossings_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  uint32_T c3_debug_family_var_map[3];
  real_T c3_nargin = 0.0;
  real_T c3_nargout = 1.0;
  boolean_T c3_out;
  real_T c3_b_nargin = 0.0;
  real_T c3_b_nargout = 1.0;
  boolean_T c3_b_out;
  real_T c3_c_nargin = 0.0;
  real_T c3_c_nargout = 1.0;
  boolean_T c3_c_out;
  real_T c3_d_nargin = 0.0;
  real_T c3_d_nargout = 1.0;
  boolean_T c3_d_out;
  real_T c3_e_nargin = 0.0;
  real_T c3_e_nargout = 1.0;
  boolean_T c3_e_out;
  real_T c3_f_nargin = 0.0;
  real_T c3_f_nargout = 1.0;
  boolean_T c3_f_out;
  real_T c3_g_nargin = 0.0;
  real_T c3_g_nargout = 1.0;
  boolean_T c3_g_out;
  real_T c3_h_nargin = 0.0;
  real_T c3_h_nargout = 1.0;
  boolean_T c3_h_out;
  real_T c3_i_nargin = 0.0;
  real_T c3_i_nargout = 1.0;
  boolean_T c3_i_out;
  real_T c3_j_nargin = 0.0;
  real_T c3_j_nargout = 1.0;
  boolean_T c3_j_out;
  real_T c3_k_nargin = 0.0;
  real_T c3_k_nargout = 1.0;
  boolean_T c3_k_out;
  real_T c3_l_nargin = 0.0;
  real_T c3_l_nargout = 1.0;
  boolean_T c3_l_out;
  real_T *c3_zcVar;
  real_T *c3_h1;
  real_T *c3_h2;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  boolean_T guard4 = false;
  boolean_T guard5 = false;
  boolean_T guard6 = false;
  boolean_T guard7 = false;
  boolean_T guard8 = false;
  boolean_T guard9 = false;
  boolean_T guard10 = false;
  boolean_T guard11 = false;
  boolean_T guard12 = false;
  c3_h2 = (real_T *)(ssGetContStates(chartInstance->S) + 1);
  c3_h1 = (real_T *)(ssGetContStates(chartInstance->S) + 0);
  c3_zcVar = (real_T *)(ssGetNonsampledZCs(chartInstance->S) + 0);
  _sfTime_ = sf_get_time(chartInstance->S);
  if (chartInstance->c3_lastMajorTime == _sfTime_) {
    *c3_zcVar = -1.0;
  } else {
    chartInstance->c3_stateChanged = false;
    if (chartInstance->c3_is_active_c3_hw8_two_tanks == 0U) {
      chartInstance->c3_stateChanged = true;
    } else {
      switch (chartInstance->c3_is_c3_hw8_two_tanks) {
       case c3_IN_Balancing:
        CV_CHART_EVAL(0, 0, 1);
        _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_j_debug_family_names,
          c3_debug_family_var_map);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 0U, c3_sf_marshallOut,
          c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 1U, c3_sf_marshallOut,
          c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_out, 2U, c3_b_sf_marshallOut,
          c3_b_sf_marshallIn);
        c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
        guard12 = false;
        if (CV_EML_COND(0, 0, 0, *c3_h1 < c3_const_h)) {
          if (CV_EML_COND(0, 0, 1, *c3_h2 < c3_const_h)) {
            CV_EML_MCDC(0, 0, 0, true);
            CV_EML_IF(0, 0, 0, true);
            c3_out = true;
          } else {
            guard12 = true;
          }
        } else {
          c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
          guard12 = true;
        }

        if (guard12 == true) {
          CV_EML_MCDC(0, 0, 0, false);
          CV_EML_IF(0, 0, 0, false);
          c3_out = false;
        }

        _SFD_SYMBOL_SCOPE_POP();
        if (c3_out) {
          chartInstance->c3_stateChanged = true;
        } else {
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_k_debug_family_names,
            c3_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargin, 0U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargout, 1U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_out, 2U,
            c3_b_sf_marshallOut, c3_b_sf_marshallIn);
          c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
          guard11 = false;
          if (CV_EML_COND(1, 0, 0, *c3_h1 > c3_const_h)) {
            if (CV_EML_COND(1, 0, 1, *c3_h2 < c3_const_h)) {
              CV_EML_MCDC(1, 0, 0, true);
              CV_EML_IF(1, 0, 0, true);
              c3_b_out = true;
            } else {
              guard11 = true;
            }
          } else {
            c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
            guard11 = true;
          }

          if (guard11 == true) {
            CV_EML_MCDC(1, 0, 0, false);
            CV_EML_IF(1, 0, 0, false);
            c3_b_out = false;
          }

          _SFD_SYMBOL_SCOPE_POP();
          if (c3_b_out) {
            chartInstance->c3_stateChanged = true;
          } else {
            _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_m_debug_family_names,
              c3_debug_family_var_map);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_nargin, 0U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_nargout, 1U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_out, 2U,
              c3_b_sf_marshallOut, c3_b_sf_marshallIn);
            c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
            guard10 = false;
            if (CV_EML_COND(3, 0, 0, *c3_h1 < c3_const_h)) {
              if (CV_EML_COND(3, 0, 1, *c3_h2 > c3_const_h)) {
                CV_EML_MCDC(3, 0, 0, true);
                CV_EML_IF(3, 0, 0, true);
                c3_c_out = true;
              } else {
                guard10 = true;
              }
            } else {
              c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
              guard10 = true;
            }

            if (guard10 == true) {
              CV_EML_MCDC(3, 0, 0, false);
              CV_EML_IF(3, 0, 0, false);
              c3_c_out = false;
            }

            _SFD_SYMBOL_SCOPE_POP();
            if (c3_c_out) {
              chartInstance->c3_stateChanged = true;
            }
          }
        }
        break;

       case c3_IN_From_tank_1:
        CV_CHART_EVAL(0, 0, 2);
        _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_j_debug_family_names,
          c3_debug_family_var_map);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_nargin, 0U, c3_sf_marshallOut,
          c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_nargout, 1U,
          c3_sf_marshallOut, c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_out, 2U, c3_b_sf_marshallOut,
          c3_b_sf_marshallIn);
        c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
        guard9 = false;
        if (CV_EML_COND(0, 0, 0, *c3_h1 < c3_const_h)) {
          if (CV_EML_COND(0, 0, 1, *c3_h2 < c3_const_h)) {
            CV_EML_MCDC(0, 0, 0, true);
            CV_EML_IF(0, 0, 0, true);
            c3_d_out = true;
          } else {
            guard9 = true;
          }
        } else {
          c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
          guard9 = true;
        }

        if (guard9 == true) {
          CV_EML_MCDC(0, 0, 0, false);
          CV_EML_IF(0, 0, 0, false);
          c3_d_out = false;
        }

        _SFD_SYMBOL_SCOPE_POP();
        if (c3_d_out) {
          chartInstance->c3_stateChanged = true;
        } else {
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_m_debug_family_names,
            c3_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_nargin, 0U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_nargout, 1U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_out, 2U,
            c3_b_sf_marshallOut, c3_b_sf_marshallIn);
          c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
          guard8 = false;
          if (CV_EML_COND(3, 0, 0, *c3_h1 < c3_const_h)) {
            if (CV_EML_COND(3, 0, 1, *c3_h2 > c3_const_h)) {
              CV_EML_MCDC(3, 0, 0, true);
              CV_EML_IF(3, 0, 0, true);
              c3_e_out = true;
            } else {
              guard8 = true;
            }
          } else {
            c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
            guard8 = true;
          }

          if (guard8 == true) {
            CV_EML_MCDC(3, 0, 0, false);
            CV_EML_IF(3, 0, 0, false);
            c3_e_out = false;
          }

          _SFD_SYMBOL_SCOPE_POP();
          if (c3_e_out) {
            chartInstance->c3_stateChanged = true;
          } else {
            _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_l_debug_family_names,
              c3_debug_family_var_map);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_nargin, 0U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_nargout, 1U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_out, 2U,
              c3_b_sf_marshallOut, c3_b_sf_marshallIn);
            c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
            guard7 = false;
            if (CV_EML_COND(2, 0, 0, *c3_h1 > c3_const_h)) {
              if (CV_EML_COND(2, 0, 1, *c3_h2 > c3_const_h)) {
                CV_EML_MCDC(2, 0, 0, true);
                CV_EML_IF(2, 0, 0, true);
                c3_f_out = true;
              } else {
                guard7 = true;
              }
            } else {
              c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
              guard7 = true;
            }

            if (guard7 == true) {
              CV_EML_MCDC(2, 0, 0, false);
              CV_EML_IF(2, 0, 0, false);
              c3_f_out = false;
            }

            _SFD_SYMBOL_SCOPE_POP();
            if (c3_f_out) {
              chartInstance->c3_stateChanged = true;
            }
          }
        }
        break;

       case c3_IN_From_tank_2:
        CV_CHART_EVAL(0, 0, 3);
        _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_j_debug_family_names,
          c3_debug_family_var_map);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_nargin, 0U, c3_sf_marshallOut,
          c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_nargout, 1U,
          c3_sf_marshallOut, c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_out, 2U, c3_b_sf_marshallOut,
          c3_b_sf_marshallIn);
        c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
        guard6 = false;
        if (CV_EML_COND(0, 0, 0, *c3_h1 < c3_const_h)) {
          if (CV_EML_COND(0, 0, 1, *c3_h2 < c3_const_h)) {
            CV_EML_MCDC(0, 0, 0, true);
            CV_EML_IF(0, 0, 0, true);
            c3_g_out = true;
          } else {
            guard6 = true;
          }
        } else {
          c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
          guard6 = true;
        }

        if (guard6 == true) {
          CV_EML_MCDC(0, 0, 0, false);
          CV_EML_IF(0, 0, 0, false);
          c3_g_out = false;
        }

        _SFD_SYMBOL_SCOPE_POP();
        if (c3_g_out) {
          chartInstance->c3_stateChanged = true;
        } else {
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_k_debug_family_names,
            c3_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_nargin, 0U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_nargout, 1U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_out, 2U,
            c3_b_sf_marshallOut, c3_b_sf_marshallIn);
          c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
          guard5 = false;
          if (CV_EML_COND(1, 0, 0, *c3_h1 > c3_const_h)) {
            if (CV_EML_COND(1, 0, 1, *c3_h2 < c3_const_h)) {
              CV_EML_MCDC(1, 0, 0, true);
              CV_EML_IF(1, 0, 0, true);
              c3_h_out = true;
            } else {
              guard5 = true;
            }
          } else {
            c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
            guard5 = true;
          }

          if (guard5 == true) {
            CV_EML_MCDC(1, 0, 0, false);
            CV_EML_IF(1, 0, 0, false);
            c3_h_out = false;
          }

          _SFD_SYMBOL_SCOPE_POP();
          if (c3_h_out) {
            chartInstance->c3_stateChanged = true;
          } else {
            _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_l_debug_family_names,
              c3_debug_family_var_map);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_nargin, 0U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_nargout, 1U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_out, 2U,
              c3_b_sf_marshallOut, c3_b_sf_marshallIn);
            c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
            guard4 = false;
            if (CV_EML_COND(2, 0, 0, *c3_h1 > c3_const_h)) {
              if (CV_EML_COND(2, 0, 1, *c3_h2 > c3_const_h)) {
                CV_EML_MCDC(2, 0, 0, true);
                CV_EML_IF(2, 0, 0, true);
                c3_i_out = true;
              } else {
                guard4 = true;
              }
            } else {
              c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
              guard4 = true;
            }

            if (guard4 == true) {
              CV_EML_MCDC(2, 0, 0, false);
              CV_EML_IF(2, 0, 0, false);
              c3_i_out = false;
            }

            _SFD_SYMBOL_SCOPE_POP();
            if (c3_i_out) {
              chartInstance->c3_stateChanged = true;
            }
          }
        }
        break;

       case c3_IN_Separated:
        CV_CHART_EVAL(0, 0, 4);
        _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_k_debug_family_names,
          c3_debug_family_var_map);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_nargin, 0U, c3_sf_marshallOut,
          c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_nargout, 1U,
          c3_sf_marshallOut, c3_sf_marshallIn);
        _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_out, 2U, c3_b_sf_marshallOut,
          c3_b_sf_marshallIn);
        c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
        guard3 = false;
        if (CV_EML_COND(1, 0, 0, *c3_h1 > c3_const_h)) {
          if (CV_EML_COND(1, 0, 1, *c3_h2 < c3_const_h)) {
            CV_EML_MCDC(1, 0, 0, true);
            CV_EML_IF(1, 0, 0, true);
            c3_j_out = true;
          } else {
            guard3 = true;
          }
        } else {
          c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
          guard3 = true;
        }

        if (guard3 == true) {
          CV_EML_MCDC(1, 0, 0, false);
          CV_EML_IF(1, 0, 0, false);
          c3_j_out = false;
        }

        _SFD_SYMBOL_SCOPE_POP();
        if (c3_j_out) {
          chartInstance->c3_stateChanged = true;
        } else {
          _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_l_debug_family_names,
            c3_debug_family_var_map);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_nargin, 0U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_nargout, 1U,
            c3_sf_marshallOut, c3_sf_marshallIn);
          _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_out, 2U,
            c3_b_sf_marshallOut, c3_b_sf_marshallIn);
          c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
          guard2 = false;
          if (CV_EML_COND(2, 0, 0, *c3_h1 > c3_const_h)) {
            if (CV_EML_COND(2, 0, 1, *c3_h2 > c3_const_h)) {
              CV_EML_MCDC(2, 0, 0, true);
              CV_EML_IF(2, 0, 0, true);
              c3_k_out = true;
            } else {
              guard2 = true;
            }
          } else {
            c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
            guard2 = true;
          }

          if (guard2 == true) {
            CV_EML_MCDC(2, 0, 0, false);
            CV_EML_IF(2, 0, 0, false);
            c3_k_out = false;
          }

          _SFD_SYMBOL_SCOPE_POP();
          if (c3_k_out) {
            chartInstance->c3_stateChanged = true;
          } else {
            _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_m_debug_family_names,
              c3_debug_family_var_map);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_nargin, 0U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_nargout, 1U,
              c3_sf_marshallOut, c3_sf_marshallIn);
            _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_out, 2U,
              c3_b_sf_marshallOut, c3_b_sf_marshallIn);
            c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
            guard1 = false;
            if (CV_EML_COND(3, 0, 0, *c3_h1 < c3_const_h)) {
              if (CV_EML_COND(3, 0, 1, *c3_h2 > c3_const_h)) {
                CV_EML_MCDC(3, 0, 0, true);
                CV_EML_IF(3, 0, 0, true);
                c3_l_out = true;
              } else {
                guard1 = true;
              }
            } else {
              c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
              guard1 = true;
            }

            if (guard1 == true) {
              CV_EML_MCDC(3, 0, 0, false);
              CV_EML_IF(3, 0, 0, false);
              c3_l_out = false;
            }

            _SFD_SYMBOL_SCOPE_POP();
            if (c3_l_out) {
              chartInstance->c3_stateChanged = true;
            }
          }
        }
        break;

       default:
        CV_CHART_EVAL(0, 0, 0);
        chartInstance->c3_is_c3_hw8_two_tanks = c3_IN_NO_ACTIVE_CHILD;
        break;
      }
    }

    if (chartInstance->c3_stateChanged) {
      *c3_zcVar = 1.0;
    } else {
      *c3_zcVar = -1.0;
    }
  }
}

static void derivatives_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  uint32_T c3_debug_family_var_map[2];
  real_T c3_nargin = 0.0;
  real_T c3_nargout = 0.0;
  uint32_T c3_b_debug_family_var_map[3];
  real_T c3_b_nargin = 0.0;
  real_T c3_b_nargout = 1.0;
  real_T c3_q;
  real_T c3_hoistedGlobal;
  real_T c3_b_hoistedGlobal;
  real_T c3_ha;
  real_T c3_hb;
  uint32_T c3_c_debug_family_var_map[6];
  real_T c3_delta;
  real_T c3_c_nargin = 2.0;
  real_T c3_c_nargout = 1.0;
  real_T c3_b_q;
  real_T c3_d34;
  real_T c3_d35;
  real_T c3_c_hoistedGlobal;
  real_T c3_d_hoistedGlobal;
  real_T c3_b_ha;
  real_T c3_b_hb;
  real_T c3_b_delta;
  real_T c3_d_nargin = 2.0;
  real_T c3_d_nargout = 1.0;
  real_T c3_c_q;
  real_T c3_d36;
  real_T c3_d37;
  real_T c3_e_nargin = 0.0;
  real_T c3_e_nargout = 1.0;
  real_T c3_d_q;
  real_T c3_d38;
  real_T c3_d39;
  real_T c3_f_nargin = 0.0;
  real_T c3_f_nargout = 1.0;
  real_T c3_e_q;
  real_T c3_e_hoistedGlobal;
  real_T c3_f_hoistedGlobal;
  real_T c3_c_ha;
  real_T c3_c_hb;
  real_T c3_c_delta;
  real_T c3_g_nargin = 2.0;
  real_T c3_g_nargout = 1.0;
  real_T c3_f_q;
  real_T c3_d40;
  real_T c3_d41;
  real_T c3_h_nargin = 0.0;
  real_T c3_h_nargout = 1.0;
  real_T c3_g_q;
  real_T c3_d42;
  real_T c3_d43;
  real_T c3_i_nargin = 0.0;
  real_T c3_i_nargout = 0.0;
  real_T c3_j_nargin = 0.0;
  real_T c3_j_nargout = 1.0;
  real_T c3_h_q;
  real_T c3_g_hoistedGlobal;
  real_T c3_d_ha;
  real_T c3_d_hb;
  real_T c3_d_delta;
  real_T c3_k_nargin = 2.0;
  real_T c3_k_nargout = 1.0;
  real_T c3_qa;
  real_T c3_d44;
  real_T c3_d45;
  real_T c3_h_hoistedGlobal;
  real_T c3_e_ha;
  real_T c3_e_hb;
  real_T c3_e_delta;
  real_T c3_l_nargin = 2.0;
  real_T c3_l_nargout = 1.0;
  real_T c3_b_qa;
  real_T c3_d46;
  real_T c3_d47;
  real_T c3_m_nargin = 0.0;
  real_T c3_m_nargout = 1.0;
  real_T c3_i_q;
  real_T c3_d48;
  real_T c3_d49;
  real_T c3_n_nargin = 0.0;
  real_T c3_n_nargout = 1.0;
  real_T c3_j_q;
  real_T c3_i_hoistedGlobal;
  real_T c3_f_ha;
  real_T c3_f_hb;
  real_T c3_f_delta;
  real_T c3_o_nargin = 2.0;
  real_T c3_o_nargout = 1.0;
  real_T c3_c_qa;
  real_T c3_d50;
  real_T c3_d51;
  real_T c3_p_nargin = 0.0;
  real_T c3_p_nargout = 1.0;
  real_T c3_k_q;
  real_T c3_d52;
  real_T c3_d53;
  real_T c3_q_nargin = 0.0;
  real_T c3_q_nargout = 0.0;
  real_T c3_r_nargin = 0.0;
  real_T c3_r_nargout = 1.0;
  real_T c3_l_q;
  real_T c3_j_hoistedGlobal;
  real_T c3_g_ha;
  real_T c3_g_hb;
  real_T c3_g_delta;
  real_T c3_s_nargin = 2.0;
  real_T c3_s_nargout = 1.0;
  real_T c3_d_qa;
  real_T c3_d54;
  real_T c3_d55;
  real_T c3_k_hoistedGlobal;
  real_T c3_h_ha;
  real_T c3_h_hb;
  real_T c3_h_delta;
  real_T c3_t_nargin = 2.0;
  real_T c3_t_nargout = 1.0;
  real_T c3_e_qa;
  real_T c3_d56;
  real_T c3_d57;
  real_T c3_u_nargin = 0.0;
  real_T c3_u_nargout = 1.0;
  real_T c3_m_q;
  real_T c3_d58;
  real_T c3_d59;
  real_T c3_v_nargin = 0.0;
  real_T c3_v_nargout = 1.0;
  real_T c3_n_q;
  real_T c3_l_hoistedGlobal;
  real_T c3_i_ha;
  real_T c3_i_hb;
  real_T c3_i_delta;
  real_T c3_w_nargin = 2.0;
  real_T c3_w_nargout = 1.0;
  real_T c3_f_qa;
  real_T c3_d60;
  real_T c3_d61;
  real_T c3_x_nargin = 0.0;
  real_T c3_x_nargout = 1.0;
  real_T c3_o_q;
  real_T c3_d62;
  real_T c3_d63;
  real_T c3_y_nargin = 0.0;
  real_T c3_y_nargout = 0.0;
  real_T c3_ab_nargin = 0.0;
  real_T c3_ab_nargout = 1.0;
  real_T c3_p_q;
  real_T c3_bb_nargin = 0.0;
  real_T c3_bb_nargout = 1.0;
  real_T c3_q_q;
  real_T c3_d64;
  real_T c3_d65;
  real_T c3_cb_nargin = 0.0;
  real_T c3_cb_nargout = 1.0;
  real_T c3_r_q;
  real_T c3_db_nargin = 0.0;
  real_T c3_db_nargout = 1.0;
  real_T c3_s_q;
  real_T c3_d66;
  real_T c3_d67;
  real_T *c3_h1_dot;
  real_T *c3_h2_dot;
  real_T *c3_mode_out;
  real_T *c3_Qin;
  real_T *c3_Qa;
  real_T *c3_Qout;
  real_T *c3_h1_out;
  real_T *c3_h2_out;
  real_T *c3_flow;
  real_T *c3_h1;
  real_T *c3_h2;
  c3_flow = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c3_Qout = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c3_h2_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c3_Qa = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c3_h1_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c3_Qin = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c3_h2_dot = (real_T *)(ssGetdX(chartInstance->S) + 1);
  c3_h2 = (real_T *)(ssGetContStates(chartInstance->S) + 1);
  c3_h1_dot = (real_T *)(ssGetdX(chartInstance->S) + 0);
  c3_h1 = (real_T *)(ssGetContStates(chartInstance->S) + 0);
  c3_mode_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  *c3_h1_dot = 0.0;
  _SFD_DATA_RANGE_CHECK(*c3_h1_dot, 1U);
  *c3_h2_dot = 0.0;
  _SFD_DATA_RANGE_CHECK(*c3_h2_dot, 3U);
  _sfTime_ = sf_get_time(chartInstance->S);
  switch (chartInstance->c3_is_c3_hw8_two_tanks) {
   case c3_IN_Balancing:
    CV_CHART_EVAL(0, 0, 1);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_f_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_hoistedGlobal = *c3_h1;
    c3_b_hoistedGlobal = *c3_h2;
    c3_ha = c3_hoistedGlobal;
    c3_hb = c3_b_hoistedGlobal;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_g_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_q, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(4, 0);
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 3);
    c3_delta = c3_ha - c3_hb;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 4);
    c3_d34 = c3_delta;
    c3_b_sign(chartInstance, &c3_d34);
    c3_d35 = 9810.0 * c3_abs(chartInstance, c3_delta);
    c3_b_sqrt(chartInstance, &c3_d35);
    c3_b_q = c3_d34 * c3_const_Ka * c3_d35;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, -4);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_h1_dot = c3_mrdivide(chartInstance, c3_q - c3_b_q, c3_const_A);
    _SFD_DATA_RANGE_CHECK(*c3_h1_dot, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_c_hoistedGlobal = *c3_h1;
    c3_d_hoistedGlobal = *c3_h2;
    c3_b_ha = c3_c_hoistedGlobal;
    c3_b_hb = c3_d_hoistedGlobal;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_g_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_q, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(4, 0);
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 3);
    c3_b_delta = c3_b_ha - c3_b_hb;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 4);
    c3_d36 = c3_b_delta;
    c3_b_sign(chartInstance, &c3_d36);
    c3_d37 = 9810.0 * c3_abs(chartInstance, c3_b_delta);
    c3_b_sqrt(chartInstance, &c3_d37);
    c3_c_q = c3_d36 * c3_const_Ka * c3_d37;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, -4);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d38 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d38);
    c3_d39 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d39);
    c3_d_q = c3_d38 * c3_const_Kout * c3_d39;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_h2_dot = c3_mrdivide(chartInstance, c3_c_q - c3_d_q, c3_const_A);
    _SFD_DATA_RANGE_CHECK(*c3_h2_dot, 3U);
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_e_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_e_hoistedGlobal = *c3_h1;
    c3_f_hoistedGlobal = *c3_h2;
    c3_c_ha = c3_e_hoistedGlobal;
    c3_c_hb = c3_f_hoistedGlobal;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_g_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_q, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(4, 0);
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 3);
    c3_c_delta = c3_c_ha - c3_c_hb;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 4);
    c3_d40 = c3_c_delta;
    c3_b_sign(chartInstance, &c3_d40);
    c3_d41 = 9810.0 * c3_abs(chartInstance, c3_c_delta);
    c3_b_sqrt(chartInstance, &c3_d41);
    c3_f_q = c3_d40 * c3_const_Ka * c3_d41;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, -4);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d42 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d42);
    c3_d43 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d43);
    c3_g_q = c3_d42 * c3_const_Kout * c3_d43;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
    break;

   case c3_IN_From_tank_1:
    CV_CHART_EVAL(0, 0, 2);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_h_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_g_hoistedGlobal = *c3_h1;
    c3_d_ha = c3_g_hoistedGlobal;
    c3_d_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_d_delta = c3_d_ha - c3_d_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d44 = c3_d_delta;
    c3_b_sign(chartInstance, &c3_d44);
    c3_d45 = 9810.0 * c3_abs(chartInstance, c3_d_delta);
    c3_b_sqrt(chartInstance, &c3_d45);
    c3_qa = c3_d44 * c3_const_Ka * c3_d45;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_h1_dot = c3_mrdivide(chartInstance, c3_h_q - c3_qa, c3_const_A);
    _SFD_DATA_RANGE_CHECK(*c3_h1_dot, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_h_hoistedGlobal = *c3_h1;
    c3_e_ha = c3_h_hoistedGlobal;
    c3_e_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_e_delta = c3_e_ha - c3_e_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d46 = c3_e_delta;
    c3_b_sign(chartInstance, &c3_d46);
    c3_d47 = 9810.0 * c3_abs(chartInstance, c3_e_delta);
    c3_b_sqrt(chartInstance, &c3_d47);
    c3_b_qa = c3_d46 * c3_const_Ka * c3_d47;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_m_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_m_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d48 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d48);
    c3_d49 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d49);
    c3_i_q = c3_d48 * c3_const_Kout * c3_d49;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_h2_dot = c3_mrdivide(chartInstance, c3_b_qa - c3_i_q, c3_const_A);
    _SFD_DATA_RANGE_CHECK(*c3_h2_dot, 3U);
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_n_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_n_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_j_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_i_hoistedGlobal = *c3_h1;
    c3_f_ha = c3_i_hoistedGlobal;
    c3_f_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_o_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_o_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_f_delta = c3_f_ha - c3_f_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d50 = c3_f_delta;
    c3_b_sign(chartInstance, &c3_d50);
    c3_d51 = 9810.0 * c3_abs(chartInstance, c3_f_delta);
    c3_b_sqrt(chartInstance, &c3_d51);
    c3_c_qa = c3_d50 * c3_const_Ka * c3_d51;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d52 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d52);
    c3_d53 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d53);
    c3_k_q = c3_d52 * c3_const_Kout * c3_d53;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c3_sfEvent);
    break;

   case c3_IN_From_tank_2:
    CV_CHART_EVAL(0, 0, 3);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_h_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_l_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_j_hoistedGlobal = *c3_h2;
    c3_g_ha = c3_j_hoistedGlobal;
    c3_g_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_s_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_s_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_g_delta = c3_g_ha - c3_g_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d54 = c3_g_delta;
    c3_b_sign(chartInstance, &c3_d54);
    c3_d55 = 9810.0 * c3_abs(chartInstance, c3_g_delta);
    c3_b_sqrt(chartInstance, &c3_d55);
    c3_d_qa = c3_d54 * c3_const_Ka * c3_d55;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_h1_dot = c3_mrdivide(chartInstance, c3_l_q - c3_d_qa, c3_const_A);
    _SFD_DATA_RANGE_CHECK(*c3_h1_dot, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_k_hoistedGlobal = *c3_h2;
    c3_h_ha = c3_k_hoistedGlobal;
    c3_h_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_h_delta = c3_h_ha - c3_h_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d56 = c3_h_delta;
    c3_b_sign(chartInstance, &c3_d56);
    c3_d57 = 9810.0 * c3_abs(chartInstance, c3_h_delta);
    c3_b_sqrt(chartInstance, &c3_d57);
    c3_e_qa = c3_d56 * c3_const_Ka * c3_d57;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_u_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_u_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_m_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d58 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d58);
    c3_d59 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d59);
    c3_m_q = c3_d58 * c3_const_Kout * c3_d59;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_h2_dot = c3_mrdivide(chartInstance, c3_e_qa - c3_m_q, c3_const_A);
    _SFD_DATA_RANGE_CHECK(*c3_h2_dot, 3U);
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_v_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_v_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_n_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_n_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_l_hoistedGlobal = *c3_h2;
    c3_i_ha = c3_l_hoistedGlobal;
    c3_i_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_w_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_w_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_i_delta = c3_i_ha - c3_i_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d60 = c3_i_delta;
    c3_b_sign(chartInstance, &c3_d60);
    c3_d61 = 9810.0 * c3_abs(chartInstance, c3_i_delta);
    c3_b_sqrt(chartInstance, &c3_d61);
    c3_f_qa = c3_d60 * c3_const_Ka * c3_d61;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_x_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_x_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_o_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d62 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d62);
    c3_d63 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d63);
    c3_o_q = c3_d62 * c3_const_Kout * c3_d63;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
    break;

   case c3_IN_Separated:
    CV_CHART_EVAL(0, 0, 4);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 3U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_e_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_y_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_y_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ab_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ab_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_p_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_h1_dot = c3_mrdivide(chartInstance, c3_p_q, c3_const_A);
    _SFD_DATA_RANGE_CHECK(*c3_h1_dot, 1U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_bb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_bb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d64 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d64);
    c3_d65 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d65);
    c3_q_q = c3_d64 * c3_const_Kout * c3_d65;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_h2_dot = c3_mrdivide(chartInstance, -c3_q_q, c3_const_A);
    _SFD_DATA_RANGE_CHECK(*c3_h2_dot, 3U);
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_cb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_cb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_r_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_db_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_db_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_s_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d66 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d66);
    c3_d67 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d67);
    c3_s_q = c3_d66 * c3_const_Kout * c3_d67;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 3U, chartInstance->c3_sfEvent);
    break;

   default:
    CV_CHART_EVAL(0, 0, 0);
    break;
  }
}

static void outputs_c3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  uint32_T c3_debug_family_var_map[2];
  real_T c3_nargin = 0.0;
  real_T c3_nargout = 0.0;
  uint32_T c3_b_debug_family_var_map[3];
  real_T c3_b_nargin = 0.0;
  real_T c3_b_nargout = 1.0;
  real_T c3_q;
  real_T c3_hoistedGlobal;
  real_T c3_b_hoistedGlobal;
  real_T c3_ha;
  real_T c3_hb;
  uint32_T c3_c_debug_family_var_map[6];
  real_T c3_delta;
  real_T c3_c_nargin = 2.0;
  real_T c3_c_nargout = 1.0;
  real_T c3_b_q;
  real_T c3_d68;
  real_T c3_d69;
  real_T c3_c_hoistedGlobal;
  real_T c3_d_hoistedGlobal;
  real_T c3_b_ha;
  real_T c3_b_hb;
  real_T c3_b_delta;
  real_T c3_d_nargin = 2.0;
  real_T c3_d_nargout = 1.0;
  real_T c3_c_q;
  real_T c3_d70;
  real_T c3_d71;
  real_T c3_e_nargin = 0.0;
  real_T c3_e_nargout = 1.0;
  real_T c3_d_q;
  real_T c3_d72;
  real_T c3_d73;
  real_T c3_f_nargin = 0.0;
  real_T c3_f_nargout = 1.0;
  real_T c3_e_q;
  real_T c3_e_hoistedGlobal;
  real_T c3_f_hoistedGlobal;
  real_T c3_c_ha;
  real_T c3_c_hb;
  real_T c3_c_delta;
  real_T c3_g_nargin = 2.0;
  real_T c3_g_nargout = 1.0;
  real_T c3_f_q;
  real_T c3_d74;
  real_T c3_d75;
  real_T c3_h_nargin = 0.0;
  real_T c3_h_nargout = 1.0;
  real_T c3_g_q;
  real_T c3_d76;
  real_T c3_d77;
  real_T c3_i_nargin = 0.0;
  real_T c3_i_nargout = 0.0;
  real_T c3_j_nargin = 0.0;
  real_T c3_j_nargout = 1.0;
  real_T c3_h_q;
  real_T c3_g_hoistedGlobal;
  real_T c3_d_ha;
  real_T c3_d_hb;
  real_T c3_d_delta;
  real_T c3_k_nargin = 2.0;
  real_T c3_k_nargout = 1.0;
  real_T c3_qa;
  real_T c3_d78;
  real_T c3_d79;
  real_T c3_h_hoistedGlobal;
  real_T c3_e_ha;
  real_T c3_e_hb;
  real_T c3_e_delta;
  real_T c3_l_nargin = 2.0;
  real_T c3_l_nargout = 1.0;
  real_T c3_b_qa;
  real_T c3_d80;
  real_T c3_d81;
  real_T c3_m_nargin = 0.0;
  real_T c3_m_nargout = 1.0;
  real_T c3_i_q;
  real_T c3_d82;
  real_T c3_d83;
  real_T c3_n_nargin = 0.0;
  real_T c3_n_nargout = 1.0;
  real_T c3_j_q;
  real_T c3_i_hoistedGlobal;
  real_T c3_f_ha;
  real_T c3_f_hb;
  real_T c3_f_delta;
  real_T c3_o_nargin = 2.0;
  real_T c3_o_nargout = 1.0;
  real_T c3_c_qa;
  real_T c3_d84;
  real_T c3_d85;
  real_T c3_p_nargin = 0.0;
  real_T c3_p_nargout = 1.0;
  real_T c3_k_q;
  real_T c3_d86;
  real_T c3_d87;
  real_T c3_q_nargin = 0.0;
  real_T c3_q_nargout = 0.0;
  real_T c3_r_nargin = 0.0;
  real_T c3_r_nargout = 1.0;
  real_T c3_l_q;
  real_T c3_j_hoistedGlobal;
  real_T c3_g_ha;
  real_T c3_g_hb;
  real_T c3_g_delta;
  real_T c3_s_nargin = 2.0;
  real_T c3_s_nargout = 1.0;
  real_T c3_d_qa;
  real_T c3_d88;
  real_T c3_d89;
  real_T c3_k_hoistedGlobal;
  real_T c3_h_ha;
  real_T c3_h_hb;
  real_T c3_h_delta;
  real_T c3_t_nargin = 2.0;
  real_T c3_t_nargout = 1.0;
  real_T c3_e_qa;
  real_T c3_d90;
  real_T c3_d91;
  real_T c3_u_nargin = 0.0;
  real_T c3_u_nargout = 1.0;
  real_T c3_m_q;
  real_T c3_d92;
  real_T c3_d93;
  real_T c3_v_nargin = 0.0;
  real_T c3_v_nargout = 1.0;
  real_T c3_n_q;
  real_T c3_l_hoistedGlobal;
  real_T c3_i_ha;
  real_T c3_i_hb;
  real_T c3_i_delta;
  real_T c3_w_nargin = 2.0;
  real_T c3_w_nargout = 1.0;
  real_T c3_f_qa;
  real_T c3_d94;
  real_T c3_d95;
  real_T c3_x_nargin = 0.0;
  real_T c3_x_nargout = 1.0;
  real_T c3_o_q;
  real_T c3_d96;
  real_T c3_d97;
  real_T c3_y_nargin = 0.0;
  real_T c3_y_nargout = 0.0;
  real_T c3_ab_nargin = 0.0;
  real_T c3_ab_nargout = 1.0;
  real_T c3_p_q;
  real_T c3_bb_nargin = 0.0;
  real_T c3_bb_nargout = 1.0;
  real_T c3_q_q;
  real_T c3_d98;
  real_T c3_d99;
  real_T c3_cb_nargin = 0.0;
  real_T c3_cb_nargout = 1.0;
  real_T c3_r_q;
  real_T c3_db_nargin = 0.0;
  real_T c3_db_nargout = 1.0;
  real_T c3_s_q;
  real_T c3_d100;
  real_T c3_d101;
  real_T *c3_mode_out;
  real_T *c3_Qin;
  real_T *c3_Qa;
  real_T *c3_Qout;
  real_T *c3_h1_out;
  real_T *c3_h1;
  real_T *c3_h2_out;
  real_T *c3_h2;
  real_T *c3_flow;
  c3_flow = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c3_Qout = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c3_h2_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c3_Qa = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c3_h1_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c3_Qin = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c3_h2 = (real_T *)(ssGetContStates(chartInstance->S) + 1);
  c3_h1 = (real_T *)(ssGetContStates(chartInstance->S) + 0);
  c3_mode_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  _sfTime_ = sf_get_time(chartInstance->S);
  switch (chartInstance->c3_is_c3_hw8_two_tanks) {
   case c3_IN_Balancing:
    CV_CHART_EVAL(0, 0, 1);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_f_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_hoistedGlobal = *c3_h1;
    c3_b_hoistedGlobal = *c3_h2;
    c3_ha = c3_hoistedGlobal;
    c3_hb = c3_b_hoistedGlobal;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_g_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_q, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(4, 0);
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 3);
    c3_delta = c3_ha - c3_hb;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 4);
    c3_d68 = c3_delta;
    c3_b_sign(chartInstance, &c3_d68);
    c3_d69 = 9810.0 * c3_abs(chartInstance, c3_delta);
    c3_b_sqrt(chartInstance, &c3_d69);
    c3_b_q = c3_d68 * c3_const_Ka * c3_d69;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, -4);
    _SFD_SYMBOL_SCOPE_POP();
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_c_hoistedGlobal = *c3_h1;
    c3_d_hoistedGlobal = *c3_h2;
    c3_b_ha = c3_c_hoistedGlobal;
    c3_b_hb = c3_d_hoistedGlobal;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_g_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_q, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(4, 0);
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 3);
    c3_b_delta = c3_b_ha - c3_b_hb;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 4);
    c3_d70 = c3_b_delta;
    c3_b_sign(chartInstance, &c3_d70);
    c3_d71 = 9810.0 * c3_abs(chartInstance, c3_b_delta);
    c3_b_sqrt(chartInstance, &c3_d71);
    c3_c_q = c3_d70 * c3_const_Ka * c3_d71;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, -4);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d72 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d72);
    c3_d73 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d73);
    c3_d_q = c3_d72 * c3_const_Kout * c3_d73;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_mode_out = 2.0;
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_e_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qin = c3_e_q;
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_e_hoistedGlobal = *c3_h1;
    c3_f_hoistedGlobal = *c3_h2;
    c3_c_ha = c3_e_hoistedGlobal;
    c3_c_hb = c3_f_hoistedGlobal;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_g_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_q, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(4, 0);
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 3);
    c3_c_delta = c3_c_ha - c3_c_hb;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, 4);
    c3_d74 = c3_c_delta;
    c3_b_sign(chartInstance, &c3_d74);
    c3_d75 = 9810.0 * c3_abs(chartInstance, c3_c_delta);
    c3_b_sqrt(chartInstance, &c3_d75);
    c3_f_q = c3_d74 * c3_const_Ka * c3_d75;
    _SFD_EML_CALL(4U, chartInstance->c3_sfEvent, -4);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qa = c3_f_q;
    c3_updateDataWrittenToVector(chartInstance, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d76 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d76);
    c3_d77 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d77);
    c3_g_q = c3_d76 * c3_const_Kout * c3_d77;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qout = c3_g_q;
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    *c3_h1_out = *c3_h1;
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    *c3_h2_out = *c3_h2;
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, chartInstance->c3_sfEvent);
    break;

   case c3_IN_From_tank_1:
    CV_CHART_EVAL(0, 0, 2);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 1U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_h_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_g_hoistedGlobal = *c3_h1;
    c3_d_ha = c3_g_hoistedGlobal;
    c3_d_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_d_delta = c3_d_ha - c3_d_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d78 = c3_d_delta;
    c3_b_sign(chartInstance, &c3_d78);
    c3_d79 = 9810.0 * c3_abs(chartInstance, c3_d_delta);
    c3_b_sqrt(chartInstance, &c3_d79);
    c3_qa = c3_d78 * c3_const_Ka * c3_d79;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_h_hoistedGlobal = *c3_h1;
    c3_e_ha = c3_h_hoistedGlobal;
    c3_e_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_b_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_e_delta = c3_e_ha - c3_e_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d80 = c3_e_delta;
    c3_b_sign(chartInstance, &c3_d80);
    c3_d81 = 9810.0 * c3_abs(chartInstance, c3_e_delta);
    c3_b_sqrt(chartInstance, &c3_d81);
    c3_b_qa = c3_d80 * c3_const_Ka * c3_d81;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_m_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_m_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d82 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d82);
    c3_d83 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d83);
    c3_i_q = c3_d82 * c3_const_Kout * c3_d83;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_mode_out = 1.0;
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_n_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_n_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_j_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_j_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qin = c3_j_q;
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    c3_i_hoistedGlobal = *c3_h1;
    c3_f_ha = c3_i_hoistedGlobal;
    c3_f_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_o_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_o_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_c_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_f_delta = c3_f_ha - c3_f_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d84 = c3_f_delta;
    c3_b_sign(chartInstance, &c3_d84);
    c3_d85 = 9810.0 * c3_abs(chartInstance, c3_f_delta);
    c3_b_sqrt(chartInstance, &c3_d85);
    c3_c_qa = c3_d84 * c3_const_Ka * c3_d85;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qa = c3_c_qa;
    c3_updateDataWrittenToVector(chartInstance, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_k_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d86 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d86);
    c3_d87 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d87);
    c3_k_q = c3_d86 * c3_const_Kout * c3_d87;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qout = c3_k_q;
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    *c3_h1_out = *c3_h1;
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    *c3_h2_out = *c3_h2;
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, chartInstance->c3_sfEvent);
    break;

   case c3_IN_From_tank_2:
    CV_CHART_EVAL(0, 0, 3);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_h_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_l_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_l_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    c3_j_hoistedGlobal = *c3_h2;
    c3_g_ha = c3_j_hoistedGlobal;
    c3_g_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_s_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_s_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_g_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_d_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_g_delta = c3_g_ha - c3_g_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d88 = c3_g_delta;
    c3_b_sign(chartInstance, &c3_d88);
    c3_d89 = 9810.0 * c3_abs(chartInstance, c3_g_delta);
    c3_b_sqrt(chartInstance, &c3_d89);
    c3_d_qa = c3_d88 * c3_const_Ka * c3_d89;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_k_hoistedGlobal = *c3_h2;
    c3_h_ha = c3_k_hoistedGlobal;
    c3_h_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_t_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_h_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_e_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_h_delta = c3_h_ha - c3_h_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d90 = c3_h_delta;
    c3_b_sign(chartInstance, &c3_d90);
    c3_d91 = 9810.0 * c3_abs(chartInstance, c3_h_delta);
    c3_b_sqrt(chartInstance, &c3_d91);
    c3_e_qa = c3_d90 * c3_const_Ka * c3_d91;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_u_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_u_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_m_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d92 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d92);
    c3_d93 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d93);
    c3_m_q = c3_d92 * c3_const_Kout * c3_d93;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_mode_out = -1.0;
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_v_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_v_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_n_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_n_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qin = c3_n_q;
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_l_hoistedGlobal = *c3_h2;
    c3_i_ha = c3_l_hoistedGlobal;
    c3_i_hb = c3_const_h;
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c3_c_debug_family_names,
      c3_c_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_delta, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_w_nargin, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_w_nargout, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_ha, 3U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_i_hb, 4U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_f_qa, 5U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(5, 0);
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 5);
    c3_i_delta = c3_i_ha - c3_i_hb;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, 6);
    c3_d94 = c3_i_delta;
    c3_b_sign(chartInstance, &c3_d94);
    c3_d95 = 9810.0 * c3_abs(chartInstance, c3_i_delta);
    c3_b_sqrt(chartInstance, &c3_d95);
    c3_f_qa = c3_d94 * c3_const_Ka * c3_d95;
    _SFD_EML_CALL(5U, chartInstance->c3_sfEvent, -6);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qa = c3_f_qa;
    c3_updateDataWrittenToVector(chartInstance, 5U);
    _SFD_DATA_RANGE_CHECK(*c3_Qa, 6U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_x_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_x_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_o_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d96 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d96);
    c3_d97 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d97);
    c3_o_q = c3_d96 * c3_const_Kout * c3_d97;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qout = c3_o_q;
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    *c3_h1_out = *c3_h1;
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    *c3_h2_out = *c3_h2;
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 2U, chartInstance->c3_sfEvent);
    break;

   case c3_IN_Separated:
    CV_CHART_EVAL(0, 0, 4);
    _SFD_CS_CALL(STATE_ENTER_DURING_FUNCTION_TAG, 3U, chartInstance->c3_sfEvent);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 2U, 2U, c3_e_debug_family_names,
      c3_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_y_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_y_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ab_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_ab_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_p_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_p_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_bb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_bb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_q_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d98 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d98);
    c3_d99 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d99);
    c3_q_q = c3_d98 * c3_const_Kout * c3_d99;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_mode_out = 0.0;
    c3_updateDataWrittenToVector(chartInstance, 0U);
    _SFD_DATA_RANGE_CHECK(*c3_mode_out, 0U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_b_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_cb_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_cb_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_r_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(6, 0);
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, 2);
    c3_r_q = c3_const_Kin * *c3_flow;
    _SFD_EML_CALL(6U, chartInstance->c3_sfEvent, -2);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qin = c3_r_q;
    c3_updateDataWrittenToVector(chartInstance, 3U);
    _SFD_DATA_RANGE_CHECK(*c3_Qin, 4U);
    _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 3U, 3U, c3_d_debug_family_names,
      c3_b_debug_family_var_map);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_db_nargin, 0U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_db_nargout, 1U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c3_s_q, 2U, c3_sf_marshallOut,
      c3_sf_marshallIn);
    CV_EML_FCN(7, 0);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, 5);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    c3_d100 = *c3_h2;
    c3_b_sign(chartInstance, &c3_d100);
    c3_d101 = 9810.0 * c3_abs(chartInstance, *c3_h2);
    c3_b_sqrt(chartInstance, &c3_d101);
    c3_s_q = c3_d100 * c3_const_Kout * c3_d101;
    c3_updateDataWrittenToVector(chartInstance, 2U);
    _SFD_EML_CALL(7U, chartInstance->c3_sfEvent, -5);
    _SFD_SYMBOL_SCOPE_POP();
    *c3_Qout = c3_s_q;
    c3_updateDataWrittenToVector(chartInstance, 7U);
    _SFD_DATA_RANGE_CHECK(*c3_Qout, 9U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 1U, 1U);
    *c3_h1_out = *c3_h1;
    c3_updateDataWrittenToVector(chartInstance, 4U);
    _SFD_DATA_RANGE_CHECK(*c3_h1_out, 5U);
    c3_errorIfDataNotWrittenToFcn(chartInstance, 2U, 3U);
    *c3_h2_out = *c3_h2;
    c3_updateDataWrittenToVector(chartInstance, 6U);
    _SFD_DATA_RANGE_CHECK(*c3_h2_out, 7U);
    _SFD_SYMBOL_SCOPE_POP();
    _SFD_CS_CALL(EXIT_OUT_OF_FUNCTION_TAG, 3U, chartInstance->c3_sfEvent);
    break;

   default:
    CV_CHART_EVAL(0, 0, 0);
    break;
  }
}

static void initSimStructsc3_hw8_two_tanks(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c3_eml_ini_fcn_to_be_inlined_23(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c3_eml_term_fcn_to_be_inlined_23(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static real_T c3_sign(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                      c3_x)
{
  real_T c3_b_x;
  c3_b_x = c3_x;
  c3_b_sign(chartInstance, &c3_b_x);
  return c3_b_x;
}

static void c3_isBuiltInNumeric(SFc3_hw8_two_tanksInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c3_eml_scalar_sign(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T c3_x)
{
  real_T c3_b_x;
  c3_b_x = c3_x;
  c3_b_eml_scalar_sign(chartInstance, &c3_b_x);
  return c3_b_x;
}

static real_T c3_abs(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                     c3_x)
{
  return c3_eml_scalar_abs(chartInstance, c3_x);
}

static real_T c3_eml_scalar_abs(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T c3_x)
{
  (void)chartInstance;
  return muDoubleScalarAbs(c3_x);
}

static real_T c3_sqrt(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                      c3_x)
{
  real_T c3_b_x;
  c3_b_x = c3_x;
  c3_b_sqrt(chartInstance, &c3_b_x);
  return c3_b_x;
}

static void c3_eml_error(SFc3_hw8_two_tanksInstanceStruct *chartInstance)
{
  static char_T c3_cv0[4] = { 's', 'q', 'r', 't' };

  static char_T c3_cv1[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  sf_mex_call_debug(sfGlobalDebugInstanceStruct, "error", 0U, 1U, 14,
                    sf_mex_call_debug(sfGlobalDebugInstanceStruct, "message", 1U,
    2U, 14, c3_emlrt_marshallOut(chartInstance, c3_cv1), 14,
    c3_b_emlrt_marshallOut(chartInstance, c3_cv0)));
}

static real_T c3_eml_scalar_sqrt(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T c3_x)
{
  real_T c3_b_x;
  c3_b_x = c3_x;
  c3_b_eml_scalar_sqrt(chartInstance, &c3_b_x);
  return c3_b_x;
}

static real_T c3_mrdivide(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T c3_b_A, real_T c3_B)
{
  return c3_rdivide(chartInstance, c3_b_A, c3_B);
}

static void c3_assert(SFc3_hw8_two_tanksInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static real_T c3_rdivide(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
  c3_x, real_T c3_y)
{
  return c3_eml_div(chartInstance, c3_x, c3_y);
}

static void c3_eml_scalexp_compatible(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static real_T c3_eml_div(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
  c3_x, real_T c3_y)
{
  return c3_div(chartInstance, c3_x, c3_y);
}

static real_T c3_div(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                     c3_x, real_T c3_y)
{
  (void)chartInstance;
  return c3_x / c3_y;
}

static void init_script_number_translation(uint32_T c3_machineNumber, uint32_T
  c3_chartNumber, uint32_T c3_instanceNumber)
{
  (void)c3_machineNumber;
  (void)c3_chartNumber;
  (void)c3_instanceNumber;
}

static const mxArray *c3_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const char_T c3_u[30])
{
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 30), false);
  return c3_y;
}

static const mxArray *c3_b_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const char_T c3_u[4])
{
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 10, 0U, 1U, 0U, 2, 1, 4), false);
  return c3_y;
}

static const mxArray *c3_c_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const real_T c3_u)
{
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 0, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static const mxArray *c3_sf_marshallOut(void *chartInstanceVoid, void *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  SFc3_hw8_two_tanksInstanceStruct *chartInstance;
  chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  sf_mex_assign(&c3_mxArrayOutData, c3_c_emlrt_marshallOut(chartInstance,
    *(real_T *)c3_inData), false);
  return c3_mxArrayOutData;
}

static real_T c3_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_nargout, const char_T *c3_identifier)
{
  real_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_nargout), &c3_thisId);
  sf_mex_destroy(&c3_nargout);
  return c3_y;
}

static real_T c3_b_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  real_T c3_y;
  real_T c3_d102;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_d102, 1, 0, 0U, 0, 0U, 0);
  c3_y = c3_d102;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  SFc3_hw8_two_tanksInstanceStruct *chartInstance;
  chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)chartInstanceVoid;
  *(real_T *)c3_outData = c3_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_mxArrayInData), c3_varName);
  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_d_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const boolean_T c3_u)
{
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 11, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static const mxArray *c3_b_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  SFc3_hw8_two_tanksInstanceStruct *chartInstance;
  chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  sf_mex_assign(&c3_mxArrayOutData, c3_d_emlrt_marshallOut(chartInstance,
    *(boolean_T *)c3_inData), false);
  return c3_mxArrayOutData;
}

static boolean_T c3_c_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_sf_internal_predicateOutput, const char_T
  *c3_identifier)
{
  boolean_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_d_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_sf_internal_predicateOutput), &c3_thisId);
  sf_mex_destroy(&c3_sf_internal_predicateOutput);
  return c3_y;
}

static boolean_T c3_d_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  boolean_T c3_y;
  boolean_T c3_b0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_b0, 1, 11, 0U, 0, 0U, 0);
  c3_y = c3_b0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  SFc3_hw8_two_tanksInstanceStruct *chartInstance;
  chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)chartInstanceVoid;
  *(boolean_T *)c3_outData = c3_c_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_mxArrayInData), c3_varName);
  sf_mex_destroy(&c3_mxArrayInData);
}

const mxArray *sf_c3_hw8_two_tanks_get_eml_resolved_functions_info(void)
{
  const mxArray *c3_nameCaptureInfo = NULL;
  c3_nameCaptureInfo = NULL;
  sf_mex_assign(&c3_nameCaptureInfo, sf_mex_createstruct("structure", 2, 16, 1),
                false);
  c3_info_helper(&c3_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c3_nameCaptureInfo);
  return c3_nameCaptureInfo;
}

static void c3_info_helper(const mxArray **c3_info)
{
  const mxArray *c3_rhs0 = NULL;
  const mxArray *c3_lhs0 = NULL;
  const mxArray *c3_rhs1 = NULL;
  const mxArray *c3_lhs1 = NULL;
  const mxArray *c3_rhs2 = NULL;
  const mxArray *c3_lhs2 = NULL;
  const mxArray *c3_rhs3 = NULL;
  const mxArray *c3_lhs3 = NULL;
  const mxArray *c3_rhs4 = NULL;
  const mxArray *c3_lhs4 = NULL;
  const mxArray *c3_rhs5 = NULL;
  const mxArray *c3_lhs5 = NULL;
  const mxArray *c3_rhs6 = NULL;
  const mxArray *c3_lhs6 = NULL;
  const mxArray *c3_rhs7 = NULL;
  const mxArray *c3_lhs7 = NULL;
  const mxArray *c3_rhs8 = NULL;
  const mxArray *c3_lhs8 = NULL;
  const mxArray *c3_rhs9 = NULL;
  const mxArray *c3_lhs9 = NULL;
  const mxArray *c3_rhs10 = NULL;
  const mxArray *c3_lhs10 = NULL;
  const mxArray *c3_rhs11 = NULL;
  const mxArray *c3_lhs11 = NULL;
  const mxArray *c3_rhs12 = NULL;
  const mxArray *c3_lhs12 = NULL;
  const mxArray *c3_rhs13 = NULL;
  const mxArray *c3_lhs13 = NULL;
  const mxArray *c3_rhs14 = NULL;
  const mxArray *c3_lhs14 = NULL;
  const mxArray *c3_rhs15 = NULL;
  const mxArray *c3_lhs15 = NULL;
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("sign"), "name", "name", 0);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1363735456U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c3_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1363736156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c3_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sign.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("eml_scalar_sign"), "name",
                  "name", 2);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sign.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1356566694U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c3_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(""), "context", "context", 3);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("abs"), "name", "name", 3);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1363735452U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c3_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1363736156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c3_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 5);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1286843912U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c3_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(""), "context", "context", 6);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("sqrt"), "name", "name", 6);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1343855586U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c3_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("eml_error"), "name", "name",
                  7);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1343855558U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c3_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 8);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1286843938U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c3_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(""), "context", "context", 9);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("mrdivide"), "name", "name",
                  9);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1388485296U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1370035086U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c3_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 10);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 10);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1363736156U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c3_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 11);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("rdivide"), "name", "name",
                  11);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 11);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1363735480U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c3_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 12);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1363736156U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c3_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 13);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 13);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1286843996U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c3_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 14);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("eml_div"), "name", "name",
                  14);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1376005888U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c3_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 15);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c3_info, c3_e_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(1389333120U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c3_info, c3_f_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c3_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c3_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c3_info, sf_mex_duplicatearraysafe(&c3_lhs15), "lhs", "lhs",
                  15);
  sf_mex_destroy(&c3_rhs0);
  sf_mex_destroy(&c3_lhs0);
  sf_mex_destroy(&c3_rhs1);
  sf_mex_destroy(&c3_lhs1);
  sf_mex_destroy(&c3_rhs2);
  sf_mex_destroy(&c3_lhs2);
  sf_mex_destroy(&c3_rhs3);
  sf_mex_destroy(&c3_lhs3);
  sf_mex_destroy(&c3_rhs4);
  sf_mex_destroy(&c3_lhs4);
  sf_mex_destroy(&c3_rhs5);
  sf_mex_destroy(&c3_lhs5);
  sf_mex_destroy(&c3_rhs6);
  sf_mex_destroy(&c3_lhs6);
  sf_mex_destroy(&c3_rhs7);
  sf_mex_destroy(&c3_lhs7);
  sf_mex_destroy(&c3_rhs8);
  sf_mex_destroy(&c3_lhs8);
  sf_mex_destroy(&c3_rhs9);
  sf_mex_destroy(&c3_lhs9);
  sf_mex_destroy(&c3_rhs10);
  sf_mex_destroy(&c3_lhs10);
  sf_mex_destroy(&c3_rhs11);
  sf_mex_destroy(&c3_lhs11);
  sf_mex_destroy(&c3_rhs12);
  sf_mex_destroy(&c3_lhs12);
  sf_mex_destroy(&c3_rhs13);
  sf_mex_destroy(&c3_lhs13);
  sf_mex_destroy(&c3_rhs14);
  sf_mex_destroy(&c3_lhs14);
  sf_mex_destroy(&c3_rhs15);
  sf_mex_destroy(&c3_lhs15);
}

static const mxArray *c3_e_emlrt_marshallOut(const char * c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", c3_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c3_u)), false);
  return c3_y;
}

static const mxArray *c3_f_emlrt_marshallOut(const uint32_T c3_u)
{
  const mxArray *c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 7, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static const mxArray *c3_g_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const int32_T c3_u)
{
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 6, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static const mxArray *c3_c_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  SFc3_hw8_two_tanksInstanceStruct *chartInstance;
  chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  sf_mex_assign(&c3_mxArrayOutData, c3_g_emlrt_marshallOut(chartInstance,
    *(int32_T *)c3_inData), false);
  return c3_mxArrayOutData;
}

static int32_T c3_e_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_b_sfEvent, const char_T *c3_identifier)
{
  int32_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_sfEvent),
    &c3_thisId);
  sf_mex_destroy(&c3_b_sfEvent);
  return c3_y;
}

static int32_T c3_f_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  int32_T c3_y;
  int32_T c3_i0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_i0, 1, 6, 0U, 0, 0U, 0);
  c3_y = c3_i0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  SFc3_hw8_two_tanksInstanceStruct *chartInstance;
  chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)chartInstanceVoid;
  *(int32_T *)c3_outData = c3_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_mxArrayInData), c3_varName);
  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_h_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const uint8_T c3_u)
{
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_create("y", &c3_u, 3, 0U, 0U, 0U, 0), false);
  return c3_y;
}

static const mxArray *c3_d_sf_marshallOut(void *chartInstanceVoid, void
  *c3_inData)
{
  const mxArray *c3_mxArrayOutData = NULL;
  SFc3_hw8_two_tanksInstanceStruct *chartInstance;
  chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)chartInstanceVoid;
  c3_mxArrayOutData = NULL;
  sf_mex_assign(&c3_mxArrayOutData, c3_h_emlrt_marshallOut(chartInstance,
    *(uint8_T *)c3_inData), false);
  return c3_mxArrayOutData;
}

static uint8_T c3_g_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_b_tp_From_tank_1, const char_T
  *c3_identifier)
{
  uint8_T c3_y;
  emlrtMsgIdentifier c3_thisId;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  c3_y = c3_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c3_b_tp_From_tank_1),
    &c3_thisId);
  sf_mex_destroy(&c3_b_tp_From_tank_1);
  return c3_y;
}

static uint8_T c3_h_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  uint8_T c3_y;
  uint8_T c3_u0;
  (void)chartInstance;
  sf_mex_import(c3_parentId, sf_mex_dup(c3_u), &c3_u0, 1, 3, 0U, 0, 0U, 0);
  c3_y = c3_u0;
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c3_mxArrayInData, const char_T *c3_varName, void *c3_outData)
{
  SFc3_hw8_two_tanksInstanceStruct *chartInstance;
  chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)chartInstanceVoid;
  *(uint8_T *)c3_outData = c3_g_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_mxArrayInData), c3_varName);
  sf_mex_destroy(&c3_mxArrayInData);
}

static const mxArray *c3_i_emlrt_marshallOut(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  const mxArray *c3_y;
  real_T *c3_Qa;
  real_T *c3_Qin;
  real_T *c3_Qout;
  real_T *c3_h1_out;
  real_T *c3_h2_out;
  real_T *c3_mode_out;
  real_T *c3_h1;
  real_T *c3_h2;
  c3_Qout = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c3_h2_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c3_Qa = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c3_h1_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c3_Qin = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c3_h2 = (real_T *)(ssGetContStates(chartInstance->S) + 1);
  c3_h1 = (real_T *)(ssGetContStates(chartInstance->S) + 0);
  c3_mode_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c3_y = NULL;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_createcellmatrix(10, 1), false);
  sf_mex_setcell(c3_y, 0, c3_c_emlrt_marshallOut(chartInstance, *c3_Qa));
  sf_mex_setcell(c3_y, 1, c3_c_emlrt_marshallOut(chartInstance, *c3_Qin));
  sf_mex_setcell(c3_y, 2, c3_c_emlrt_marshallOut(chartInstance, *c3_Qout));
  sf_mex_setcell(c3_y, 3, c3_c_emlrt_marshallOut(chartInstance, *c3_h1_out));
  sf_mex_setcell(c3_y, 4, c3_c_emlrt_marshallOut(chartInstance, *c3_h2_out));
  sf_mex_setcell(c3_y, 5, c3_c_emlrt_marshallOut(chartInstance, *c3_mode_out));
  sf_mex_setcell(c3_y, 6, c3_c_emlrt_marshallOut(chartInstance, *c3_h1));
  sf_mex_setcell(c3_y, 7, c3_c_emlrt_marshallOut(chartInstance, *c3_h2));
  sf_mex_setcell(c3_y, 8, c3_h_emlrt_marshallOut(chartInstance,
    chartInstance->c3_is_active_c3_hw8_two_tanks));
  sf_mex_setcell(c3_y, 9, c3_h_emlrt_marshallOut(chartInstance,
    chartInstance->c3_is_c3_hw8_two_tanks));
  return c3_y;
}

static void c3_i_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u)
{
  real_T *c3_Qa;
  real_T *c3_Qin;
  real_T *c3_Qout;
  real_T *c3_h1_out;
  real_T *c3_h2_out;
  real_T *c3_mode_out;
  real_T *c3_h1;
  real_T *c3_h2;
  c3_Qout = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
  c3_h2_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
  c3_Qa = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
  c3_h1_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c3_Qin = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c3_h2 = (real_T *)(ssGetContStates(chartInstance->S) + 1);
  c3_h1 = (real_T *)(ssGetContStates(chartInstance->S) + 0);
  c3_mode_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  *c3_Qa = c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 0)),
    "Qa");
  *c3_Qin = c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 1)),
    "Qin");
  *c3_Qout = c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u,
    2)), "Qout");
  *c3_h1_out = c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u,
    3)), "h1_out");
  *c3_h2_out = c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u,
    4)), "h2_out");
  *c3_mode_out = c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c3_u, 5)), "mode_out");
  *c3_h1 = c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 6)),
    "h1");
  *c3_h2 = c3_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 7)),
    "h2");
  chartInstance->c3_is_active_c3_hw8_two_tanks = c3_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c3_u, 8)),
     "is_active_c3_hw8_two_tanks");
  chartInstance->c3_is_c3_hw8_two_tanks = c3_g_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c3_u, 9)), "is_c3_hw8_two_tanks");
  sf_mex_assign(&chartInstance->c3_setSimStateSideEffectsInfo,
                c3_j_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c3_u, 10)), "setSimStateSideEffectsInfo"), true);
  sf_mex_destroy(&c3_u);
}

static const mxArray *c3_j_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_b_setSimStateSideEffectsInfo, const char_T
  *c3_identifier)
{
  const mxArray *c3_y = NULL;
  emlrtMsgIdentifier c3_thisId;
  c3_y = NULL;
  c3_thisId.fIdentifier = c3_identifier;
  c3_thisId.fParent = NULL;
  sf_mex_assign(&c3_y, c3_k_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c3_b_setSimStateSideEffectsInfo), &c3_thisId), false);
  sf_mex_destroy(&c3_b_setSimStateSideEffectsInfo);
  return c3_y;
}

static const mxArray *c3_k_emlrt_marshallIn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, const mxArray *c3_u, const emlrtMsgIdentifier *c3_parentId)
{
  const mxArray *c3_y = NULL;
  (void)chartInstance;
  (void)c3_parentId;
  c3_y = NULL;
  sf_mex_assign(&c3_y, sf_mex_duplicatearraysafe(&c3_u), false);
  sf_mex_destroy(&c3_u);
  return c3_y;
}

static void c3_updateDataWrittenToVector(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, uint32_T c3_vectorIndex)
{
  (void)chartInstance;
  c3_dataWrittenToVector[(uint32_T)_SFD_EML_ARRAY_BOUNDS_CHECK(0, (int32_T)
    c3_vectorIndex, 0, 11, 1, 0)] = true;
}

static void c3_errorIfDataNotWrittenToFcn(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance, uint32_T c3_vectorIndex, uint32_T c3_dataNumber)
{
  (void)chartInstance;
  _SFD_DATA_READ_BEFORE_WRITE_CHECK(c3_dataNumber, c3_dataWrittenToVector
    [(uint32_T)_SFD_EML_ARRAY_BOUNDS_CHECK(0, (int32_T)c3_vectorIndex, 0, 11, 1,
    0)]);
}

static void c3_b_sign(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                      *c3_x)
{
  c3_b_eml_scalar_sign(chartInstance, c3_x);
}

static void c3_b_eml_scalar_sign(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T *c3_x)
{
  (void)chartInstance;
  *c3_x = muDoubleScalarSign(*c3_x);
}

static void c3_b_sqrt(SFc3_hw8_two_tanksInstanceStruct *chartInstance, real_T
                      *c3_x)
{
  if (*c3_x < 0.0) {
    c3_eml_error(chartInstance);
  }

  c3_b_eml_scalar_sqrt(chartInstance, c3_x);
}

static void c3_b_eml_scalar_sqrt(SFc3_hw8_two_tanksInstanceStruct *chartInstance,
  real_T *c3_x)
{
  (void)chartInstance;
  *c3_x = muDoubleScalarSqrt(*c3_x);
}

static void init_dsm_address_info(SFc3_hw8_two_tanksInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c3_hw8_two_tanks_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3991860406U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(6513286U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(551922536U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2333334286U);
}

mxArray *sf_c3_hw8_two_tanks_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("xvrckM2qr7mZGHzQyGZD3G");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,6,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxData);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c3_hw8_two_tanks_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c3_hw8_two_tanks_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c3_hw8_two_tanks(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x10'type','srcId','name','auxInfo'{{M[1],M[53],T\"Qa\",},{M[1],M[42],T\"Qin\",},{M[1],M[44],T\"Qout\",},{M[1],M[55],T\"h1_out\",},{M[1],M[56],T\"h2_out\",},{M[1],M[83],T\"mode_out\",},{M[5],M[38],T\"h1\",},{M[5],M[40],T\"h2\",},{M[8],M[0],T\"is_active_c3_hw8_two_tanks\",},{M[9],M[0],T\"is_c3_hw8_two_tanks\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 10, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c3_hw8_two_tanks_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc3_hw8_two_tanksInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _hw8_two_tanksMachineNumber_,
           3,
           8,
           17,
           0,
           26,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize ist own list of scripts */
        init_script_number_translation(_hw8_two_tanksMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_hw8_two_tanksMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _hw8_two_tanksMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,2,0,1,"mode_out");
          _SFD_SET_DATA_PROPS(1,0,0,0,"h1");
          _SFD_SET_DATA_PROPS(2,7,0,0,"h");
          _SFD_SET_DATA_PROPS(3,0,0,0,"h2");
          _SFD_SET_DATA_PROPS(4,2,0,1,"Qin");
          _SFD_SET_DATA_PROPS(5,2,0,1,"h1_out");
          _SFD_SET_DATA_PROPS(6,2,0,1,"Qa");
          _SFD_SET_DATA_PROPS(7,2,0,1,"h2_out");
          _SFD_SET_DATA_PROPS(8,7,0,0,"A");
          _SFD_SET_DATA_PROPS(9,2,0,1,"Qout");
          _SFD_SET_DATA_PROPS(10,1,1,0,"x1");
          _SFD_SET_DATA_PROPS(11,1,1,0,"x2");
          _SFD_SET_DATA_PROPS(12,7,0,0,"rho");
          _SFD_SET_DATA_PROPS(13,7,0,0,"g");
          _SFD_SET_DATA_PROPS(14,7,0,0,"Kin");
          _SFD_SET_DATA_PROPS(15,7,0,0,"Ka");
          _SFD_SET_DATA_PROPS(16,7,0,0,"Kout");
          _SFD_SET_DATA_PROPS(17,1,1,0,"flow");
          _SFD_SET_DATA_PROPS(18,9,0,0,"");
          _SFD_SET_DATA_PROPS(19,9,0,0,"");
          _SFD_SET_DATA_PROPS(20,8,0,0,"");
          _SFD_SET_DATA_PROPS(21,8,0,0,"");
          _SFD_SET_DATA_PROPS(22,9,0,0,"");
          _SFD_SET_DATA_PROPS(23,8,0,0,"");
          _SFD_SET_DATA_PROPS(24,8,0,0,"");
          _SFD_SET_DATA_PROPS(25,9,0,0,"");
          _SFD_STATE_INFO(0,0,0);
          _SFD_STATE_INFO(1,0,0);
          _SFD_STATE_INFO(2,0,0);
          _SFD_STATE_INFO(3,0,0);
          _SFD_STATE_INFO(4,0,2);
          _SFD_STATE_INFO(5,0,2);
          _SFD_STATE_INFO(6,0,2);
          _SFD_STATE_INFO(7,0,2);
          _SFD_CH_SUBSTATE_COUNT(4);
          _SFD_CH_SUBSTATE_DECOMP(0);
          _SFD_CH_SUBSTATE_INDEX(0,0);
          _SFD_CH_SUBSTATE_INDEX(1,1);
          _SFD_CH_SUBSTATE_INDEX(2,2);
          _SFD_CH_SUBSTATE_INDEX(3,3);
          _SFD_ST_SUBSTATE_COUNT(0,0);
          _SFD_ST_SUBSTATE_COUNT(1,0);
          _SFD_ST_SUBSTATE_COUNT(2,0);
          _SFD_ST_SUBSTATE_COUNT(3,0);
        }

        _SFD_CV_INIT_CHART(4,1,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        {
          _SFD_CV_INIT_STATE(1,0,0,0,0,0,NULL,NULL);
        }

        {
          _SFD_CV_INIT_STATE(2,0,0,0,0,0,NULL,NULL);
        }

        {
          _SFD_CV_INIT_STATE(3,0,0,0,0,0,NULL,NULL);
        }

        {
          _SFD_CV_INIT_STATE(4,0,0,0,0,0,NULL,NULL);
        }

        {
          _SFD_CV_INIT_STATE(5,0,0,0,0,0,NULL,NULL);
        }

        {
          _SFD_CV_INIT_STATE(6,0,0,0,0,0,NULL,NULL);
        }

        {
          _SFD_CV_INIT_STATE(7,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(16,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(4,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(7,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(8,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(5,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(6,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(1,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(9,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(13,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(10,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(11,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(14,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(2,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(12,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(3,0,NULL,NULL,0,NULL);
        _SFD_CV_INIT_TRANS(15,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(6,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(6,0,"qi",0,-1,34);
        _SFD_CV_INIT_EML(7,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(7,0,"qo",0,-1,181);
        _SFD_CV_INIT_EML(5,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(5,0,"qa_from_tank",0,-1,234);
        _SFD_CV_INIT_EML(4,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(4,0,"qa_balancing",0,-1,101);
        _SFD_CV_INIT_EML(1,1,0,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML(3,1,0,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML(0,1,0,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML(2,1,0,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML(16,0,0,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML(0,0,0,1,0,0,0,0,0,2,1);
        _SFD_CV_INIT_EML_IF(0,0,0,1,17,1,17);

        {
          static int condStart[] = { 1, 11 };

          static int condEnd[] = { 7, 17 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,0,0,1,17,2,0,&(condStart[0]),&(condEnd[0]),3,
                                &(pfixExpr[0]));
        }

        _SFD_CV_INIT_EML(1,0,0,1,0,0,0,0,0,2,1);
        _SFD_CV_INIT_EML_IF(1,0,0,1,17,1,17);

        {
          static int condStart[] = { 1, 11 };

          static int condEnd[] = { 7, 17 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(1,0,0,1,17,2,0,&(condStart[0]),&(condEnd[0]),3,
                                &(pfixExpr[0]));
        }

        _SFD_CV_INIT_EML(2,0,0,1,0,0,0,0,0,2,1);
        _SFD_CV_INIT_EML_IF(2,0,0,1,17,1,17);

        {
          static int condStart[] = { 1, 11 };

          static int condEnd[] = { 7, 17 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(2,0,0,1,17,2,0,&(condStart[0]),&(condEnd[0]),3,
                                &(pfixExpr[0]));
        }

        _SFD_CV_INIT_EML(3,0,0,1,0,0,0,0,0,2,1);
        _SFD_CV_INIT_EML_IF(3,0,0,1,17,1,17);

        {
          static int condStart[] = { 1, 11 };

          static int condEnd[] = { 7, 17 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(3,0,0,1,17,2,0,&(condStart[0]),&(condEnd[0]),3,
                                &(pfixExpr[0]));
        }

        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)c3_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)c3_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)c3_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)c3_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)c3_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)c3_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)c3_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)c3_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(12,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(13,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(14,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(15,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(16,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(17,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c3_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 4294967295;
          _SFD_SET_DATA_COMPILED_PROPS(18,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)NULL,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4294967295;
          _SFD_SET_DATA_COMPILED_PROPS(19,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)NULL,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4294967295;
          _SFD_SET_DATA_COMPILED_PROPS(20,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)NULL,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4294967295;
          _SFD_SET_DATA_COMPILED_PROPS(21,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)NULL,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4294967295;
          _SFD_SET_DATA_COMPILED_PROPS(22,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)NULL,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4294967295;
          _SFD_SET_DATA_COMPILED_PROPS(23,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)NULL,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4294967295;
          _SFD_SET_DATA_COMPILED_PROPS(24,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)NULL,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4294967295;
          _SFD_SET_DATA_COMPILED_PROPS(25,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)NULL,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_VALUE_PTR(18,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(19,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(20,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(21,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(22,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(23,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(24,(void *)(NULL));
        _SFD_SET_DATA_VALUE_PTR(25,(void *)(NULL));

        {
          real_T *c3_mode_out;
          real_T *c3_h1;
          real_T *c3_h2;
          real_T *c3_Qin;
          real_T *c3_h1_out;
          real_T *c3_Qa;
          real_T *c3_h2_out;
          real_T *c3_Qout;
          real_T *c3_x1;
          real_T *c3_x2;
          real_T *c3_flow;
          c3_flow = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c3_x2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c3_x1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          c3_Qout = (real_T *)ssGetOutputPortSignal(chartInstance->S, 6);
          c3_h2_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 5);
          c3_Qa = (real_T *)ssGetOutputPortSignal(chartInstance->S, 4);
          c3_h1_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
          c3_Qin = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c3_h2 = (real_T *)(ssGetContStates(chartInstance->S) + 1);
          c3_h1 = (real_T *)(ssGetContStates(chartInstance->S) + 0);
          c3_mode_out = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          _SFD_SET_DATA_VALUE_PTR(0U, c3_mode_out);
          _SFD_SET_DATA_VALUE_PTR(1U, c3_h1);
          _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c3_h);
          _SFD_SET_DATA_VALUE_PTR(3U, c3_h2);
          _SFD_SET_DATA_VALUE_PTR(4U, c3_Qin);
          _SFD_SET_DATA_VALUE_PTR(5U, c3_h1_out);
          _SFD_SET_DATA_VALUE_PTR(6U, c3_Qa);
          _SFD_SET_DATA_VALUE_PTR(7U, c3_h2_out);
          _SFD_SET_DATA_VALUE_PTR(8U, &chartInstance->c3_A);
          _SFD_SET_DATA_VALUE_PTR(9U, c3_Qout);
          _SFD_SET_DATA_VALUE_PTR(10U, c3_x1);
          _SFD_SET_DATA_VALUE_PTR(11U, c3_x2);
          _SFD_SET_DATA_VALUE_PTR(12U, &chartInstance->c3_rho);
          _SFD_SET_DATA_VALUE_PTR(13U, &chartInstance->c3_g);
          _SFD_SET_DATA_VALUE_PTR(14U, &chartInstance->c3_Kin);
          _SFD_SET_DATA_VALUE_PTR(15U, &chartInstance->c3_Ka);
          _SFD_SET_DATA_VALUE_PTR(16U, &chartInstance->c3_Kout);
          _SFD_SET_DATA_VALUE_PTR(17U, c3_flow);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _hw8_two_tanksMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "ciBPslj5pOvjlkQKUFKZlB";
}

static void sf_opaque_initialize_c3_hw8_two_tanks(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc3_hw8_two_tanksInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*)
    chartInstanceVar);
  initialize_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c3_hw8_two_tanks(void *chartInstanceVar)
{
  enable_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c3_hw8_two_tanks(void *chartInstanceVar)
{
  disable_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_zeroCrossings_c3_hw8_two_tanks(void *chartInstanceVar)
{
  zeroCrossings_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_outputs_c3_hw8_two_tanks(void *chartInstanceVar)
{
  outputs_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_derivatives_c3_hw8_two_tanks(void *chartInstanceVar)
{
  derivatives_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c3_hw8_two_tanks(void *chartInstanceVar)
{
  sf_gateway_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c3_hw8_two_tanks(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c3_hw8_two_tanks
    ((SFc3_hw8_two_tanksInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c3_hw8_two_tanks();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c3_hw8_two_tanks(SimStruct* S, const
  mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[3];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxDuplicateArray(st);      /* high level simctx */
  prhs[2] = (mxArray*) sf_get_sim_state_info_c3_hw8_two_tanks();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 3, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c3_hw8_two_tanks(SimStruct* S)
{
  return sf_internal_get_sim_state_c3_hw8_two_tanks(S);
}

static void sf_opaque_set_sim_state_c3_hw8_two_tanks(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c3_hw8_two_tanks(S, st);
}

static void sf_opaque_terminate_c3_hw8_two_tanks(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc3_hw8_two_tanksInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_hw8_two_tanks_optimization_info();
    }

    finalize_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*)
      chartInstanceVar);
    utFree((void *)chartInstanceVar);
    if (crtInfo != NULL) {
      utFree((void *)crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c3_hw8_two_tanks(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c3_hw8_two_tanks((SFc3_hw8_two_tanksInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c3_hw8_two_tanks(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_hw8_two_tanks_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,3);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,3,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,3,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,3);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,3,3);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,3,6);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=6; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 3; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,3);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3118984429U));
  ssSetChecksum1(S,(2805444328U));
  ssSetChecksum2(S,(1348978089U));
  ssSetChecksum3(S,(2965630929U));
  ssSetNumContStates(S,2);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c3_hw8_two_tanks(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Stateflow");
  }
}

static void mdlStart_c3_hw8_two_tanks(SimStruct *S)
{
  SFc3_hw8_two_tanksInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc3_hw8_two_tanksInstanceStruct *)utMalloc(sizeof
    (SFc3_hw8_two_tanksInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc3_hw8_two_tanksInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 0;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c3_hw8_two_tanks;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c3_hw8_two_tanks;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c3_hw8_two_tanks;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c3_hw8_two_tanks;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c3_hw8_two_tanks;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c3_hw8_two_tanks;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c3_hw8_two_tanks;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c3_hw8_two_tanks;
  chartInstance->chartInfo.zeroCrossings =
    sf_opaque_zeroCrossings_c3_hw8_two_tanks;
  chartInstance->chartInfo.outputs = sf_opaque_outputs_c3_hw8_two_tanks;
  chartInstance->chartInfo.derivatives = sf_opaque_derivatives_c3_hw8_two_tanks;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c3_hw8_two_tanks;
  chartInstance->chartInfo.mdlStart = mdlStart_c3_hw8_two_tanks;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c3_hw8_two_tanks;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c3_hw8_two_tanks_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c3_hw8_two_tanks(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c3_hw8_two_tanks(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c3_hw8_two_tanks(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c3_hw8_two_tanks_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
