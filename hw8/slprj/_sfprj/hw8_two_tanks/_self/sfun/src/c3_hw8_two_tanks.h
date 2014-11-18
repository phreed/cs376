#ifndef __c3_hw8_two_tanks_h__
#define __c3_hw8_two_tanks_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc3_hw8_two_tanksInstanceStruct
#define typedef_SFc3_hw8_two_tanksInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c3_sfEvent;
  uint8_T c3_tp_From_tank_1;
  uint8_T c3_tp_Separated;
  uint8_T c3_tp_Balancing;
  uint8_T c3_tp_From_tank_2;
  boolean_T c3_isStable;
  boolean_T c3_stateChanged;
  real_T c3_lastMajorTime;
  uint8_T c3_is_active_c3_hw8_two_tanks;
  uint8_T c3_is_c3_hw8_two_tanks;
  real_T c3_h;
  real_T c3_A;
  real_T c3_rho;
  real_T c3_g;
  real_T c3_Kin;
  real_T c3_Ka;
  real_T c3_Kout;
  uint8_T c3_doSetSimStateSideEffects;
  const mxArray *c3_setSimStateSideEffectsInfo;
} SFc3_hw8_two_tanksInstanceStruct;

#endif                                 /*typedef_SFc3_hw8_two_tanksInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c3_hw8_two_tanks_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c3_hw8_two_tanks_get_check_sum(mxArray *plhs[]);
extern void c3_hw8_two_tanks_method_dispatcher(SimStruct *S, int_T method, void *
  data);

#endif
