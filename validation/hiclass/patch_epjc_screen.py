#!/usr/bin/env python3
from pathlib import Path

root=Path('HICLASS')
header=root/'include/background.h'
model=root/'gravity_smg/gravity_models_smg.c'

s=header.read_text()
old='''enum gravity_model {propto_omega, propto_scale,\n    constant_alphas,\n    eft_alphas_power_law, eft_gammas_power_law, eft_gammas_exponential,\n    galileon, nkgb,\n    brans_dicke,\n    quintessence_monomial, quintessence_tracker,\n    alpha_attractor_canonical\n};'''
new='''enum gravity_model {propto_omega, propto_scale,\n    constant_alphas,\n    epjc_curvature_screen,\n    eft_alphas_power_law, eft_gammas_power_law, eft_gammas_exponential,\n    galileon, nkgb,\n    brans_dicke,\n    quintessence_monomial, quintessence_tracker,\n    alpha_attractor_canonical\n};'''
assert old in s
header.write_text(s.replace(old,new))

s=model.read_text()
anchor='''  if (strcmp(string1,"constant_alphas") == 0) {\n     pba->gravity_model_smg = constant_alphas;\n     pba->field_evolution_smg = _FALSE_;\n     pba->M2_evolution_smg = _TRUE_;\n     flag2=_TRUE_;\n     pba->parameters_2_size_smg = 5;\n     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);\n   }\n'''
insert=anchor+'''\n  /* EPJ C curvature-screened f(R,T_chi) designer embedding.\n   * parameters_smg = alpha_hat, xi_tilde, R_s/H0^2.\n   * This is the fully dynamical Horndeski metric/scalaron sector; the\n   * selective WKB dark-sector exchange is audited separately. */\n  if (strcmp(string1,"epjc_curvature_screen") == 0) {\n     pba->gravity_model_smg = epjc_curvature_screen;\n     pba->field_evolution_smg = _FALSE_;\n     pba->M2_evolution_smg = _FALSE_;\n     flag2=_TRUE_;\n     pba->parameters_2_size_smg = 3;\n     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);\n   }\n'''
assert anchor in s
s=s.replace(anchor,insert,1)

anchor2='''  else if (pba->gravity_model_smg == constant_alphas) {\n\n    double c_k = pba->parameters_2_smg[0];\n    double c_b = pba->parameters_2_smg[1];\n    double c_m = pba->parameters_2_smg[2];\n    double c_t = pba->parameters_2_smg[3];\n\n    pvecback[pba->index_bg_kineticity_smg] = c_k;\n    pvecback[pba->index_bg_braiding_smg] = c_b;\n    pvecback[pba->index_bg_tensor_excess_smg] = c_t;\n    pvecback[pba->index_bg_M2_running_smg] = c_m;\n    pvecback[pba->index_bg_delta_M2_smg] = delta_M2; //M2-1\n    pvecback[pba->index_bg_M2_smg] = 1.+delta_M2;\n  }\n'''
insert2=anchor2+'''\n  else if (pba->gravity_model_smg == epjc_curvature_screen) {\n    /* Dimensionless definitions: alpha_hat=alpha H0^2,\n     * xi_tilde=xi rho_crit,0, rs=R_s/H0^2. For the validation\n     * baseline all neutrinos are ultra-relativistic, so the trace is\n     * carried by baryons+CDM+Lambda and r_N is exact. */\n    const double ah = pba->parameters_2_smg[0];\n    const double xt = pba->parameters_2_smg[1];\n    const double rs = pba->parameters_2_smg[2];\n    const double H02 = pba->H0*pba->H0;\n    const double rho_m = pvecback[pba->index_bg_rho_b] + pvecback[pba->index_bg_rho_cdm];\n    const double rho_c = pvecback[pba->index_bg_rho_cdm];\n    const double rho_L = pvecback[pba->index_bg_rho_smg];\n    const double r = 3.0*(rho_m + 4.0*rho_L)/H02;\n    const double rN = -9.0*rho_m/H02;\n    const double t = -rho_c/H02;\n    const double tN = -3.0*t;\n    const double x = r/rs;\n    const double x2 = x*x;\n    const double om = 1.0 + x2;\n    const double g1 = (1.0-3.0*x2)/(om*om*om);\n    const double g2h = 12.0*x*(x2-1.0)/(rs*om*om*om*om);\n    const double Ft = 1.0 + 2.0*ah*r/sqrt(om) + xt*t*g1;\n    const double FtN = 2.0*ah*rN/pow(om,1.5) + xt*(tN*g1 + t*g2h*rN);\n    const double norm = 1.0/(1.0+2.0*ah*rs);\n    const double M2 = norm*Ft;\n    const double aM = FtN/Ft;\n\n    pvecback[pba->index_bg_kineticity_smg] = 0.0;\n    pvecback[pba->index_bg_braiding_smg] = -aM;\n    pvecback[pba->index_bg_tensor_excess_smg] = 0.0;\n    pvecback[pba->index_bg_M2_running_smg] = aM;\n    pvecback[pba->index_bg_delta_M2_smg] = M2-1.0;\n    pvecback[pba->index_bg_M2_smg] = M2;\n  }\n'''
assert anchor2 in s
s=s.replace(anchor2,insert2,1)

# Add a no-op initial-condition switch case (M2 is algebraic, not integrated).
needle='''\t  case constant_alphas:\n\t\t\tpvecback_integration[pba->index_bi_delta_M2_smg] = pba->parameters_2_smg[4]-1.;\n\t\t\tbreak;\n'''
repl=needle+'''\n\t  case epjc_curvature_screen:\n\t\t\t/* M2 is supplied algebraically in gravity_models_get_alphas_par_smg. */\n\t\t\tbreak;\n'''
assert needle in s
s=s.replace(needle,repl,1)

# Add stdout diagnostics.
needle2='''    case constant_alphas:\n      printf("Modified gravity: constant_alphas with parameters: \\n");\n      printf(" -> c_K = %g, c_B = %g, c_M = %g, c_T = %g, M_*^2_init = %g \\n",\n\t      pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],\n\t      pba->parameters_2_smg[4]);\n    break;\n'''
repl2=needle2+'''\n    case epjc_curvature_screen:\n      printf("Modified gravity: EPJC curvature-screened designer embedding \\n");\n      printf(" -> alpha_hat = %g, xi_tilde = %g, R_s/H0^2 = %g \\n",\n             pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2]);\n    break;\n'''
assert needle2 in s
s=s.replace(needle2,repl2,1)

# Extend error message only for discoverability (not physics).
s=s.replace("'propto_omega', 'propto_scale', 'constant_alphas',", "'propto_omega', 'propto_scale', 'constant_alphas', 'epjc_curvature_screen',")
model.write_text(s)
print('patched hi_class with epjc_curvature_screen')
