#!/usr/bin/env python3
from pathlib import Path

root=Path('HICLASS')
header=root/'include/background.h'
model=root/'gravity_smg/gravity_models_smg.c'

# Add covariant model enum.
s=header.read_text()
old='''enum gravity_model {propto_omega, propto_scale,\n    constant_alphas,\n    eft_alphas_power_law, eft_gammas_power_law, eft_gammas_exponential,\n    galileon, nkgb,\n    brans_dicke,\n    quintessence_monomial, quintessence_tracker,\n    alpha_attractor_canonical\n};'''
new='''enum gravity_model {propto_omega, propto_scale,\n    constant_alphas,\n    eft_alphas_power_law, eft_gammas_power_law, eft_gammas_exponential,\n    galileon, nkgb,\n    brans_dicke, epjc_fr_covariant,\n    quintessence_monomial, quintessence_tracker,\n    alpha_attractor_canonical\n};'''
assert old in s
header.write_text(s.replace(old,new,1))

s=model.read_text()
# Register model immediately after Brans-Dicke properties block.
needle='''  if (strcmp(string1,"nkgb") == 0 || strcmp(string1,"n-kgb") == 0 || strcmp(string1,"N-KGB") == 0 || strcmp(string1,"nKGB") == 0) {'''
insert='''  if (strcmp(string1,"epjc_fr_covariant") == 0) {\n    /* Exact covariant scalar-tensor representation of the pole-free\n     * curvature completion. parameters_smg =\n     * [Lambda/H0^2, alpha*H0^2, R_s/H0^2, R_reg/H0^2].\n     * Lambda is the closure-shooting parameter. */\n    pba->gravity_model_smg = epjc_fr_covariant;\n    pba->field_evolution_smg = _TRUE_;\n    flag2=_TRUE_;\n    pba->parameters_size_smg = 4;\n    class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);\n    if (has_tuning_index_smg == _FALSE_) pba->tuning_index_smg = 0;\n    if (has_dxdy_guess_smg == _FALSE_) pba->tuning_dxdy_guess_smg = 3.0;\n    class_test(pba->parameters_smg[1] <= 0., errmsg, "epjc_fr_covariant requires alpha_hat>0");\n    class_test(pba->parameters_smg[2] <= 0., errmsg, "epjc_fr_covariant requires R_s/H0^2>0");\n    class_test(pba->parameters_smg[3] <= pba->parameters_smg[2], errmsg, "epjc_fr_covariant requires R_reg>R_s");\n  }\n\n'''+needle
assert needle in s
s=s.replace(needle,insert,1)

# Add exact Horndeski functions after Brans-Dicke Gs block.
needle2='''  else if(pba->gravity_model_smg == nkgb){'''
insert2='''  else if(pba->gravity_model_smg == epjc_fr_covariant){\n    /* Legendre map for\n     * f=A[R-2 Lambda + 2 alpha R_s^2(sqrt(1+(R/R_s)^2)-1)].\n     * phi=F, G4=phi/2, G2=-U(phi)/2.\n     * The high-curvature inversion is regularized only above R_reg;\n     * convergence in R_reg is tested at spectrum level. */\n    const double lam = pba->parameters_smg[0];\n    const double ah  = pba->parameters_smg[1];\n    const double rs  = pba->parameters_smg[2];\n    const double rreg= pba->parameters_smg[3];\n    const double A   = 1.0/(1.0+2.0*ah*rs);\n    const double b   = 2.0*ah*rs;\n    const double ureg= rreg/sqrt(rreg*rreg+rs*rs);\n    double u=(phi/A-1.0)/b;\n    /* Tiny excursions above the asymptotic branch are mapped to the\n     * chosen high-curvature regulator. Do not disable stability tests. */\n    if (u > ureg) u=ureg;\n    if (u < -ureg) u=-ureg;\n    const double omu2=fmax(1.0-u*u,1.e-30);\n    const double sq=sqrt(omu2);\n    const double V = pow(pba->H0,2)*A*(lam + ah*rs*rs*(1.0-sq));\n    const double Rhalf = 0.5*pow(pba->H0,2)*rs*u/sq;\n    const double Vpp = pow(pba->H0,2)/(4.0*A*ah)*pow(omu2,-1.5);\n\n    pgf->G2 = -V;\n    pgf->G2_phi = -Rhalf;\n    pgf->G2_phiphi = -Vpp;\n    pgf->DG4 = (phi-1.0)/2.0;\n    pgf->G4 = phi/2.0;\n    pgf->G4_phi = 1.0/2.0;\n  }\n\n'''+needle2
assert needle2 in s
s=s.replace(needle2,insert2,1)

# Initial conditions: heavy scalaron on the regulated high-curvature minimum.
needle3='''    case brans_dicke:\n\t\t\tpvecback_integration[pba->index_bi_phi_smg] = pba->parameters_smg[2];\n\t\t\tpvecback_integration[pba->index_bi_phi_prime_smg] = pba->parameters_smg[3];\n\t\t\tbreak;\n'''
insert3=needle3+'''\n    case epjc_fr_covariant:\n      {\n        const double ah=pba->parameters_smg[1];\n        const double rs=pba->parameters_smg[2];\n        const double rreg=pba->parameters_smg[3];\n        const double A=1.0/(1.0+2.0*ah*rs);\n        const double ureg=rreg/sqrt(rreg*rreg+rs*rs);\n        pvecback_integration[pba->index_bi_phi_smg] = A*(1.0+2.0*ah*rs*ureg);\n        pvecback_integration[pba->index_bi_phi_prime_smg] = 0.0;\n      }\n      break;\n'''
assert needle3 in s
s=s.replace(needle3,insert3,1)

# Extend discoverability error string.
s=s.replace("'brans_dicke', 'galileon'", "'brans_dicke', 'epjc_fr_covariant', 'galileon'")
model.write_text(s)
print('patched hi_class with exact covariant EPJC f(R) completion')
