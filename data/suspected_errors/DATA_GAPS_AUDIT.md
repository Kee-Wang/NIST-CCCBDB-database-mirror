# Data Gap Audit

*Generated: 2026-02-19 09:35:20*
*Species audited: 2186*

## Executive Summary

| Category | Count | Description |
|----------|------:|-------------|
| **Missing data** | 213 | Data doesn't exist on CCCBDB |
| **Data quality** | 2 | Data exists but is wrong |
| **Total** | 215 | |

### Breakdown by Code

| Code | Count | Meaning |
|------|------:|---------|
| `cas_typo` | 1 | Invalid CAS on species list |
| `no_calc_geometry` | 66 | Geometry page genuinely empty on CCCBDB |
| `no_energy_data` | 27 | Energy page genuinely empty on CCCBDB |
| `no_spin_data` | 27 | Spin page genuinely empty on CCCBDB |
| `single_atom` | 92 | Atoms have no molecular structure |
| `suspect_energy` | 1 | Identical energies for different basis sets |
| `wrong_spin` | 1 | Parity violation or known wrong spin |

## Missing Data

Data genuinely does not exist on CCCBDB, or could not be fetched.

### `cas_typo` (1 issue)

**metadata** (1):

CAS 137019869 appears to be a typo of CAS 13709869 (extra digit '1'). This entry has no energy data and no geometry data. The valid entry is BF3H3N_13709869 with full data.

[`BF3H3N_137019869`](../molecules/BF3H3N_137019869/expgeom2x.html)

### `no_calc_geometry` (66 issues)

**geometry** (66):

Geometry page loaded with molecule name but no coordinates

[`BCl3H3N_49860186`](../molecules/BCl3H3N_49860186/geom3x.html), [`BF3H3N_137019869`](../molecules/BF3H3N_137019869/geom3x.html), [`BrO3-_32062149`](../molecules/BrO3-_32062149/geom3x.html), [`C13H10O_119619`](../molecules/C13H10O_119619/geom3x.html)
[`C18H12_92240`](../molecules/C18H12_92240/geom3x.html), [`C2Cl3F_359295`](../molecules/C2Cl3F_359295/geom3x.html), [`C2H3N3_288880`](../molecules/C2H3N3_288880/geom3x.html), [`C2HCl5_76017`](../molecules/C2HCl5_76017/geom3x.html)
[`C3H2Cl2_83682320`](../molecules/C3H2Cl2_83682320/geom3x.html), [`C3H2O_2961800`](../molecules/C3H2O_2961800/geom3x.html), [`C3H4O3_96491`](../molecules/C3H4O3_96491/geom3x.html), [`C3H5F_1184607`](../molecules/C3H5F_1184607/geom3x.html)
[`C3H6Cl2_142289`](../molecules/C3H6Cl2_142289/geom3x.html), [`C3H6Cl2_78999`](../molecules/C3H6Cl2_78999/geom3x.html), [`C3H6F2_420451`](../molecules/C3H6F2_420451/geom3x.html), [`C3H7+_2143615`](../molecules/C3H7+_2143615/geom3x.html)
[`C3H7N_503297`](../molecules/C3H7N_503297/geom3x.html), [`C3H7O+_43022033`](../molecules/C3H7O+_43022033/geom3x.html), [`C3H9B_593908`](../molecules/C3H9B_593908/geom3x.html), [`C3H9O+_43625656`](../molecules/C3H9O+_43625656/geom3x.html)
[`C4H4O2_497234`](../molecules/C4H4O2_497234/geom3x.html), [`C4H5N_109751`](../molecules/C4H5N_109751/geom3x.html), [`C4H6S_627510`](../molecules/C4H6S_627510/geom3x.html), [`C4H8Cl2_2211678`](../molecules/C4H8Cl2_2211678/geom3x.html)
[`C4H8O2S_126330`](../molecules/C4H8O2S_126330/geom3x.html), [`C4H8S_110010`](../molecules/C4H8S_110010/geom3x.html), [`C4H9F_2366521`](../molecules/C4H9F_2366521/geom3x.html), [`C4H9F_359013`](../molecules/C4H9F_359013/geom3x.html)
[`C4O4--_28737408`](../molecules/C4O4--_28737408/geom3x.html), [`C5H12S_1679090`](../molecules/C5H12S_1679090/geom3x.html), [`C5H5N_2180689`](../molecules/C5H5N_2180689/geom3x.html), [`C5H8_627214`](../molecules/C5H8_627214/geom3x.html)
[`C6H12O_66251`](../molecules/C6H12O_66251/geom3x.html), [`C6H4F2_372189`](../molecules/C6H4F2_372189/geom3x.html), [`C6H5_2396012`](../molecules/C6H5_2396012/geom3x.html), [`C6H6S_108985`](../molecules/C6H6S_108985/geom3x.html)
[`C6H6_5291907`](../molecules/C6H6_5291907/geom3x.html), [`C7H10_26456633`](../molecules/C7H10_26456633/geom3x.html), [`C7H7_3551277`](../molecules/C7H7_3551277/geom3x.html), [`C7H8_278068`](../molecules/C7H8_278068/geom3x.html)
[`C7H8_544252`](../molecules/C7H8_544252/geom3x.html), [`C8H10_106423`](../molecules/C8H10_106423/geom3x.html), [`CBr2F2_75616`](../molecules/CBr2F2_75616/geom3x.html), [`CBr3I_14349805`](../molecules/CBr3I_14349805/geom3x.html)
[`CH2S3_594081`](../molecules/CH2S3_594081/geom3x.html), [`CH3N+_35430172`](../molecules/CH3N+_35430172/geom3x.html), [`CH4NO+_9000924`](../molecules/CH4NO+_9000924/geom3x.html), [`CH7N2+_62901706`](../molecules/CH7N2+_62901706/geom3x.html)
[`CHCl2I_594047`](../molecules/CHCl2I_594047/geom3x.html), [`CHI3_75478`](../molecules/CHI3_75478/geom3x.html), [`CHN2_20813325`](../molecules/CHN2_20813325/geom3x.html), [`CHO2_2564865`](../molecules/CHO2_2564865/geom3x.html)
[`CHS2+_87552861`](../molecules/CHS2+_87552861/geom3x.html), [`Cl3Si_19165345`](../molecules/Cl3Si_19165345/geom3x.html), [`ClH2N_110599903`](../molecules/ClH2N_110599903/geom3x.html), [`ClH3N+_12432483`](../molecules/ClH3N+_12432483/geom3x.html)
[`ClO3_13932100`](../molecules/ClO3_13932100/geom3x.html), [`F7I_16921963`](../molecules/F7I_16921963/geom3x.html), [`FZn+_9000139`](../molecules/FZn+_9000139/geom3x.html), [`H2KO+_9000918`](../molecules/H2KO+_9000918/geom3x.html)
[`H2NO_13408292`](../molecules/H2NO_13408292/geom3x.html), [`HN2O_107190772`](../molecules/HN2O_107190772/geom3x.html), [`HN2_36882130`](../molecules/HN2_36882130/geom3x.html), [`Li2O2_12031800`](../molecules/Li2O2_12031800/geom3x.html)
[`LiO2_120472248`](../molecules/LiO2_120472248/geom3x.html), [`Na2O_1313593`](../molecules/Na2O_1313593/geom3x.html)

### `no_energy_data` (27 issues)

**energy** (27):

Energy page loaded with molecule name but no data

[`BF3H3N_137019869`](../molecules/BF3H3N_137019869/energy2x.html), [`C3H5F_1184607`](../molecules/C3H5F_1184607/energy2x.html), [`C3H6Cl2_78999`](../molecules/C3H6Cl2_78999/energy2x.html), [`C3H6F2_420451`](../molecules/C3H6F2_420451/energy2x.html)
[`C3H7+_2143615`](../molecules/C3H7+_2143615/energy2x.html), [`C4H5N_109751`](../molecules/C4H5N_109751/energy2x.html), [`C4H8O2S_126330`](../molecules/C4H8O2S_126330/energy2x.html), [`C4H9F_359013`](../molecules/C4H9F_359013/energy2x.html)
[`C6H12_18931710`](../molecules/C6H12_18931710/energy2x.html), [`C6H6_5291907`](../molecules/C6H6_5291907/energy2x.html), [`CBr3I_14349805`](../molecules/CBr3I_14349805/energy2x.html), [`CH4NO+_9000924`](../molecules/CH4NO+_9000924/energy2x.html)
[`CHCl2I_594047`](../molecules/CHCl2I_594047/energy2x.html), [`CHI3_75478`](../molecules/CHI3_75478/energy2x.html), [`CHS2+_87552861`](../molecules/CHS2+_87552861/energy2x.html), [`Cl5Sb_7647189`](../molecules/Cl5Sb_7647189/energy2x.html)
[`ClH2N_110599903`](../molecules/ClH2N_110599903/energy2x.html), [`ClSb_19952126`](../molecules/ClSb_19952126/energy2x.html), [`Co_7440484`](../molecules/Co_7440484/energy2x.html), [`Cr_7440473`](../molecules/Cr_7440473/energy2x.html)
[`F7I_16921963`](../molecules/F7I_16921963/energy2x.html), [`FSb_25285727`](../molecules/FSb_25285727/energy2x.html), [`H2KO+_9000918`](../molecules/H2KO+_9000918/energy2x.html), [`IK_7681110`](../molecules/IK_7681110/energy2x.html)
[`LiO2_120472248`](../molecules/LiO2_120472248/energy2x.html), [`Mn_7439965`](../molecules/Mn_7439965/energy2x.html), [`NSb_12333572`](../molecules/NSb_12333572/energy2x.html)

### `no_spin_data` (27 issues)

**spin** (27):

Spin page loaded but no spin data found (genuinely empty)

[`C2F3N+_353855`](../molecules/C2F3N+_353855/spin2x.html), [`C2F4+_116143`](../molecules/C2F4+_116143/spin2x.html), [`C2H2O2+_107222`](../molecules/C2H2O2+_107222/spin2x.html), [`C3H3_9000220`](../molecules/C3H3_9000220/spin2x.html)
[`C4H10O+_71363`](../molecules/C4H10O+_71363/spin2x.html), [`C4H4O+_110009`](../molecules/C4H4O+_110009/spin2x.html), [`C6H6+_71432`](../molecules/C6H6+_71432/spin2x.html), [`CCl4+_56235`](../molecules/CCl4+_56235/spin2x.html)
[`CH3Cl+_74873`](../molecules/CH3Cl+_74873/spin2x.html), [`CH4NO+_50785803`](../molecules/CH4NO+_50785803/spin2x.html), [`CH4NO+_9000924`](../molecules/CH4NO+_9000924/spin2x.html), [`CHN2_20813325`](../molecules/CHN2_20813325/spin2x.html)
[`CaH+_14452756`](../molecules/CaH+_14452756/spin2x.html), [`ClO3_13932100`](../molecules/ClO3_13932100/spin2x.html), [`ClSb_19952126`](../molecules/ClSb_19952126/spin2x.html), [`Co_7440484`](../molecules/Co_7440484/spin2x.html)
[`Cr_7440473`](../molecules/Cr_7440473/spin2x.html), [`F2H_18130740`](../molecules/F2H_18130740/spin2x.html), [`FSb_25285727`](../molecules/FSb_25285727/spin2x.html), [`HN2O_107190772`](../molecules/HN2O_107190772/spin2x.html)
[`LiO2_120472248`](../molecules/LiO2_120472248/spin2x.html), [`Mn_7439965`](../molecules/Mn_7439965/spin2x.html), [`NTi_25583204`](../molecules/NTi_25583204/spin2x.html), [`Ni_7440020`](../molecules/Ni_7440020/spin2x.html)
[`OTe_13451177`](../molecules/OTe_13451177/spin2x.html), [`Ti+_7440326`](../molecules/Ti+_7440326/spin2x.html), [`Ti-_7440326`](../molecules/Ti-_7440326/spin2x.html)

### `single_atom` (92 issues)

**geometry** (92):

Single atom — no molecular geometry possible

`Al+_7429905`, `Al-_7429905`, `Al_7429905`, `As+_7440382`, `As-_7440382`, `As_7440382`, `B+_7440428`, `B-_7440428`, `B_7440428`, `Be+_7440417`
`Be-_7440417`, `Be_7440417`, `C+_7440440`, `C-_7440440`, `C_7440440`, `Ca+_7440702`, `Ca-_7440702`, `Ca_7440702`, `Cl+_22537151`, `Cl-_22537151`
`Cl_22537151`, `Co_7440484`, `Cr_7440473`, `D+_16873179`, `D-_16873179`, `D_16873179`, `F+_14762948`, `F-_14762948`, `F_14762948`, `Fe_7439896`
`Ga+_7440553`, `Ga-_7440553`, `Ga_7440553`, `Ge+_7440564`, `Ge-_7440564`, `Ge_7440564`, `H+_12385136`, `H-_12385136`, `H_12385136`, `I+_14362448`
`I-_14362448`, `I_14362448`, `K+_7440097`, `K-_7440097`, `K_7440097`, `Kr+_7439909`, `Kr_7439909`, `Li+_7439932`, `Li-_7439932`, `Li_7439932`
`Mg+_7439954`, `Mg-_7439954`, `Mg_7439954`, `Mn_7439965`, `N+_17778880`, `N-_17778880`, `N_17778880`, `Na+_7440235`, `Na-_7440235`, `Na_7440235`
`Ni_7440020`, `O+_17778802`, `O-_17778802`, `O_17778802`, `P+_7723140`, `P-_7723140`, `P_7723140`, `S+_7704349`, `S-_7704349`, `S_7704349`
`Sb_7440360`, `Sc+_7440202`, `Sc_7440202`, `Se+_7782492`, `Se-_7782492`, `Se_7782492`, `Si+_7440213`, `Si-_7440213`, `Si_7440213`, `Sn+_7440315`
`Sn_7440315`, `Te+_13494809`, `Te_13494809`, `Ti+_7440326`, `Ti-_7440326`, `Ti_7440326`, `V+_7440622`, `V_7440622`, `Xe_7440633`, `Zn+_7440666`
`Zn-_7440666`, `Zn_7440666`

## Data Quality

Data exists on CCCBDB but is known to be wrong.

### `suspect_energy` (1 issue)

**energy** (1):

CCSD(T)=FULL/STO-3G and CCSD(T)=FULL/6-31G* report identical energy (-39.806897 hartree). Different basis sets should produce different energies — likely a data copy error in CCCBDB.

[`CH4_74828`](../molecules/CH4_74828/energy2x.html)

### `wrong_spin` (1 issue)

**spin** (1):

Fe reported as closed shell (multiplicity=1, S^2=0). Atomic iron ground state is 5D4 (multiplicity=5, four unpaired 3d electrons). Archived spin2x.html confirms CCCBDB genuinely reports 'Fe (Iron atom) is closed shell. S^2=0' — upstream data error.

[`Fe_7439896`](../molecules/Fe_7439896/spin2x.html)

## Recovery

```bash
# Heal session failures (scan, clear bad HTML, re-fetch)
python scripts/cccbdb_mirror.py --heal

# Re-fetch all missing pages
python scripts/cccbdb_mirror.py --pages all --missing-only

# Re-run audit to track progress
python scripts/audit_gaps.py --diff
```
