# File Categorization for MATLAB to Python Transpilation

## Overview

Total MATLAB files: 83
- To DELETE: 15 files (dead/debug code)
- To CONSOLIDATE: 25 files (duplicates/variants)
- To KEEP/REFACTOR: 43 files

---

## Files to DELETE (15 files)

### Debug/Test Scripts
| File | Reason |
|------|--------|
| `sub_task_3/BUGS.m` | Debug script, no function definition |
| `sub_task_3/KING_TESTER.m` | Test script |
| `sub_task_3/tctest.m` | Test script |
| `sub_task_3/tctest2.m` | Test script |
| `sub_task_3/SPEC_TEST.m` | Test script |
| `sub_task_3/TestImage.m` | Test script |

### Unnamed/Temporary Scripts
| File | Reason |
|------|--------|
| `sub_task_3/Untitled.m` | Unnamed temporary script |
| `sub_task_3/Untitled4.m` | Unnamed temporary script |
| `sub_task_3/Untitled13.m` | Unnamed temporary script |
| `sub_task_3/temp_box.m` | Temporary script |
| `sub_task_3/uptown_funk.m` | Non-descriptive name, likely test |
| `sub_task_3/using_hg2.m` | Graphics test script |

### External Utilities (Replace with Python equivalents)
| File | Reason |
|------|--------|
| `sub_task_3/export_fig.m` | External utility, use matplotlib instead |

### Superseded Files
| File | Reason |
|------|--------|
| `sub_task_3/FixingMBN.m` | One-off fix script |
| `sub_task_3/RAD_COMP.m` | Obsolete comparison script |

---

## Files to CONSOLIDATE (25 files → 8 modules)

### 1. Myoglobin Readers → `myoglobin_reader.py`
| File | Notes |
|------|-------|
| `sub_task_1/MyoglobinRead.m` | Base version |
| `sub_task_1/MyoglobinRead_1.m` | Variant 1 |
| `sub_task_1/MyoglobinRead_2.m` | Variant 2 |
| `sub_task_1/MyoglobinRead_Revised.m` | Latest version - use as base |

**Python approach**: Single parameterized function with options

### 2. Mode Mapping → `mode_mapper.py`
| File | Notes |
|------|-------|
| `sub_task_3/Make_map.m` | Generic version |
| `sub_task_3/Make_map_MBO.m` | MBO-specific |
| `sub_task_3/Make_map_MBN.m` | MBN-specific |

**Python approach**: Single function with `molecule_type` parameter

### 3. 2D Spectral Plots → `spectral_plots.py`
| File | Notes |
|------|-------|
| `sub_task_3/SD2DPlots.m` | Base version |
| `sub_task_3/SD2DPlots_mbd.m` | MBD variant |
| `sub_task_3/SD2DPlots_Alpha.m` | Alpha variant |
| `sub_task_3/SD2DPlots_Alpha_mbd.m` | Alpha + MBD variant |
| `sub_task_3/SDs_ANGS.m` | Angle variant |

**Python approach**: Single class with configuration options

### 4. Trajectory Generation → `trajectory_generator.py`
| File | Notes |
|------|-------|
| `sub_task_3/Traj_Gen.m` | Base version |
| `sub_task_3/Tran_Gen_2.m` | Iteration 2 |
| `sub_task_3/Tran_Gen_3.m` | Iteration 3 |
| `sub_task_3/Tran_Gen_4.m` | Iteration 4 |
| `sub_task_3/Tran_Gen_5.m` | Iteration 5 - use as base |
| `sub_task_3/Model_Traj_Gen.m` | Model-based version |
| `sub_task_3/Traj_Phase_Adjust.m` | Phase adjustment |

**Python approach**: Single module with multiple functions

### 5. Normal Mode Analysis → `normal_mode_analysis.py`
| File | Notes |
|------|-------|
| `sub_task_3/Gaussian_NM_Dot.m` | Generic version |
| `sub_task_3/Gaussian_NM_Dot_MBO.m` | MBO-specific |
| `sub_task_3/Gaussian_NM_Dot_MBN.m` | MBN-specific |
| `sub_task_3/Gaussian_Extract_NM_DotProd.m` | Extraction variant |

**Python approach**: Single function with molecule parameter

### 6. Gaussian PDB Extract → `gaussian_pdb_extract.py`
| File | Notes |
|------|-------|
| `sub_task_3/Gaussian_Gromacs_PDB_extract.m` | Generic |
| `sub_task_3/Gaussian_Gromacs_PDB_extract_MBO.m` | MBO-specific |
| `sub_task_3/Gaussian_Gromacs_PDB_extract_MBN.m` | MBN-specific |

**Python approach**: Single function with molecule parameter

### 7. Big Data Processing → `batch_processor.py`
| File | Notes |
|------|-------|
| `sub_task_3/BIG_DATA.m` | Original version |
| `sub_task_3/BIG_DATA_BETTER.m` | Improved version - use as base |

**Python approach**: Single module with latest implementation

### 8. Histogram Plotting → `histogram_plots.py`
| File | Notes |
|------|-------|
| `sub_task_3/mbc_hist.m` | MBC histogram |
| `sub_task_3/mbd_hist.m` | MBD histogram |

**Python approach**: Single parameterized function

---

## Files to KEEP and REFACTOR (43 files)

### Core Utilities (8 files) → `tdrr/core/`
| File | Python Target | Priority |
|------|---------------|----------|
| `sub_task_3/apodfun.m` | `spectral.py::apodfun()` | HIGH |
| `sub_task_3/lagwind.m` | `spectral.py::lag_window()` | HIGH |
| `sub_task_3/bisp3cum.m` | `spectral.py::bispectrum()` | HIGH |
| `sub_task_3/lowfilter.m` | `filters.py::lowpass_filter()` | HIGH |
| `sub_task_3/toep.m` | `filters.py::toeplitz_matrix()` | HIGH |
| `sub_task_3/wign.m` | `filters.py::wigner_distribution()` | MEDIUM |
| `sub_task_3/specfitexpn.m` | `fitting.py::spectral_fit_exp()` | MEDIUM |
| `sub_task_3/expnfitt.m` | `fitting.py::exponential_fit()` | MEDIUM |

### Data I/O (8 files) → `tdrr/io/`
| File | Python Target | Priority |
|------|---------------|----------|
| `sub_task_1/GauData.m` | `gaussian.py::GaussianData` | HIGH |
| `sub_task_1/rdEngGau.m` | `gaussian.py::read_energy()` | HIGH |
| `sub_task_1/combData.m` | `gaussian.py::combine_data()` | HIGH |
| `sub_task_1/NormalDataScript.m` | `gaussian.py::normalize_data()` | MEDIUM |
| `sub_task_3/Gaussian_Extract.m` | `gaussian.py::extract_log()` | HIGH |
| `sub_task_3/Gaussian_Write.m` | `gaussian.py::write_input()` | HIGH |
| `sub_task_3/pdb2mat.m` | `pdb.py::parse_pdb()` | HIGH |
| `sub_task_2/sub_task_2_2/XVG_Extract.m` | `gromacs.py::extract_xvg()` | MEDIUM |

### Force Analysis (8 files) → `tdrr/analysis/force.py`
| File | Python Target | Priority |
|------|---------------|----------|
| `sub_task_2/sub_task_2_1/FORCE_PLOTTER.m` | `force_plots.py::plot_forces()` | MEDIUM |
| `sub_task_2/sub_task_2_2/FORCE_PLOTTER.m` | `force_plots.py::plot_forces_v2()` | MEDIUM |
| `sub_task_2/sub_task_2_2/FORCE_VECS.m` | `force.py::compute_vectors()` | MEDIUM |
| `sub_task_2/sub_task_2_2/FORCE_DATA_WRITER.m` | `force.py::write_data()` | LOW |
| `sub_task_2/sub_task_2_2/Energy_Plot.m` | `energy.py::plot_energy()` | MEDIUM |
| `sub_task_2/sub_task_2_2/B_FACTOR_WRITER.m` | `pdb.py::write_b_factors()` | LOW |
| `sub_task_2/sub_task_2_2/GRO_MAKER.m` | `gromacs.py::write_gro()` | LOW |
| `sub_task_2/sub_task_2_2/RES_RANGE_GETTER.m` | `force.py::get_residue_range()` | LOW |

### Gaussian Processing (7 files) → `tdrr/analysis/gaussian.py`
| File | Python Target | Priority |
|------|---------------|----------|
| `sub_task_3/Gaussian_Generate.m` | `gaussian.py::generate_input()` | MEDIUM |
| `sub_task_3/Gaussian_Generate_Grid.m` | `gaussian.py::generate_grid()` | MEDIUM |
| `sub_task_3/Gaussian_Extract_FrcCnstMat.m` | `gaussian.py::extract_force_matrix()` | MEDIUM |
| `sub_task_3/Gaussian_LOG_Tamper.m` | `gaussian.py::modify_log()` | LOW |
| `sub_task_3/Gaussian_LOG_adjuster.m` | `gaussian.py::adjust_log()` | LOW |
| `sub_task_3/MBC_GAU_GJF.m` | `gaussian.py::create_gjf()` | LOW |
| `sub_task_3/MAT2FORTPDB.m` | `converters.py::mat_to_fortran_pdb()` | LOW |

### Spectroscopy Analysis (8 files) → `tdrr/analysis/`
| File | Python Target | Priority |
|------|---------------|----------|
| `sub_task_3/PROCESSOR.m` | `processor.py::run_pipeline()` | HIGH |
| `sub_task_3/COMP_SD2D.m` | `spectral.py::compute_sd2d()` | HIGH |
| `sub_task_3/Kcorr.m` | `correlation.py::compute_kcorr()` | HIGH |
| `sub_task_3/KCORRARR.m` | `correlation.py::compute_kcorr_array()` | MEDIUM |
| `sub_task_3/Tri_Corr_start.m` | `correlation.py::triple_correlation()` | MEDIUM |
| `sub_task_3/DINH.m` | `dihedral.py::compute_dihedral()` | MEDIUM |
| `sub_task_3/GROMACS_RD_WR_ANGS.m` | `gromacs.py::read_write_angles()` | LOW |
| `sub_task_3/MAT2TXT.m` | `converters.py::mat_to_txt()` | LOW |

### Visualization (4 files) → `tdrr/plotting/`
| File | Python Target | Priority |
|------|---------------|----------|
| `sub_task_3/Fig6Plotter.m` | `figures.py::plot_figure_6()` | LOW |
| `sub_task_3/Fig6Writer.m` | `figures.py::write_figure_6()` | LOW |
| `sub_task_3/Solu_Solv_myo.m` | `myoglobin_plots.py::plot_solvation()` | LOW |
| `sub_task_3/WAT_READER.m` | `water.py::read_water_data()` | LOW |

---

## Dependencies Map

### Core Dependencies (must transpile first)
```
apodfun.m ← used by many spectral analysis files
lagwind.m ← used by bisp3cum.m
toep.m ← used by bisp3cum.m
```

### I/O Dependencies
```
GauData.m ← data structure used throughout
pdb2mat.m ← used by force analysis
rdEngGau.m ← used by data extraction scripts
```

### Analysis Dependencies
```
bisp3cum.m ← uses lagwind.m, toep.m
PROCESSOR.m ← uses apodfun.m
SD2DPlots.m ← uses apodfun.m
Make_map.m ← uses Gaussian_Extract.m
```

---

## Refactoring Notes

### Remove Interactive Input
Files with interactive `input()` calls to refactor:
- `lagwind.m` - lines 44-97
- `bisp3cum.m` - lines 80-234

**Action**: Remove all `input()` calls, require parameters

### Extract Hardcoded Paths
Files with hardcoded Windows paths:
- `Gaussian_Extract.m` - lines 2-10
- All Gaussian processing files

**Action**: Extract to configuration or function parameters

### Split I/O from Computation
Files mixing I/O with computation:
- `bisp3cum.m` - has plotting code at end (lines 340-394)
- `Gaussian_Extract.m` - full script with I/O

**Action**: Separate into pure computation and I/O wrapper functions

### Remove Debug Blocks
Files with `if 0` or `if false` blocks:
- `Gaussian_Extract.m` - lines 130-140
- `pdb2mat.m` - commented out sections

**Action**: Remove all dead code blocks

---

## Priority Order for Transpilation

### Phase 1: Core (must do first)
1. `apodfun.m`
2. `lagwind.m`
3. `toep.m`
4. `bisp3cum.m`

### Phase 2: Data Structures
1. `GauData.m`
2. `pdb2mat.m`
3. `rdEngGau.m`

### Phase 3: I/O
1. `Gaussian_Extract.m`
2. `Gaussian_Write.m`
3. `XVG_Extract.m`
4. `GRO_MAKER.m`

### Phase 4: Analysis
1. `PROCESSOR.m` (refactor to function first)
2. `Make_map.m` (consolidated)
3. `SD2DPlots.m` (consolidated)
4. `Traj_Gen.m` (consolidated)

### Phase 5: Secondary Analysis
1. `Kcorr.m`
2. `COMP_SD2D.m`
3. Force analysis files

### Phase 6: Visualization
1. All plotting functions

---

## Success Metrics

- [ ] All 15 dead files deleted
- [ ] 25 files consolidated into 8 modules
- [ ] 43 files refactored with proper docstrings
- [ ] All hardcoded paths extracted
- [ ] All interactive input removed
- [ ] All functions are pure (no side effects)
- [ ] 100% test coverage for core utilities
- [ ] Python transpilation complete
- [ ] Python tests passing
