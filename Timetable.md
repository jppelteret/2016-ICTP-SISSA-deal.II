# Course schedule

## Monday 19.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | Introduction<br>First steps [Steps 1,2] | JPP |
| 11:15 | 1.25 hours | Introduction to FEM | LH |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | Solving Poisson's equation [Steps 3,4] | JPP |
| 15:45 | 1.25 hours | Exercises, Q&A | JPP, LH |

#### Exercises
- 2017 Lab 01-04

## Tuesday 20.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 11:15 | 1.25 hours | Local refinement + hanging nodes [Quasi step 6] | JPP |
| 09:30 | 1.25 hours | Local (adaptive); computing errors [Steps 6,7] | LH |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | Shared memory parallelisation (TBB) | JPP |
| 15:45 | 1.25 hours | Exercises, Q&A | JPP, LH |

#### Exercises
- Lab Local refinement
- Lab adaptive refinement
- Lab TBB

## Wednesday 21.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | MPI parallelisation: Shared [Quasi step 17/18] | JPP |
| 11:15 | 1.25 hours | Exercises | JPP, LH |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | MPI parallelisation: Distributed [Step 40] | LH |
| 15:45 | 1.25 hours | Exercises, Q&A | JPP, LH |

#### Exercises
- Lab parallel::shared + PETSc
- Lab parallel::distributed + Trilinos

## Thursday 22.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | Utility classes and functions; Git workflow | JPP |
| 11:15 | 1.25 hours | Time dependent problems [Step 23]; Solution transfer | JPP |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | Automatic differentiation | JPP |
| 15:45 | 1.25 hours | Exercises, Q&A | JPP |

#### Exercises
- Lab time-dependent problem
- Lab automatic differentiation

## Friday 23.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | Q&A, work on projects | JPP, LH |
| 11:15 | 1.25 hours | Q&A, work on projects | JPP, LH |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | Q&A, work on projects | JPP, LH |
| 15:45 | 1.25 hours | Q&A, work on projects | JPP, LH |

#### Project
- Cahn-Hilliard with Sacado in parallel
  - References
    - Miehe, C.; Hildebrand, F. E. & Boger, L. Mixed variational potentials and inherent symmetries of the Cahn-Hilliard theory of diffusive phase separation Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, The Royal Society, 2014, 470, 20130641-2013064.
    DOI: 10.1098/rspa.2013.0641
    - Miehe, C.; Mauthe, S. & Ulmer, H. Formulation and numerical exploitation of mixed variational principles for coupled problems of Cahn-Hilliard-type and standard diffusion in elastic solids International Journal for Numerical Methods in Engineering, Wiley-Blackwell, 2014, 99, 737-762.
    DOI: 10.1002/nme.4700
  - Three (?) different possible assembly modes
    - Energy functional (Miehe2014)
    - Residual
    - Full residual + linearisation (?)
- Make PR to deal.II library
