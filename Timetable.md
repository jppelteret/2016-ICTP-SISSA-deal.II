# Course schedule

## Monday 19.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | Introduction<br>First steps [Steps 1,2] | JPP |
| 11:15 | 1.25 hours | Introduction to FEM | LH |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | Solving Poisson's equation [Step 3] <br> Dimension independent programming [Step 4] | JPP |
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
- 2017 Lab 05,07
- *TODO*: Lab TBB

## Wednesday 21.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | MPI parallelisation: Shared [Quasi step 17/18] | JPP |
| 11:15 | 1.25 hours | Exercises | JPP, LH |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | MPI parallelisation: Distributed [Step 40] | LH |
| 15:45 | 1.25 hours | Exercises, Q&A | JPP, LH |

#### Exercises
- *TODO*: Lab parallel::shared + PETSc
- *TODO*: Lab parallel::distributed + Trilinos

## Thursday 22.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | Parameter handler / acceptor; Git workflow; Linear operator | JPP |
| 11:15 | 1.25 hours | Time dependent problems [Step 23]; Solution transfer | JPP |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | Automatic differentiation | JPP |
| 15:45 | 1.25 hours | Exercises | JPP |

#### Exercises
- *TODO*: Lab time-dependent problem
- *TODO*: Lab automatic differentiation

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

------------

Information
- http://indico.ictp.it/event/7751/overview
- http://indico.ictp.it/event/7751/other-view?view=ictptimetable

Summary
- DAY 1: Compressed version of 2016 day 1, 2 (step 1-3 + ddim programming) ; Laplace, dim indep.
  - [JP] Step 1,2
  - Break
  - [LUCA] FEM intro
  - Dim inp. programming; Error adaptivity

- DAY 2:
  - ... step 6 [up to end of morning]
  - parallelisation strategies
    - TBB
    - Shared

- DAY 3:
  - Full day: MPI parallelisation (Thus 2016 lectures)
    - Parallel version of Laplace (step 40)

- DAY 4 & 5:
  - Advanced topics:
    - Parameter acceptor
    - Local refinement
    - Manifolds
    - Git workflow for complex libraries
    - Nonlinear problem
    - AD? - Yes (Thurs morn / afternoon)
      - Laplace
        - AD;
      - Elastic energy; compute res + tangent
      - Needs development version (docker)
    - Building up coupled problems? - No; too complex
  - Goal:
    - Cahn-Hilliard with Sacado in parallel
    - Make PR to deal.II library

Introduction
- Who the students are
- What they hope to achieve
