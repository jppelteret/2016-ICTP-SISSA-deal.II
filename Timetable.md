# Course schedule

## Monday 19.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | Introduction<br>First steps [Steps 1,2] | JPP |
| 11:15 | 1.25 hours | Introduction to FEM | LH |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | Solving Poisson's equation [Step 3] <br> Dimension independent programming | JPP |
| 15:45 | 1.25 hours | Exercises, Q&A | JPP, LH |

## Tuesday 20.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | Computing errors | JPP |
| 11:15 | 1.25 hours | Local (adaptive) refinement [Step 6] | JPP |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | Shared memory parallelisation (TBB) | JPP |
| 15:45 | 1.25 hours | Exercises, Q&A | JPP, LH |

## Wednesday 21.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | MPI parallelisation: Shared [quasi step 17/18] | JPP |
| 11:15 | 1.25 hours | Exercises | JPP, LH |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | MPI parallelisation: Distributed [step-40] | LH |
| 15:45 | 1.25 hours | Exercises, Q&A | JPP, LH |

## Thursday 22.03

| Time | Duration | Title | Speaker  |
|:-----|:---------|:------|:---------|
| 09:30 | 1.25 hours | Parameter handler / acceptor; Git workflow | JPP |
| 11:15 | 1.25 hours | Exercises | JPP, LH |
|<td colspan=3>LUNCH</td>|
| 14:00 | 1.25 hours | Automatic differentiation | JPP |
| 15:45 | 1.25 hours | Exercises | JPP |

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
  - Three different possible assembly modes
    - Energy functional (Miehe2014)
    - Residual
    - Full residual + linearisation
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
