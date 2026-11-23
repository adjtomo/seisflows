# Bug: a failed L-BFGS line search restarts forever

Internal note (found 2026-11-20 in a student's Miyabi inversion log). Kept
next to `NONBLOCKING_MASTER_JOB_CHANGELOG.md` so it can be tracked, and so it
can be turned into a GitHub issue / PR description for upstream.

**Status**: fixed in the working tree, not committed.
**Affects**: any inversion with `optimize: LBFGS` whose line search fails.
`main` has the same bug, it is not caused by the non-blocking work.
**Cost in the observed run**: from the first restart (2026-11-18 13:20) to the
end of the log (2026-11-21 06:11) the same line search ran three times. The
second and third repeats took about 41 hours (20 and 21) for no new
information, and it would never have ended.

## Symptom

Iteration 4 of the student's inversion (`/Users/prof/Scratch/sflog.txt`):

1. L-BFGS gives a search direction (`sflog.txt:965-969`, 53.57 degrees from the
   previous one) and the first line search fails after 5 steps
   (`sflog.txt:1122`, "backtracking line search has failed").
2. The workflow restarts the line search (`sflog.txt:1124-1126`). The restarted
   line search then fails after 5 steps (`sflog.txt:1271`, "bracketing line
   search has failed") and restarts again (`sflog.txt:1274`), and again at
   `sflog.txt:1423`, and a fourth line search is running when the log ends.
3. Every restarted line search is identical. The step lengths and misfits
   (f = 4.847, 4.710, 4.706, 4.665 for alpha = 1.271, 0.636, 0.318, 0.159,
   against 4.659 for the starting model) repeat to four digits at
   `sflog.txt:1195, 1231, 1267`, `1344, 1380, 1416` and `1493, 1529, 1565`.
4. Every failure logs the same angle between the search direction and the
   gradient, `theta: 0.935` (radians, 53.6 degrees) at `sflog.txt:1123, 1272,
   1421`, and then "search direction above threshold 0.001, attempting
   restart".

## Cause

The line search can only be restarted if the search direction is not already
steepest descent, because a steepest descent line search that fails cannot be
helped by restarting it. The code that decides this and the code that restarts
do not agree on what a restart does to the search direction:

| Step | Code | What it does |
|---|---|---|
| Line search fails | `Inversion.update_line_search` FAIL branch, `seisflows/workflow/inversion.py:703-729` | Asks the optimizer, then restarts |
| Can we restart? | `Gradient.attempt_line_search_restart`, `seisflows/optimize/gradient.py:588-619` | Loads `g_new` and `p_new` from disk and allows a restart if the angle between `p_new` and `-g_new` is above 1e-3 |
| Restart | `LBFGS.restart`, `seisflows/optimize/LBFGS.py:195-234` | Clears the L-BFGS memory and the line search history |
| Set up the line search again | `Inversion.initialize_line_search` -> `Gradient.initialize_search` / `calculate_step_length` / `compute_trial_model`, `seisflows/optimize/gradient.py:365, 414, 506` | Reads `p_new` **from disk** and takes a step along it |

`LBFGS.restart` never touched `p_new` on disk. Its docstring said that "the
workflow will need to call `compute_direction`" afterwards, and the workflow
does not. `compute_direction` is only called in
`Inversion.evaluate_gradient_from_kernels` (`inversion.py:542`), once per
iteration, before the line search. So after a restart:

- `p_new` on disk is still the L-BFGS direction that just failed,
- `initialize_line_search` makes the same trial model (same `p_new`, and the
  same first step length from `step_len_init`, 1.271 in the log), so the
  restarted line search repeats the one that failed,
- it fails the same way, and `attempt_line_search_restart` loads the same
  `p_new`, finds it is 53.6 degrees from `-g_new`, and allows another restart.

Nothing in the loop changes between restarts, so it never ends. In a blocking
workflow this would recurse until the walltime ran out, and in a non-blocking
one it resubmits jobs.

NLCG does it correctly: `NLCG.restart` (`seisflows/optimize/NLCG.py:151-162`)
saves `p_new = -g_new` before returning. After that
`attempt_line_search_restart` sees steepest descent, refuses a second restart,
and the workflow stops with "has run out of fallback options"
(`inversion.py:721-729`). `LBFGS.restart` was missing the same thing.

## Fix

`LBFGS.restart` now writes `p_new = -g_new` to disk, as NLCG does. It also sets
`_LBFGS_iter = 1`, so the restarted iteration counts as the first L-BFGS
iteration (a steepest descent one) and the next iteration builds the memory
again, instead of taking a second steepest descent step. NLCG sets its
counter the same way. With this:

1. the restarted line search takes a step along `-g_new`, a different search
   from the one that failed,
2. if that fails as well, `attempt_line_search_restart` finds the angle is 0,
   returns False and the workflow stops cleanly.

The gradient is not preconditioned in a restart, the same as the first L-BFGS
iteration (`LBFGS.compute_direction`, `LBFGS.py:155-159`) and NLCG.

Files changed:

- `seisflows/optimize/LBFGS.py`: `restart` (and its docstring)
- `seisflows/tests/test_optimize.py`:
  `test_optimize_restart_resets_search_direction` (L-BFGS and NLCG: after a
  restart `p_new == -g_new` and a second restart is refused),
  `test_lbfgs_restart_clears_memory`

Not changed: `Inversion.update_line_search`. It has no way to know what an
optimizer's `restart` does to `p_new`, so the contract is on the optimizer.

## Things this does not fix

- **Why the L-BFGS direction fails.** Even alpha = 0.159 raises the misfit
  (4.665 against 4.659). `LBFGS._check_status` (`LBFGS.py:343-365`) does accept
  the direction as a descent direction (the angle between `g` and the
  inverse-Hessian-times-`g` is 53.6 degrees, under 90), so a small enough step
  should reduce the misfit if the gradient is accurate. That it does not
  suggests the gradient, or the misfit along it, is noisy (e.g., the windows
  or measurements changing between models, or a gradient that is not
  consistent with the misfit). The restart fix gives steepest descent a
  chance but does not address that.
- **A stalled inversion.** Misfits of 4.786, 4.684, 4.665, 4.659 over the first
  four iterations: 2.1%, then 0.4%, then 0.13%. After the fix the run is
  likely to stop after the steepest descent line search fails too, which is
  then a modeling decision (windows, regularization, step), not a code one.
- **The restart repeats `step_len_init`.** The restarted line search is
  treated as a first iteration (`step_lens.count(0) == 1` in
  `Gradient.calculate_step_length`, `gradient.py:436`), so it applies
  `step_len_init` (0.05 in the student's file) to the new direction. This
  is the right behavior for steepest descent, and is just worth knowing.
- **A periodic restart** (`_LBFGS_iter > lbfgs_max`,
  `LBFGS.compute_direction`, `LBFGS.py:162-167`) used to leave the counter at
  0, so the next iteration was also a "first iteration" and also steepest
  descent. It now leaves it at 1. `lbfgs_max` defaults to infinity, so this
  only matters when it is set.

## The running job

The student's job keeps submitting line searches (each is 54 simulations).
To stop it, delete the resume job first and then the compute jobs, otherwise
the resume job starts when the compute job ends and submits the next one:

```bash
qdel <resume job id>
qdel <compute array job id>
```

Then update SeisFlows and run `seisflows submit` in the working directory. The
dispatcher treats tasks that did not finish as failed and runs them again
(`rerun` times). This resume path after a manual `qdel` has not been
tested. The `p_new` on disk is still the failed L-BFGS direction until the next
failure, which now restarts with steepest descent, so the first resumed line
search will fail again before the fix takes effect.
