# Claude contexts

Working notes written during Claude Code sessions on SeisFlows, kept so that
the context behind some changes is not lost. They are not part of the
built documentation (`claude_contexts` is in `exclude_patterns` of `conf.py`),
and they are not maintained as docs: parts of them are dated.

- `NONBLOCKING_MASTER_JOB_CHANGELOG.md`: the tracker for the non-blocking main
  job work (design decisions, what was done in which order, what was validated
  on which cluster, and what is not done), and later work that came out of it:
  partition parameters, pack mode, `seisflows clean`.
- `BUG_LBFGS_RESTART_LOOP.md`: write-up of the L-BFGS line search that
  restarted forever, with log and code references.
- `TESTING_FRAMEWORK_STATUS.md`: the status of the testing framework project
  that the non-blocking work continues from.

Features that are designed but not implemented are in
[`docs/planned_features.rst`](../planned_features.rst).
