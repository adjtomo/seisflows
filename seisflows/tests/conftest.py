"""
Shared pytest configuration for the SeisFlows test suite.
"""
import os
import pytest

from seisflows import ROOT_DIR, logger

# Parent directory of the `seisflows` package as imported by *this* process
# (i.e. this checkout, whichever git worktree/branch it happens to be)
_PACKAGE_ROOT = os.path.dirname(ROOT_DIR)


def pytest_configure(config):
    """
    Make sure subprocesses spawned during the test run import THIS checkout's
    `seisflows` package, not whatever a `pip install -e` editable install
    happens to point at.

    .. note::
        `system.cluster.Cluster.run()` (and its subclasses, e.g. `Slurm`,
        `Pbs`) dispatch each task by spawning a brand new Python process that
        runs `seisflows/system/runscripts/run` via its shebang line. That
        child process's `sys.path` is seeded from the script's own directory
        plus site-packages -- it does NOT inherit the cwd-based shadowing
        that lets a plain `import seisflows` resolve to the current checkout
        when running pytest interactively. If the environment's `seisflows`
        editable install points at a *different* checkout (e.g. you're
        running tests from a git worktree other than the one `pip install
        -e .` was run from), the spawned child silently imports and runs the
        OTHER checkout's code -- including test helper functions, for tests
        that exercise `Cluster.run()` directly. Prepending this checkout's
        root to `PYTHONPATH` (inherited by any subprocess we spawn) avoids
        that entirely.
    """
    existing = [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep)
               if p]
    if _PACKAGE_ROOT not in existing:
        os.environ["PYTHONPATH"] = os.pathsep.join(
            [_PACKAGE_ROOT] + existing)


@pytest.fixture(autouse=True)
def reset_logger_filters():
    """
    A workflow that starts out resuming hides its log until it gets back to
    the dispatch it is collecting (see `Dispatcher.begin_resume`). Tests
    don't always get that far, so don't let one test hide another's log
    """
    yield
    logger.filters.clear()
