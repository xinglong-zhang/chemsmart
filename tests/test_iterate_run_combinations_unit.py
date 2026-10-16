"""
Direct unit tests for the multiprocessing-heavy
``IterateJobRunner.run_combinations`` in ``chemsmart.jobs.iterate.runner``.

``multiprocessing.Manager`` and ``multiprocessing.Process`` are replaced
with lightweight, deterministic fakes so the timeout/terminate/kill
watchdog logic, the queue-draining loops (including their generic-
exception branches), and the missing/failed-result bookkeeping can be
exercised without spawning real subprocesses or waiting on real
timeouts.
"""

import queue
from unittest.mock import MagicMock, patch

from chemsmart.jobs.iterate.runner import (
    IterateCombination,
    IterateJobRunner,
    IterateMoleculePool,
)


class FakeQueue:
    """Stand-in for a ``multiprocessing.Manager().Queue()``.

    ``get_nowait_plan`` is a list consumed one entry per call to
    ``get_nowait``: ``None`` means "return the next queued item (or
    raise ``queue.Empty`` if none)", and any exception instance/class
    means "raise this instead". Once the plan is exhausted, falls back
    to the "return next item or raise Empty" behavior.
    """

    def __init__(self, get_nowait_plan=None):
        self._items = []
        self._plan = list(get_nowait_plan or [])

    def put(self, item):
        self._items.append(item)

    def get_nowait(self):
        if self._plan:
            effect = self._plan.pop(0)
            if effect is not None:
                raise effect
        if self._items:
            return self._items.pop(0)
        raise queue.Empty()


class FakeProcess:
    """Stand-in for a ``multiprocessing.Process``.

    ``alive`` controls every ``is_alive()`` call for this process
    (kept simple: a process is either always "still running" or
    already "finished", which is all the watchdog logic needs).
    """

    def __init__(self, target=None, args=(), daemon=False, alive=False):
        self.target = target
        self.args = args
        self.daemon = daemon
        self.pid = 4242
        self._alive = alive
        self.terminated = False
        self.killed = False
        self.joined_timeouts = []

    def start(self):
        pass

    def is_alive(self):
        return self._alive

    def join(self, timeout=None):
        self.joined_timeouts.append(timeout)

    def terminate(self):
        self.terminated = True

    def kill(self):
        self.killed = True


def _make_pool_and_combo(label_suffix="a"):
    skeleton = MagicMock()
    substituent = MagicMock()
    pool = IterateMoleculePool(
        skeletons=[skeleton], substituents=[substituent]
    )
    combo = IterateCombination(
        skeleton_idx=0,
        skeleton_label=f"skel{label_suffix}",
        skeleton_link_index=1,
        skeleton_indices=None,
        substituent_idx=0,
        substituent_label=f"sub{label_suffix}",
        substituent_link_index=1,
    )
    return pool, combo


def _patched_manager(fake_queue):
    fake_manager = MagicMock()
    fake_manager.Queue.return_value = fake_queue
    fake_manager.shutdown.return_value = None
    return fake_manager


class TestRunCombinationsQueueDrainExceptions:
    def test_generic_exception_from_get_nowait_breaks_drain_loop(self):
        """A non-Empty exception raised by the queue must be swallowed
        the same way queue.Empty is (line 341-342), not propagated."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()
        result_mol = MagicMock()

        fake_queue = FakeQueue(get_nowait_plan=[RuntimeError("transient")])
        fake_queue.put((combo.label, result_mol))

        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        # The RuntimeError on the first get_nowait() call just breaks the
        # drain loop for that pass; the result is still picked up on a
        # later drain (main loop or final double-check).
        assert dict(results)[combo.label] is result_mol

    def test_generic_exception_in_final_double_check_is_swallowed(self):
        """Same as above, but for the final "double check queue one
        last time" drain after the main loop exits (line 393-394)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()

        # First call (inside the main loop's single pass): queue.Empty
        # (nothing arrived yet). Process is already "finished", so the
        # main loop exits after one pass. Second call (the final
        # double-check): a generic exception, which must be swallowed
        # rather than propagated.
        fake_queue = FakeQueue(
            get_nowait_plan=[queue.Empty(), RuntimeError("boom")]
        )
        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        # Never arrived on the queue at all -> counted as a crashed/
        # missing worker.
        assert dict(results)[combo.label] is None

    def test_final_double_check_picks_up_late_arriving_result(self):
        """A result that wasn't yet on the queue during the main loop's
        last drain, but is present after the loop exits, is still
        picked up by the final double-check (line 390)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()
        result_mol = MagicMock()

        # First get_nowait() call (main loop's only pass): Empty.
        # Second call (final double-check): the real result.
        fake_queue = FakeQueue(get_nowait_plan=[queue.Empty()])
        fake_process = FakeProcess(alive=False)

        def _late_put():
            fake_queue.put((combo.label, result_mol))

        # Simulate the result "arriving" between the main loop's drain
        # and the final double-check by seeding it right after the
        # first (Empty) call has been consumed from the plan.
        original_get_nowait = fake_queue.get_nowait
        call_count = {"n": 0}

        def _get_nowait_with_late_arrival():
            call_count["n"] += 1
            if call_count["n"] == 2:
                _late_put()
            return original_get_nowait()

        fake_queue.get_nowait = _get_nowait_with_late_arrival

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        assert dict(results)[combo.label] is result_mol


class TestRunCombinationsTimeoutWatchdog:
    def test_stuck_process_is_terminated_then_killed(self):
        """A process still alive after both the timeout window and a
        terminate() attempt is force-killed (line 365-370)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()

        fake_queue = FakeQueue()
        # Always "alive": guarantees both is_alive() checks in the
        # watchdog branch (before terminate, and after terminate+join)
        # see a still-running process.
        fake_process = FakeProcess(alive=True)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            # timeout=0 -> immediately considered overdue on the very
            # first watchdog pass.
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=0
            )

        assert fake_process.terminated is True
        assert fake_process.killed is True
        assert dict(results)[combo.label] is None

    def test_duplicate_label_explicit_none_not_appended_twice(self):
        """Two combinations that happen to share the same label, both
        resolving to an explicit (non-timed-out) None result, must
        only add that label to failed_labels once (line 413->401)."""
        runner = IterateJobRunner()
        pool, combo_a = _make_pool_and_combo()
        combo_b = IterateCombination(
            skeleton_idx=0,
            skeleton_label=combo_a.skeleton_label,
            skeleton_link_index=combo_a.skeleton_link_index,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label=combo_a.substituent_label,
            substituent_link_index=combo_a.substituent_link_index,
        )
        assert combo_a.label == combo_b.label

        fake_queue = FakeQueue()
        fake_queue.put((combo_a.label, None))
        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo_a, combo_b], nprocs=1, timeout=120
            )

        assert [mol for _, mol in results] == [None, None]


class TestRunCombinationsResultBookkeeping:
    def test_missing_result_counts_as_crashed_worker(self):
        """A process that exits without ever writing to the queue is
        reported as a crashed worker (line 403-407)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()

        fake_queue = FakeQueue()  # never gets anything put on it
        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        assert dict(results)[combo.label] is None

    def test_explicit_none_result_counts_as_failed_not_timed_out(self):
        """A worker that legitimately finishes and reports None (e.g.
        analyzer found no valid geometry) is counted as a failure, not
        a timeout (line 413-414)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()

        fake_queue = FakeQueue()
        fake_queue.put((combo.label, None))
        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        assert dict(results)[combo.label] is None
