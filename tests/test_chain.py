from chemsmart.jobs.chain import ChainJob, JobPhase
from chemsmart.jobs.job import Job


class DummyJob(Job):
    TYPE = "dummy"

    def __init__(self, label, complete=False):
        super().__init__(molecule=None, label=label, jobrunner=None)
        self._complete = complete
        self.run_count = 0

    def is_complete(self):
        return self._complete

    def _run(self, **kwargs):
        self.run_count += 1
        self._complete = True


class TestChainJob:
    def test_runs_phases_in_order_and_marks_complete(self):
        first = DummyJob("first")
        second = DummyJob("second")
        chain = ChainJob(
            molecule=None,
            label="chain",
            jobrunner=None,
            phases=[
                JobPhase(name="one", jobs=[first]),
                JobPhase(name="two", jobs=[second]),
            ],
        )

        chain._run()

        assert first.run_count == 1
        assert second.run_count == 1
        assert chain.is_complete()

    def test_stops_when_required_phase_is_incomplete(self):
        first = DummyJob("first")
        second = DummyJob("second")

        def _run_without_completing(**kwargs):
            first.run_count += 1

        first._run = _run_without_completing
        chain = ChainJob(
            molecule=None,
            label="chain",
            jobrunner=None,
            phases=[
                JobPhase(
                    name="one",
                    jobs=[first],
                    require_complete=True,
                    stop_message="one incomplete",
                ),
                JobPhase(name="two", jobs=[second]),
            ],
        )

        chain._run()

        assert first.run_count == 1
        assert second.run_count == 0
        assert not chain.is_complete()

    def test_skip_if_omits_phase_from_run_and_completeness(self):
        required = DummyJob("required")
        skipped = DummyJob("skipped")
        chain = ChainJob(
            molecule=None,
            label="chain",
            jobrunner=None,
            phases=[
                JobPhase(name="required", jobs=[required]),
                JobPhase(
                    name="optional",
                    jobs=[skipped],
                    skip_if=lambda: True,
                ),
            ],
        )

        chain._run()

        assert required.run_count == 1
        assert skipped.run_count == 0
        assert chain.is_complete()

    def test_jobs_factory_builds_later_phase_from_earlier_output(self):
        first = DummyJob("first")
        created = {}

        def _second_jobs():
            created["job"] = DummyJob(f"{first.label}_sp")
            return [created["job"]]

        chain = ChainJob(
            molecule=None,
            label="chain",
            jobrunner=None,
            phases=[
                JobPhase(name="one", jobs=[first]),
                JobPhase(name="two", jobs_factory=_second_jobs),
            ],
        )

        chain._run()

        assert created["job"].label == "first_sp"
        assert created["job"].run_count == 1
        assert chain.is_complete()
