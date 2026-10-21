# Review evidence

Sample retrieved on 16 October 2026 from 12 CHEMSMART PR discussions: #391, #475, #524, #563, #677, #679, #681, #687, #697, #699, #700, and #703. This is a selected historical sample, not an exhaustive analysis of all reviews. The excerpts below are authored by `xinglong-zhang`; bot reports and contributor replies are not evidence of his preferences. PR #700 contained approval but no substantive standards; #679 also requested ignoring editor-local configuration. Read the current thread and code before applying an example literally.

## Contents

- [Interpretation safeguards](#interpretation-safeguards)
- [User clarification: class layout](#user-clarification-class-layout)
- [PR #391](#pr-391)
- [PR #475](#pr-475)
- [PR #524](#pr-524)
- [PR #563](#pr-563)
- [PR #677](#pr-677)
- [PR #681](#pr-681)
- [PR #687](#pr-687)
- [PR #697](#pr-697)
- [PR #699](#pr-699)
- [PR #703](#pr-703)

## Interpretation safeguards

- Give later corrections precedence. In #563, Xinglong withdrew a proposal to move CDXML object operations to utilities after recognising their use of `self` (r3529267520).
- In the same PR he explicitly corrected an earlier objection to class methods returning multiple objects (r3536478367). Python does not constrain classmethod return types. Preserve the valid concern about incomplete output objects constructed from settings, not the withdrawn blanket return-type restriction.
- Public APIs are a recurring preference, not a requirement to expose every internal hook. Inspect callers and the supported interface.
- General safety checks, prioritisation, and output format in SKILL.md are operational synthesis. They are not verbatim user rules. Historical names/defaults require current verification.

## User clarification: class layout

In the conversation on 17 October 2026 (Asia/Singapore), the user asked whether the private static helpers preceding `GaussianMECPJob.__init__` should be moved below it. Apply constructor-first ordering as the requested readability preference. The suggested broader ordering of factories, properties, public operations, and private helpers is review guidance, not a Python language rule. Keep legitimate static methods; do not infer that object-oriented design requires instance methods for every pure calculation.

## PR #391

[Source](https://github.com/xinglong-zhang/chemsmart/pull/391#discussion_r2797193286) — `chemsmart/jobs/orca/runner.py`

> why is this copied as part of the _write_input() private method? if orca jobs other than NEB reference specific .xyz file, this would fail? The original logic of copying outside should be better.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/391#discussion_r2797196687) — `chemsmart/jobs/orca/settings.py`

> should be method of ORCAInputWriter class?

[Source](https://github.com/xinglong-zhang/chemsmart/pull/391#discussion_r2797205423) — `chemsmart/jobs/settings.py`

> these pertain to QMMM settings, you removed these so they are handled correctly?

## PR #475

[Source](https://github.com/xinglong-zhang/chemsmart/pull/475#discussion_r2797277368) — `chemsmart/jobs/grouper/connectivity.py`

> here and in other places, this writing appear multiple times, possible to have this as a separate ResultsRecoder class or a class method for Grouper.record_results(output_format='.xlsx') for example? will reduce code redundancy

[Source](https://github.com/xinglong-zhang/chemsmart/pull/475#discussion_r2798152825) — `tests/test_groupers.py`

> pytest fixtures should be in conftest.py

[Source](https://github.com/xinglong-zhang/chemsmart/pull/475#discussion_r2808397616) — `chemsmart/cli/grouper/grouper.py`

> should be obtainable from Folder object and use general output_files so that Gaussian and ORCA can be equally used.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/475#discussion_r3090307083) — `chemsmart/io/gaussian/output.py`

> are there tests added for these two new properties? can add the tests directly to the existing ones that are already testing the Gaussian log files.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/475#discussion_r3090367090) — `chemsmart/jobs/grouper/energy.py`

> is it possible to exclude the molecules without energy instead of raising an error, since in a list of output files, maybe some did not complete properly but we still want to proceed with grouping for those that completed normally

[Source](https://github.com/xinglong-zhang/chemsmart/pull/475#discussion_r3090405130) — `chemsmart/jobs/grouper/isomorphism.py`

> this is what i mean - each grouper has a grouper.record() method that calls the recorder to record the results, each grouper can have a private method for different recording if needed.

## PR #524

[Source](https://github.com/xinglong-zhang/chemsmart/pull/524#discussion_r3006065481) — `chemsmart/cli/config.py`

> @copilot , could this `_update_shell_config` be moved to utils.io.py simiarly for `windows_update_env`?

[Source](https://github.com/xinglong-zhang/chemsmart/pull/524#pullrequestreview-4026680615)

> @copilot, please check if option `--conda-path` is still used? Seems that this has already been removed in the codes but not updated in RTD

[Source](https://github.com/xinglong-zhang/chemsmart/pull/524#issuecomment-4153134315)

> @ZihanL1 and @Lewen-WangCC, this branch now works well for installation on gitbash and conda powershell, please check and once okay, i will merge into main. I have tested this on Windows virtual machine on Mac and it works well.

## PR #563

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3259537157) — `chemsmart/cli/pka.py`

> If every time this (and other related functions are called, it requires a pKa output table as an argument, then perhaps it's best to construct as a class?

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3259707859) — `chemsmart/jobs/gaussian/pka.py`

> seems to have too many things you are dealing with, the project_settings and job_settings are only merged in CLI. the Job class only requires the molecule and settings.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3259756651) — `chemsmart/jobs/gaussian/pka.py`

> instead of having this here, perhaps we can have a out.validate_imaginary_frequencies() method, as this may be used elsewhere too.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3259790095) — `chemsmart/jobs/gaussian/settings.py`

> this value is in water and will it be used wrongly if value is needed in other solvents?

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3440489346) — `chemsmart/cli/orca/pka.py`

> can be a @classmethod of ORCApKaJobSettings class since ORCApKaJobSettings is returned.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3529191654) — `chemsmart/cli/pka.py`

> use filepath consistently instead of file_path. check everywhere!

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3529235191) — `chemsmart/io/orca/output.py`

> as previously. I assume that when you address Gaussian codes, you will address ORCA codes too, so i am less meticulous when checking ORCA codes

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3529267520) — `chemsmart/io/file.py`

> never mind the first comment - the more i read the methods, the more i realise they belong to operations on the object self, which is a CDXML file.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3535253580) — `chemsmart/io/gaussian/output.py`

> IN this case, it should not be a @classmethod, since this will return the Gaussian16pKaOutput object, if i then want to access the normal attributes in this object, it cannot be found, if it is built from settings instead of from the log file.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/563#discussion_r3536478367) — `chemsmart/io/gaussian/output.py`

> According to GPT, this is okay and I misunderstood:
>
> > it is okay for a @classmethod to create and return multiple class objects instead of a single object. Python does not require a classmethod to return cls(...); it can return any useful object: one instance, a list, a dict, a tuple, etc.
> >
> > In your case, for_pka_species is reasonable because it uses cls(...), so subclasses of ORCApKaOutput would also work correctly.

## PR #677

[Source](https://github.com/xinglong-zhang/chemsmart/pull/677#discussion_r3901452730) — `chemsmart/jobs/grouper/base.py`

> I would imagine this to be useful as a public method? if you pass this grouper object to something else, wouldn't this be useful:
>
> > class SomeObject:
> >     def self.__init__(grouper)
> >         self.grouper = grouper
> >
> >     def some_method():
> >         groups = self.grouper.group_by_threshod()

[Source](https://github.com/xinglong-zhang/chemsmart/pull/677#discussion_r3901463877) — `chemsmart/jobs/grouper/connectivity.py`

> This looks like it should be a Molecule class function/method, so that you will have:
>
> `graph = Molecule.to_rdkit_connectivity_graph()`
>
> or a Molecule property so that:
>
> `graph = Molecule.rdkit_connectivity_graph`

[Source](https://github.com/xinglong-zhang/chemsmart/pull/677#discussion_r3901472147) — `chemsmart/jobs/grouper/formula.py`

> This is a private method being called outside the class, consider converting it into public method.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/677#discussion_r3901475068) — `chemsmart/jobs/grouper/isomorphism.py`

> Can't we use Molecule.to_rdkit() that already exists? or modify that one?

[Source](https://github.com/xinglong-zhang/chemsmart/pull/677#discussion_r3912474184) — `chemsmart/jobs/grouper/isomorphism.py`

> Please do not use `ChemSmart`, use `CHEMSMART` instead! I have emphasized this close to 10 times already.

## PR #681

[Source](https://github.com/xinglong-zhang/chemsmart/pull/681#discussion_r3920539477) — `docs/source/pymol-electronic-structure.rst`

> Use Gaussian default of 0 corresponding to 80x80x80.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/681#discussion_r3920740159) — `chemsmart/cli/mol/mol.py`

> This help messages still says -2?

## PR #687

[Source](https://github.com/xinglong-zhang/chemsmart/pull/687#discussion_r3986832526) — `chemsmart/cli/gaussian/link.py`

> Seems that these click options can be wrapped as a decorator that can be reused.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/687#discussion_r3986852014) — `chemsmart/io/gaussian/output.py`

> regex pattern may be placed in utils.repatterns.py

[Source](https://github.com/xinglong-zhang/chemsmart/pull/687#discussion_r3986939087) — `chemsmart/jobs/gaussian/runner.py`

> scratch_parent_folder attribute is not created for the job? I prefer that usable attributes be created in the class, rather than at runtime (due to the fact that if the latter, one may not be able to access that attribute as you don't know if the function creating the attribute is run or not), thus I will usually avoid using getattr() but instead have the attribute set up in the class constructor.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/687#discussion_r3986942425) — `chemsmart/jobs/gaussian/runner.py`

> If self.scratch_dir is None, wouldn't this cause error?

[Source](https://github.com/xinglong-zhang/chemsmart/pull/687#discussion_r4013492682) — `chemsmart/cli/chemsmart`

> Is this file permission changed?

## PR #697

[Source](https://github.com/xinglong-zhang/chemsmart/pull/697#pullrequestreview-5260968266)

> I think in general this is not quite right, each Molecule object should have the distances, angles and dihedrals matrix. The IRC trajectory/output file can already extract the Molecules object, so it is needed that each Molecule object have these matrices. As discussed, there should be Molecule.distances_matrix, Molecule.angles_matrix and Molecule.dihedrals_matrix properties. @laurag172213-bit

## PR #699

[Source](https://github.com/xinglong-zhang/chemsmart/pull/699#discussion_r4092120232) — `chemsmart/analysis/thermochemistry.py`

> This looks like it can subclass the OutputFile, so that the self.contents becomes/are already available from parent class.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/699#discussion_r4092134783) — `chemsmart/analysis/thermochemistry.py`

> seems like these can be refactored as class properties, instead of simple variables from __init__.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/699#discussion_r4092139718) — `chemsmart/analysis/thermochemistry.py`

> Should this self.normal_termination always be set to True? Shouldn't this be determined from the content of the file?

[Source](https://github.com/xinglong-zhang/chemsmart/pull/699#discussion_r4092169303) — `chemsmart/analysis/thermochemistry.py`

> Why does _load_molecule() in parent class use  Molecule.from_filepath(filename) but here it uses self.file_object.molecule? If self.file_object.molecule works in the parent class, then this need function needs not be re-written again, since the parent one will be used.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/699#discussion_r4092257375) — `chemsmart/jobs/gaussian/mecp.py`

> It looks like this can be written as a more general reusable utility function, e.g., `strip_parenthesized_read(match, word)` so that word can be set as any word that we want to strip/remove.

[Source](https://github.com/xinglong-zhang/chemsmart/pull/699#discussion_r4092267377) — `chemsmart/jobs/gaussian/mecp.py`

> will this match `nosymm` or `nosymmetry` as well?

## PR #703

[Source](https://github.com/xinglong-zhang/chemsmart/pull/703#discussion_r4290310573) — `tests/test_GaussianIO.py`

> why is `stable=opt` not in `output.route_string`? It should also be captured in the output route_string isn't it?

[Source](https://github.com/xinglong-zhang/chemsmart/pull/703#discussion_r4290338471) — `chemsmart/io/gaussian/output.py`

> why "stable=opt" in route, then the `route_string` is not returned?
