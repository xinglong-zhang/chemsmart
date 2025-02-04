# import logging
# from functools import cached_property
#
# import numpy as np
#
# from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
# # from chemsmart.utils.grouper import (
# #     StructuralSelfConsistentGrouper,
# #     StructuralSequentialGrouper,
# # )
#
# logger = logging.getLogger(__name__)
#
#
# class GaussianSAOptJob(GaussianJob):
#     TYPE = "g16saopt"
#
#     def __init__(
#         self,
#         molecules,
#         settings,
#         label,
#         num_structures_to_opt=None,
#         grouper_type="seq",
#         num_grouper_processes=1,
#         proportion_structures_to_use=0.1,
#         **kwargs,
#     ):
#         super().__init__(
#             molecule=molecules, settings=settings, label=label, **kwargs
#         )
#         self.num_structures_to_opt = num_structures_to_opt
#         self.grouper_type = grouper_type
#         self.num_grouper_processes = num_grouper_processes
#
#         # proportion of the traj from last portion to obtain structures
#         self.proportion_to_opt = proportion_structures_to_use
#         last_num_structures = int(
#             round(len(molecules) * proportion_structures_to_use, 1)
#         )
#         self.molecules = molecules[-last_num_structures:]
#
#     @cached_property
#     def num_structures(self):
#         return len(self.molecules)
#
#     @cached_property
#     def all_energies(self):
#         return [structure.energy for structure in self.molecules]
#
#     @cached_property
#     def energies_indices(self):
#         """List of energy indices in ascending energies order."""
#         return np.argsort(self.all_energies)
#         # return sorted(range(len(self.all_energies)), key=lambda k: self.all_energies[k])
#
#     @cached_property
#     def sorted_energies(self):
#         """List of sorted energies in ascending order."""
#         return [self.all_energies[i] for i in self.energies_indices]
#
#     @cached_property
#     def sorted_molecules(self):
#         """List of molecules sorted by ascending energies."""
#         return [self.molecules[i] for i in self.energies_indices]
#
#     @cached_property
#     def unique_structures_sequential_grouper(self):
#         """Return list of unique structures using StructuralSequentialGrouper.
#
#         Skip_structure_reduction is set to True as structures as molecular.
#         """
#         grouper = StructuralSequentialGrouper.from_molecules(
#             molecules=self.sorted_molecules,
#             num_procs=self.num_grouper_processes,
#             skip_structure_reduction=True,
#         )
#         groups, _ = grouper.group()
#         unique_images = [group[0] for group in groups]
#         return unique_images
#
#     @cached_property
#     def unique_structures_self_consistent_grouper(self):
#         """Return list of unique structures using StructuralSelfConsistentGrouper.
#
#         Skip_structure_reduction is set to True as structures as molecular.
#         """
#         grouper = StructuralSelfConsistentGrouper.from_molecules(
#             molecules=self.sorted_molecules,
#             num_procs=self.num_grouper_processes,
#             skip_structure_reduction=True,
#         )
#         groups, group_indices = grouper.group()
#         unique_images = [group[0] for group in groups]
#
#         # sort unique images by ascending order of energies again
#         unique_images_energies = [image.energy for image in unique_images]
#         unique_images_sorted = [
#             unique_images[i] for i in np.argsort(unique_images_energies)
#         ]
#         return unique_images_sorted
#
#     @property
#     def unique_structures(self):
#         if self.grouper_type == "seq":
#             return self.unique_structures_sequential_grouper
#
#         if self.grouper_type == "sc":
#             return self.unique_structures_self_consistent_grouper
#
#         raise ValueError("Grouper types available are: 'seq' or 'sc'.")
#
#     @property
#     def unique_structures_energies(self):
#         return [structure.energy for structure in self.unique_structures]
#
#     @property
#     def num_unique_structures(self):
#         return len(self.unique_structures)
#
#     def _prepare_all_jobs(self):
#         jobs = []
#         logger.debug(
#             f"Number of structures used for optimization: {self.num_unique_structures}\n"
#         )
#         logger.debug(f"Unique structures: {self.unique_structures}")
#         logger.debug(
#             f"Unique structures energies: {self.unique_structures_energies}"
#         )
#         for i in range(self.num_unique_structures):
#             label = f"{self.label}_c{i + 1}"  # 1-indexed for structures
#             jobs += [
#                 GaussianGeneralJob(
#                     molecule=self.unique_structures[i],
#                     settings=self.settings,
#                     label=label,
#                     skip_completed=self.skip_completed,
#                 )
#             ]
#         return jobs
#
#     @property
#     def all_structures_opt_jobs(self):
#         return self._prepare_all_jobs()
#
#     @property
#     def last_run_job_index(self):
#         return self._check_last_finished_job_index()
#
#     @property
#     def incomplete_structure_opt_jobs(self):
#         return [
#             job
#             for job in self.all_structures_opt_jobs
#             if not job.is_complete()
#         ]
#
#     def _check_last_finished_job_index(self):
#         for i, job in enumerate(self.all_structures_opt_jobs):
#             if not job.is_complete():
#                 return i
#         # If all complete
#         return self.num_unique_structures
#
#     def _run_all_jobs(self, jobrunner):
#         for job in self.all_structures_opt_jobs[: self.num_structures_to_opt]:
#             job.run(jobrunner=jobrunner)
#
#     def _run(self, jobrunner, **kwargs):
#         self._run_all_jobs(jobrunner=jobrunner)
#
#     def is_complete(self):
#         return self._run_all_sa_opt_jobs_are_complete()
#
#     def _run_all_sa_opt_jobs_are_complete(self):
#         return all(
#             job.is_complete()
#             for job in self.all_structures_opt_jobs[
#                 : self.num_structures_to_opt
#             ]
#         )
