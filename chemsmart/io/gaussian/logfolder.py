import os
import logging
from .output import Gaussian16Output
from chemsmart.io.molecules.dataset import Dataset

logger = logging.getLogger(__name__)
class GaussianLogFolder:
    """Log folder containing all Gaussian log files for assembly into a full db."""

    def __init__(self, folder):
        """:param folder: Parent folder for all log files; type of str"""
        self.folder = folder

    def get_all_logs(self):
        log_paths = []
        for subdir, _dirs, files in os.walk(self.folder):
            for file in files:
                filepath = os.path.join(subdir, file)
                if filepath.endswith('.log'):
                    log_paths.append(filepath)
        return log_paths

    def assemble_full_db(self, db_filename='all.db', **kwargs):
        # write fully assembly db as 'final.db' in self.folder by default
        db_path = os.path.join(self.folder, db_filename)
        self.db_path = db_path

        log_paths = self.get_all_logs()
        for file in log_paths:
            logger.info(f'Parsing Gaussian output file {file}')
            output_file = Gaussian16Output(file)
            dataset = output_file.to_dataset(**kwargs)

            # leaves out files that are PCM SP.log and '#t' route section and premature termination
            if dataset is not None:
                dataset = dataset.remove_duplicates()
                dataset.write(db_path, append=True)

            logger.info(f'Done assembling all structures into dataset: {dataset}')

        deduplicated_db = Dataset.from_file(db_path).remove_duplicates()
        # append = False to write over itself so that duplicates are removed
        deduplicated_db.write(db_path, append=False, overwrite=True)
        return deduplicated_db

    def get_total_service_units(self, job_runtime_file='job_runtime.txt'):
        log_paths = self.get_all_logs()
        with open(job_runtime_file, 'w') as f:
            total_core_hours_in_folder = 0
            for file in log_paths:
                logger.info(f'Parsing Gaussian output file {file}')
                output_file = Gaussian16Output(file)
                core_hous = output_file.total_core_hours
                total_core_hours_in_folder += core_hous
                f.write(f'Job: {file:<130} Total time: {core_hous:6.1f} core-hours\n')
            total_core_hours_in_folder = round(total_core_hours_in_folder, 1)
            f.write(f'TOTAL core-hours in folder {self.folder} is: {total_core_hours_in_folder}\n')
        return total_core_hours_in_folder