import os

from .job_handler import JobHandler

class Psi4JobHandler(JobHandler):

    def get_job_template_path(self):
        """
        Gets the path to the template job file for this JobHandler.

        Args:
            None.

        Returns:
            Absolute path to job template file.
        """

        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "psi4_job_template.py")
