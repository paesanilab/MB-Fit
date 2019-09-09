import os

from .job_handler import JobHander

class QchemJobHandler(JobHander):

    def get_job_template_path(self):
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "qchem_job_template.py")
