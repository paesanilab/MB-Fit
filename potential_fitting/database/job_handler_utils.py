from potential_fitting.utils import SettingsReader
from potential_fitting.exceptions import NoSuchLibraryError
from .psi4_job_handler import Psi4JobHander
from .qchem_job_handler import QchemJobHandler


def get_job_handler(settings_path, logging = True):
    settings = SettingsReader(settings_path)
    if settings.get("energy_calculator", "code") == "psi4":
        return Psi4JobHander(settings_path)
    elif settings.get("energy_calculator", "code") == "qchem":
        return QchemJobHandler(settings_path)
    else:
        raise NoSuchLibraryError(settings.get("energy_calculator", "code"))
