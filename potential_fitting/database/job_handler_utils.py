from potential_fitting.utils import SettingsReader
from potential_fitting.exceptions import NoSuchLibraryError
from .psi4_job_handler import Psi4JobHandler
from .qchem_job_handler import QchemJobHandler


def get_job_handler(settings_path):
    """
    Gets a JobHandler that will make jobs for the code given in the settings file.

    Args:
        settings_path       - Local path to settings file.

    Returns:
        A new JobHandler object.
    """

    settings = SettingsReader(settings_path)
    if settings.get("energy_calculator", "code") == "psi4":
        return Psi4JobHandler(settings_path)
    elif settings.get("energy_calculator", "code") == "qchem":
        return QchemJobHandler(settings_path)
    else:
        raise NoSuchLibraryError(settings.get("energy_calculator", "code"))
