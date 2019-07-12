from .psi4_calculator import Psi4Calculator
from .qchem_calculator import QchemCalculator
from potential_fitting.utils import SettingsReader

def get_calculator(settings_path, logging = True):
    settings = SettingsReader(settings_path)
    if settings.get("energy_calculator", "code") == "psi4":
        return Psi4Calculator(settings_path, logging)
    elif settings.get("energy_calculator", "code") == "qchem":
        return QchemCalculator(settings_path, logging)
    else:
        raise NoSuchLibraryError(settings.get("energy_calculator", "code"))
