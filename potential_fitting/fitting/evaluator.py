import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib import rcParams

from potential_fitting.utils import system, files
from potential_fitting.training_set import TrainingSet

class Evaluator:

    def __init__(self, settings, path_to_eval_file):
        self.settings = settings
        self.path_to_eval_file = path_to_eval_file
        self.ts_obj = None
        self.energies = []
        self.rmsd = []

    def calculate_energies(self, parameter_file_path, training_set_file_path, correlation_file_path=None):

        if correlation_file_path is None:
            correlation_file_path = files.init_file(os.path.join(self.settings.get("files", "log_path"), "eval.dat"))

        with open(correlation_file_path, "w") as correlation_file:
            system.call(self.path_to_eval_file, parameter_file_path, training_set_file_path, out_file=correlation_file)

        with open(correlation_file_path, "r") as correlation_file:

            fit_energies = []
            correlation_file.readline()
            for line in correlation_file.readlines():
                if not line.startswith("#"):
                    fit_energies.append(float(line.split()[1]))

        self.ts_obj = TrainingSet.get_training_set_from_xyz_file(training_set_file_path, self.settings, energy_names=["binding_energy", "ref_energy"])

        self.ts_obj.add_energies("fit_energy", fit_energies)

        
    def write_correlation_file(self, correlation_file = "correlation.dat", split_energy = None):
        be = self.ts_obj.get_energies("binding_energy")
        ref = self.ts_obj.get_energies("ref_energy")
        nb = self.ts_obj.get_energies("fit_energy")

        nref = len(ref)
        nnb = len(nb)
        nbe = len(be)

        if nref != nnb or nbe != nref:
            raise PotentialFittingError("Sizes of reference ({}), binding ({}), and calculated ({}) energies don't match in the Evaluator class.".format(nref,nbe,nnb))

        correlation_file_path = files.init_file(correlation_file)

        rmsd = 0.0
        corr_str = ""
        
        low_rmsd = 0.0
        high_rmsd = 0.0

        low_str = ""
        high_str = ""

        if split_energy is None:
            split_energy = 100000.0

        low_energy_corr = []
        high_energy_corr = []
        energy_corr = []
        
        for i in range(nbe):
            line = "{0:16.8f}{1:16.8f}\n".format(ref[i],nb[i])
            number = (ref[i] - nb[i]) * (ref[i] - nb[i])
            energy_list = [ref[i],nb[i]]

            corr_str += line
            rmsd += number
            energy_corr.append(energy_list)
            
            if be[i] > split_energy:
                high_str += line
                high_energy_corr.append(energy_list)
                high_rmsd += number
            else:
                low_str += line
                low_energy_corr.append(energy_list)
                low_rmsd += number


        with open(correlation_file_path,'w') as corr:
            corr.write(corr_str)
        rmsd = math.sqrt(rmsd/float(nbe))
        print("Total RMSD = {0:16.4f}".format(rmsd))

        # Result is all[, low, high]
        self.energies.clear()
        self.energies.append(energy_corr)

        self.rmsd.clear()
        self.rmsd.append(rmsd)

        if split_energy is not None:

            low_corr = "low_" + correlation_file
            high_corr = "high_" + correlation_file
            
            self.energies.append(low_energy_corr)
            self.energies.append(high_energy_corr)

            if len(low_energy_corr) > 0:
                low_correlation_file_path = files.init_file(low_corr)
                with open(low_correlation_file_path,'w') as low_c:
                   low_c.write(low_str)
                low_rmsd = math.sqrt(low_rmsd/float(len(low_energy_corr)))
                self.rmsd.append(low_rmsd)
                print("Low Energy RMSD = {0:16.4f}".format(low_rmsd))
            else:
                self.rmsd.append(0.0)

            if len(high_energy_corr) > 0:
                high_correlation_file_path = files.init_file(high_corr)
                with open(high_correlation_file_path,'w') as high_c:
                   high_c.write(high_str)
                high_rmsd = math.sqrt(high_rmsd/float(len(high_energy_corr)))
                self.rmsd.append(high_rmsd)
                print("High Energy RMSD = {0:16.4f}".format(high_rmsd))
            else:
                self.rmsd.append(0.0)
            
        return self.energies

    def plot(self, do_ttm = False, split_energy = None, correlation_prefix = "correlation"):

        if do_ttm:
            colors = ["#FF8000","#EE220C"]
            labels = ["Reference (kcal/mol)", "TTM-nrg (kcal/mol)"]
        else:
            colors = ["#61D836","#017100"]
            labels = ["Reference (kcal/mol)", "MB-nrg (kcal/mol)"]

        x = []
        y = []

        for i in range(len(self.energies)):
            x.append([])
            y.append([])
            for j in range(len(self.energies[i])):
                x[i].append(self.energies[i][j][0])
                y[i].append(self.energies[i][j][1])

        #Creationg of figure
        fig, axs = plt.subplots(1, 1, figsize=(6,6))
        rcParams['font.family'] = 'Helvetica'
        
        axs.set_xlim([0,50])
        axs.set_ylim([0,50])
        axs.xaxis.set_minor_locator(MultipleLocator(5))
        axs.yaxis.set_minor_locator(MultipleLocator(5))
        axs.xaxis.set_major_locator(MultipleLocator(10))
        axs.yaxis.set_major_locator(MultipleLocator(10))

        axs.tick_params(top=True, bottom=True, left=True, right=True, direction='in', labelsize=15, length=6)
        axs.tick_params(which='minor', top=True, bottom=True, left=True, right=True, direction='in', labelsize=15, length=4)

        axs.set_xlabel(labels[0], fontsize=16)
        axs.set_ylabel(labels[1], fontsize=16)

        legends = []

        if split_energy is None:
            axs.scatter(x[0], y[0], marker='s',color=colors[0], facecolors='none')
            legends.append("")
        else:
            if len(x[1]) > 0:
                axs.scatter(x[1], y[1], marker='s',color=colors[0], facecolors='none')
                legends.append("BE < " + str(int(split_energy)) + " kcal/mol")
            if len(x[2]) > 0:
                axs.scatter(x[2], y[2], marker='^',color=colors[1], facecolors='none')
                legends.append("BE > " + str(int(split_energy)) + " kcal/mol")
        
        axs.legend(legends, fontsize=15, loc=4, frameon=False)
        #axs.text(3, 53, r'a) CH$_4$: TTM-nrg', fontsize=15)

        # x = y line
        xi = np.linspace(0, 100 ,1000, endpoint = True)
        yi = xi
        axs.plot(xi,yi, '--',color ='k')

        #Tight layout
        plt.tight_layout()
        
        #Save images
        figName = correlation_prefix +'.png'
        plt.savefig(figName, format='png', dpi=1000)
        figName = correlation_prefix +'.pdf'
        plt.savefig(figName, format='pdf', dpi=1000)

