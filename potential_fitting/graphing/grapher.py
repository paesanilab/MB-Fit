import matplotlib.pyplot as plt
import numpy
from collections import OrderedDict
import csv
import math

from potential_fitting.utils import files

class Grapher:

    def make_distribution_graph(self, data_sets, data_set_names, energy_name, output_file, num_bins=50, max_percent=0.97):

        plt.clf()

        energies_dict = OrderedDict()

        for data_set_name, data_set in zip(data_set_names, data_sets):
            energies_dict[data_set_name] = sorted(data_set.get_energies(energy_name))

        # the smallest energy in any of the data sets is the lower bound of our energy range.
        minimum = min([energies[0] for energies in energies_dict.values()])

        # the largest energy  at the max_percent percentile in any data set will be the energy after which
        # energies will be placed into an overflow bucket.
        # any energies > max_percent_energy will be placed into an overflow bin.
        try:
            max_percent_energy = max([energies[math.ceil(len(energies) * max_percent)] for energies in energies_dict.values()])
        except IndexError:
            max_percent_energy = max([energies[-1] for energies in energies_dict.values()])

        # the largest energy in any of the data sets is the upper bound of our energy range.
        maximum = max([energies[-1] for energies in energies_dict.values()])

        bin_width = (max_percent_energy - minimum) / num_bins

        min_bin = minimum

        max_bin = max_percent_energy

        overflow_max_bin = max_percent_energy

        # create the bins between min_bin and max_bin

        bins = list(numpy.linspace(min_bin, max_bin, num_bins))

        # the last bin will hold energies between max_bin - 1 and big_max

        bins.append(overflow_max_bin)

        # colors to use for each plot, should start with dark colors and go to lighter ones

        colors = [(0, 0, 0, 1), (0, 0, 1, 0.7), (0, 1, 0, 0.5), (1, 0, 0, 0.3)]

        for index, energies in enumerate(energies_dict.values()):
            plt.hist(energies, bins=bins, color=colors[index])

        # list of labels for data in each input file

        legend_labels = list(energies_dict.keys())

        plt.legend(legend_labels)

        # set the axis labels

        plt.xlabel("{} energy (kcal/mol)".format(energy_name))

        plt.ylabel("Number of Configurations")


        # Title of the chart

        plt.title("Energies Distribution")

        plt.xlim(min_bin, max_bin + bin_width)

        num_ticks = 10

        ticks = list(numpy.linspace(min_bin, max_bin, num_ticks))

        labels = [str(round(tick, 2)) for tick in ticks]

        ticks.append(max_bin)

        labels.append(str(round(max_bin, 2)) + "+")

        plt.xticks(ticks, labels)

        files.init_file(output_file)
        plt.savefig(output_file)

    def make_correlation_graph(self, data_set, ref_energy_name, energy_names, output_file):

        plt.clf()

        ref_energies = data_set.get_energies(ref_energy_name)

        energies_dict = OrderedDict()

        for energy_name in energy_names:
            energies_dict[energy_name] = data_set.get_energies(energy_name)

        colors = [(0, 0, 0, 1), (0, 0, 1, 0.7), (0, 1, 0, 0.5), (1, 0, 0, 0.3)]

        for index, energies in enumerate(energies_dict.values()):
            plt.scatter(ref_energies, energies,
                        color=colors[index], s=5, alpha=0.5)
        # Adding a legend
        plt.legend(energy_names)

        min_ref = min(ref_energies)
        max_ref = max(ref_energies)

        # plot reference line y=x.
        x = [min_ref, max_ref]
        y = x

        plt.plot(x, y, c='orange', alpha=0.5)


        #Adding axes titles
        plt.xlabel("{} energy [Kcal/mol]".format(ref_energy_name))
        plt.ylabel("energy [Kcal/mol]")

        files.init_file(output_file)
        plt.savefig(output_file)


    def make_error_graph(self, data_set, ref_energy_name, energy_names, output_file):

        plt.clf()

        ref_energies = data_set.get_energies(ref_energy_name)

        energies_dict = OrderedDict()

        for energy_name in energy_names:
            energies_dict[energy_name] = data_set.get_energies(energy_name)

        colors = [(0, 0, 0, 1), (0, 0, 1, 0.7), (0, 1, 0, 0.5), (1, 0, 0, 0.3)]

        for index, energies in enumerate(energies_dict.values()):
            plt.scatter(ref_energies, [energy - ref_energy for energy, ref_energy in zip(energies, ref_energies)],
                        color=colors[index], s=5, alpha=0.5)
        # Adding a legend
        plt.legend(energy_names)

        min_ref = min(ref_energies)
        max_ref = max(ref_energies)

        # plot reference line y=0
        x = [min_ref, max_ref]
        y = [0, 0]

        plt.plot(x, y, c='orange', alpha=0.5)

        #Adding axes titles
        plt.xlabel("{} energy [Kcal/mol]".format(ref_energy_name))
        plt.ylabel("delta energy [Kcal/mol]")

        files.init_file(output_file)
        plt.savefig(output_file)
