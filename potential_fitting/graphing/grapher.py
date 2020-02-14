import matplotlib.pyplot as plt
from collections import OrderedDict
import csv
import math

from potential_fitting.utils import files

class Grapher:

    def make_distribution_graph(self, training_sets, training_set_names, energy_name, output_file, num_bins=50, max_percent=0.97):

        plt.clf()

        energies_dict = OrderedDict()

        for training_set_name, training_set in zip(training_set_names, training_sets):
            energies_dict[training_set_name] = sorted(training_set.get_energies(energy_name))

        minimum = min([energies[0] for energies in energies_dict.values()])
        try:
            max_percent_energy = max([energies[math.ceil(len(energies) * max_percent)] for energies in energies_dict.values()])
        except IndexError:
            max_percent_energy = max([energies[-1] for energies in energies_dict.values()])

        maximum = max([energies[-1] for energies in energies_dict.values()])

        try:
            bin_width = math.ceil((max_percent_energy - minimum) / num_bins)
        except ZeroDivisionError:
            bin_width = 1

        try:
            min_bin = math.floor(minimum / bin_width) * bin_width # inclusive
        except ZeroDivisionError:
            min_bin = 0

        try:
            max_bin = math.ceil(max_percent_energy / bin_width) * bin_width  # exclusive
        except ZeroDivisionError:
            max_bin = 1

        try:
            overflow_max_bin = math.ceil(maximum / bin_width) * bin_width  # exclusive
        except ZeroDivisionError:
            overflow_max_bin = 1

        # create the bins between min_bin and max_bin

        bins = list(range(min_bin, max_bin + bin_width, bin_width))

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

        tick_width = math.ceil(bin_width * num_bins / num_ticks)

        ticks = list(range(min_bin, max_bin + bin_width, tick_width))

        labels = [str(tick) for tick in ticks]

        ticks.append(max_bin)

        labels.append(str(max_bin) + "+")

        plt.xticks(ticks, labels)

        files.init_file(output_file)
        plt.savefig(output_file)

    def make_correlation_graph(self, training_set, ref_energy_name, energy_names, output_file):

        plt.clf()

        ref_energies = training_set.get_energies(ref_energy_name)

        energies_dict = OrderedDict()

        for energy_name in energy_names:
            energies_dict[energy_name] = training_set.get_energies(energy_name)

        colors = [(0, 0, 0, 1), (0, 0, 1, 0.7), (0, 1, 0, 0.5), (1, 0, 0, 0.3)]

        for index, energies in enumerate(energies_dict.values()):
            plt.scatter(ref_energies, energies,
                        color=colors[index], s=5, alpha=0.5)
        # Adding a legend
        plt.legend(energy_names)

        plt.plot(ref_energies, ref_energies, c='orange', alpha=0.5)


        #Adding axes titles
        plt.xlabel("{} energy [Kcal/mol]".format(ref_energy_name))
        plt.ylabel("energy [Kcal/mol]")

        files.init_file("fit.png")

        files.init_file(output_file)
        plt.savefig(output_file)


    def make_error_graph(self, training_set, ref_energy_name, energy_names, output_file):

        plt.clf()

        ref_energies = training_set.get_energies(ref_energy_name)

        energies_dict = OrderedDict()

        for energy_name in energy_names:
            energies_dict[energy_name] = training_set.get_energies(energy_name)

        colors = [(0, 0, 0, 1), (0, 0, 1, 0.7), (0, 1, 0, 0.5), (1, 0, 0, 0.3)]

        for index, energies in enumerate(energies_dict.values()):
            plt.scatter(ref_energies, [energy - ref_energy for energy, ref_energy in zip(energies, ref_energies)],
                        color=colors[index], s=5, alpha=0.5)
        # Adding a legend
        plt.legend(energy_names)

        plt.plot(ref_energies, ref_energies, c='orange', alpha=0.5)


        #Adding axes titles
        plt.xlabel("{} energy [Kcal/mol]".format(ref_energy_name))
        plt.ylabel("delta energy [Kcal/mol]")

        files.init_file("fit.png")

        files.init_file(output_file)
        plt.savefig(output_file)
