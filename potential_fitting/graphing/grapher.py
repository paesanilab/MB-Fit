import matplotlib.pyplot as plt
from collections import OrderedDict
import csv
import math

class Grapher:

    def make_distribution_graph(self, training_sets, training_set_names, energy_name, output_file, num_bins=50, max_percent=0.97):

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

        plt.savefig(output_file)

    def make_correlation_graph(self, training_set, energy_name1, energy_name2, output_file):
        pass

    def make_error_graph(self, training_set, energy_name1, energy_name2, output_file):
        pass
