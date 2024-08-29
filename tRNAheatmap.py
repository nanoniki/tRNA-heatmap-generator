import os
import argparse
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import ast
import re

import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator

plt.rcParams['pdf.fonttype'] = 42
plt.switch_backend('agg')


class HeatmapGenerator:
    def __init__(self, input_file, output_file, lengths_file, adaptors, title, resolution, difference, bidirectional,
                 structure):
        self.input_file = input_file
        self.output_file = output_file
        self.lengths_file = lengths_file
        self.structure_file = structure
        self.adaptors = [int(adaptors.split(',')[0]), int(adaptors.split(',')[1])]
        self.title = title
        self.difference = difference
        self.bidirectional = bidirectional
        self.resolution = resolution

        self.tRNA_lengths = {}
        self.tRNA_structure = {}
        self.tRNA_labels = []  # labels for y-axis
        self.clean_labels = []
        self.tRNA_count = 0  # the number of tRNAs
        self.tRNA_max_length = 0  # the max tRNA length w/o adapters
        self.data = []

    def load_lengths(self):
        """Load tRNA lengths from a file and initialize relevant attributes. This method reads a tab-separated file
        containing tRNA names and their respective lengths. It parses each line of the file, extracting the tRNA name
        and its length, and stores them in the `tRNA_lengths` dictionary. After loading the lengths, it populates the
        `tRNA_labels` list attribute with the tRNA names and sets the `tRNA_count` attribute to the number of tRNAs
        loaded. It also calculates the maximum length among the tRNAs and subtracts the lengths of any specified
        adaptors, storing the result in the `tRNA_max_length` attribute."""

        with open(self.lengths_file, 'r') as f:
            for line in f:
                tRNA, length = line.strip().split()
                self.tRNA_lengths[tRNA] = int(length)

        # Calculate the maximum length of tRNAs after considering adaptors
        # Adjust the calculation to properly account for adaptors
        max_length_adjusted = max(self.tRNA_lengths.values()) - sum(self.adaptors)

        # Update the maximum length attribute
        if self.structure_file:
            self.tRNA_max_length = max_length_adjusted + 3
        else:
            self.tRNA_max_length = max_length_adjusted

        # If any tRNA length is less than the sum of adaptors, adjust it to 0
        for tRNA, length in self.tRNA_lengths.items():
            if length < sum(self.adaptors):
                self.tRNA_lengths[tRNA] = 0

        # Filter out tRNAs with adjusted length of 0
        self.tRNA_labels = [tRNA for tRNA, length in self.tRNA_lengths.items() if length > 0]

        for name in self.tRNA_labels:
            if 'Tyr-GTA' in name:
                self.clean_labels.append(name.split("-")[1]+name.split("-")[2])
            else:
                self.clean_labels.append(name.split("-")[1])

        # Update the tRNA count based on the filtered tRNA labels
        self.tRNA_count = len(self.tRNA_labels)

    def load_structure(self):

        with open(self.structure_file, 'r') as f:
            # Skip the header
            next(f)
            for line in f:
                tRNA = line.strip().split('\t')[0]
                if 'Tyr' in tRNA:
                    tRNA = "Saccharomyces_cerevisiae_tRNA-Tyr-GTA-1-1"
                v_loop = int(line.strip().split('\t')[-1])
                elements = line.strip().split('\t')[1:-1]
                other_v = []
                for e in elements:
                    if e.startswith('(') and e.endswith(')'):
                        other_v.append(ast.literal_eval(e))

                coordinates = other_v + [v_loop]

                self.tRNA_structure[tRNA] = coordinates

    def UtoT(self, trna):
        """This function takes a tRNA sequence as input and replaces 'U' nucleotides in the anticodon
        region with 'T'. It specifically targets the last three characters of the input sequence
        as the anticodon. If the input sequence has fewer than 3 characters, it returns None."""
        if len(trna) >= 3:
            if not re.search(r'\d', trna[-3:]):
                anticodon = trna[-3:]
                # join amino acid identity with reverse translated anti-codon
                translated_ac = trna[:-3] + anticodon.replace("U", "T")
            else:
                anticodon = trna[-7:].split('-')[0]
                if 'iMet' in trna:
                    translated_ac = trna[:len(trna) - 7] + anticodon.replace("U", "T") + trna[-4:]
                elif 'Tyr' in trna:
                    translated_ac = trna[:len(trna) - 7] + anticodon.replace("U", "T") + trna[-4:]
                else:
                    translated_ac = trna[:len(trna) - 7] + '-' + anticodon.replace("U", "T") + trna[-4:]

            return translated_ac

    def TtoU(self, lbl):
        """This function performs transcription by replacing all 'T' characters with 'U'.
        Additionally, it handles specific codon mappings where 'T' is replaced by 'U'
        according to the provided codon_mapping dictionary. If the given label does not
        match any specific codon, it replaces 'T' with 'U' in the label."""

        codon_mapping = {
            'ThrAGT': 'ThrAGU',
            'ThrCGT': 'ThrCGU',
            'ThrTGT': 'ThrUGU',
            'TrpCCA': 'TrpCCA',
            'TyrGTA': 'TyrGUA',
            'iMet': 'iMetCAU'}

        return codon_mapping.get(lbl, lbl.replace("T", "U"))

    def import_data(self):
        """Import and process data from the marginCaller generated input file. This function reads data from the
        specified input file, parses each line, and extracts relevant information based on predefined conditions.
        It populates the self.data attribute with processed data for further analysis."""

        # Open the input file for reading
        with open(self.input_file, 'r') as file:
            # Iterate over each line in the file
            for line in file:
                # Skip lines starting with '#'
                if line.startswith('#'):
                    continue

                # Skip lines with 'mito' in header
                if 'mito' in line:
                    continue

                # Split the line by tabs
                l = line.strip().split('\t')

                # Extract data from the line
                header = l[0]
                for trna, length in self.tRNA_lengths.items():
                    tRNA = self.UtoT(trna)
                    if tRNA in header:
                        new_header = trna if trna != 'MetCAT' or 'iMetCAT' in header else trna

                tRNA_len = self.tRNA_lengths[new_header]
                position = int(l[1])
                reference_nt = l[3]
                alternative_nt = l[4].split(',')
                posterior_probs = [float(x) for x in l[7].split(',')]

                if self.adaptors[0] < position <= tRNA_len - self.adaptors[1]:
                    line_data = [new_header, (position - self.adaptors[0]), reference_nt, alternative_nt,
                                 posterior_probs]
                    self.data.append(line_data)
                else:
                    pass

    def structure_data(self):
        """
        Structures the input data into position, nucleotide, and posterior probability matrices for each tRNA.

        This function initializes dictionaries to store position, nucleotide, and posterior probability data
        for each tRNA. It then iterates over each data entry, calculates the reference posterior probability,
        and populates the dictionaries accordingly. Next, it initializes a posterior probability matrix with zeros
        and fills it with data from the posterior probability dictionary. If a tRNA has more posterior probability
        values than the maximum length allowed, it prints an error message. If a tRNA has fewer values, it pads
        the list with None values to match the maximum length. Finally, it returns a list of positions and the
        posterior probability matrix.

        Returns:
        - Tuple[List[int], numpy.ndarray]: A tuple containing a list of positions and a 2D numpy array
          representing the posterior probability matrix for all tRNAs.

        Note:
        - The input data should be a list of tuples, where each tuple contains:
          (tRNA_label, position, nucleotide, posterior_probabilities)."""

        # Initialize dictionaries to store position, nucleotide, and posterior probability data for each tRNA
        pos_dict = {tRNA: [] for tRNA in self.tRNA_labels}
        nt_dict = {tRNA: [] for tRNA in self.tRNA_labels}
        pp_dict = {tRNA: [] for tRNA in self.tRNA_labels}

        # Iterate over each data entry and populate dictionaries
        for entry in self.data:
            # Calculate the reference posterior probability
            ref_pp = round(1.0 - sum(entry[4]), 3)

            # Append data to respective dictionaries
            pos_dict[entry[0]].append(entry[1])
            nt_dict[entry[0]].append(entry[2])
            pp_dict[entry[0]].append(ref_pp)

        # Initialize posterior probability matrix with zeros
        posterior_prob = np.zeros((self.tRNA_count, self.tRNA_max_length))

        # Fill in posterior probability matrix with data from pp_dict

        for i, tRNA in enumerate(self.tRNA_labels):
            posteriors = pp_dict[tRNA]

            # print(i, tRNA, self.tRNA_structure.keys())
            nan_tuples = self.tRNA_structure[tRNA][:-1]
            # print(tRNA, nan_tuples)

            for coord, count in nan_tuples:
                # print(tRNA, coord, count)
                nan_values = [np.nan for _ in range(count)]
                posteriors[coord:coord] = nan_values
                # posteriors.insert(coord, (np.nan * count))

            tRNA_len_noA = pos_dict[tRNA][-1] + len(self.tRNA_structure[tRNA][:-1])  # tRNA length no adaptor
            tmp_vll = abs(self.tRNA_max_length - tRNA_len_noA) - (sum(y for _, y in nan_tuples) - len(nan_tuples))
            if tmp_vll < 0:
                var_loop_len = 0
            else:
                var_loop_len = abs(self.tRNA_max_length - tRNA_len_noA) - (
                            sum(y for _, y in nan_tuples) - len(nan_tuples))
            var_loop_start = (self.tRNA_structure[tRNA][-1] - self.adaptors[0] - 1)  # variable loop start position

            # Check for errors in data length
            if len(posteriors) > self.tRNA_max_length:
                print('Error: The length of', tRNA, 'is longer than max length given. You may not have filtered out '
                                                    'secondary and/or supplementary alignments.')
                # check TyrGTA GluCTC or iMet labels and dictionary keys
                print('tRNA, length:', tRNA, len(posteriors))
                print('posteriors:', posteriors)
                print('max length:', self.tRNA_max_length)

            elif len(posteriors) < self.tRNA_max_length:

                if self.structure_file:
                    # Insert None for each position of the variable loop
                    posteriors[var_loop_start:var_loop_start] = [np.nan] * var_loop_len
                else:
                    # Pad the list with None values to match the max length
                    posteriors.extend([None] * (self.tRNA_max_length - len(posteriors)))

            # Fill in posterior probability matrix row
            posterior_prob[i, :] = posteriors[::-1]  # Reverse the list to match positions

        return list(range(1, self.tRNA_max_length + 1)), np.fliplr(posterior_prob)

    def plot_heatmap(self, pubs, noT_labels):
        """This function plots a heatmap of publication data using matplotlib.

        Parameters:
        - publication_data (numpy.ndarray): A 2D array representing publication data.
        - noT_labels (list): A list of labels for tRNAs."""

        # Constants
        FIG_WIDTH_CM = 18  # Width of the figure in centimeters
        FIG_HEIGHT_CM = 12  # Height of the figure in centimeters
        TICK_STEP = 4  # Step for the tick marks on the x-axis
        if self.structure_file:
            X_TICK_FONT_SIZE = 4  # Font size for xtick labels in structured heatmap
        else:
            X_TICK_FONT_SIZE = 6  # Font size for xtick labels in regular heatmap
        Y_TICK_FONT_SIZE = 4  # Font size for ytick labels
        TITLE_FONT_SIZE = 10  # Font size for the plot title
        LABEL_FONT_SIZE = 6  # Font size for axis labels
        COLOR_BAR_FONT_SIZE = 6  # Font size for color bar labels

        if self.difference:
            COLOR_BAR_LABEL = 'Change in reference match probability'  # Label for the color bar
        else:
            COLOR_BAR_LABEL = 'Reference match probability'  # Label for the color bar

        yeast_structure_xticks = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14',
                                  '15', '16', '17', '18', '19', '20', '20a ', '20b ', '21', '22', '23', '24', '25', '26',
                                  '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41',
                                  '42', '43', '44', '45', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', '46',
                                  '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61',
                                  '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76']

        # Convert centimeters to inches
        cm = 1 / 2.54

        # Create a new figure and axis object
        fig, ax = plt.subplots(figsize=(FIG_WIDTH_CM * cm, FIG_HEIGHT_CM * cm))

        # Define the color map
        if self.difference:
            if self.bidirectional:
                color_map = 'Spectral_r'
            else:
                color_map = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightgoldenrodyellow", "purple",
                                                                                     "midnightblue"])
        else:
            if self.bidirectional:
                color_map = 'Spectral_r'
            else:
                color_map = 'YlGnBu'

        # Plot the heatmap
        if self.bidirectional:
            im = ax.imshow(pubs, cmap=color_map, vmin=-0.5, vmax=0.5, aspect=1)
        else:
            im = ax.imshow(pubs, cmap=color_map, vmin=0, vmax=1, aspect=1)

        # Set tick marks and labels for the x-axis
        # Set font size and alignment for x-axis tick labels
        if self.structure_file:
            ax.set_xticks(range(self.tRNA_max_length))
            ax.set_xticklabels(yeast_structure_xticks[:self.tRNA_max_length])
            plt.setp(ax.get_xticklabels(), fontsize=X_TICK_FONT_SIZE, ha="center", va="center", rotation_mode="anchor",
                     rotation=90)

            # Overlay for None values
            for (i, j), val in np.ndenumerate(pubs):
                if np.isnan(val):
                    # Calculate the center y-coordinate of the cell
                    y_center = i
                    x_center = j
                    # Select the marker based on the row index
                    marker_clr = 'black'

                    # Draw a dot in the middle of the cell only if marker_clr is not None
                    if marker_clr is not None:
                        ax.scatter(x_center, y_center, color=marker_clr, s=0.1,
                                   marker='o')  # 's' is the size of the dot

        else:
            ax.set_xticks([i for i in range(0, (self.tRNA_max_length - 1), TICK_STEP)],
                          labels=[i + 1 for i in range(0, (self.tRNA_max_length - 1), TICK_STEP)])
            plt.setp(ax.get_xticklabels(), fontsize=X_TICK_FONT_SIZE, ha="center", rotation_mode="anchor")

        # Set tick marks and labels for the y-axis
        ax.set_yticks(np.arange(len(noT_labels)), labels=noT_labels)

        # Set font size for y-axis tick labels
        plt.setp(ax.get_yticklabels(), fontsize=Y_TICK_FONT_SIZE)

        # Set plot title
        ax.set_title(self.title, fontsize=TITLE_FONT_SIZE, pad=7)

        # Add color bar
        cbar = plt.colorbar(im, shrink=0.5)
        cbar.ax.tick_params(labelsize=COLOR_BAR_FONT_SIZE)
        cbar.ax.set_ylabel(COLOR_BAR_LABEL, fontsize=LABEL_FONT_SIZE, labelpad=10, rotation=270)

        # Set axis labels
        ax.set_xlabel("Nucleotide Position (5'→3')", fontsize=LABEL_FONT_SIZE)
        ax.set_ylabel("tRNA Isoacceptor", fontsize=LABEL_FONT_SIZE)

        # Save the plot to a file
        plt.savefig(self.output_file, transparent=True, dpi=self.resolution)


def main():
    parser = argparse.ArgumentParser(description="This code creates a heatmap of posterior probabilities along tRNA "
                                                 "isoacceptors sequenced with Nanopore. The input files must be "
                                                 "generated with marginCaller (see "
                                                 "https://github.com/benedictpaten/marginAlign-tRNA)"
                                                 "using the setting --threshold=0 so that every position is output.")
    parser.add_argument("-i", "--input", nargs='+', required=True, metavar="INPUT",
                        help="One or two input pileup .vcf files generated with marginCaller. If --difference is used, "
                             "provide two files.")
    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Path to output location with desired file extension (pdf, png, jpeg, etc).")
    parser.add_argument("-l", "--lengths", required=True, type=str,
                        help="Path to a tab separated file containing the name of each tRNA and their respective "
                             "lengths with the adaptors. The naming conventions used for tRNA in this file will "
                             "appear as y-axis labels.")
    parser.add_argument("-a", "--adaptors", required=True, type=str,
                        help="The lengths of your 5' and 3' adaptors separated by a comma.", metavar="ADAPTORS")
    parser.add_argument("-t", "--title", type=str, default=None, help="Title of plot. Default is no title.")
    parser.add_argument("-d", "--difference", required=False, action="store_true",
                        help="Calculate the difference between two input files (second input - first input).")
    parser.add_argument("-b", "--bidirectional", required=False, action="store_true",
                        help="Produces difference plot with directionality maintained (cbar normalized from [-0.5, "
                             "0.5]).")
    parser.add_argument("-s", "--structure", required=False, type=str, help="Path to a tab separated file containing "
                                                                            "the name of each tRNA and the respective "
                                                                            "coordinates of their variable positions.")
    parser.add_argument("-r", "--resolution", default=300, required=False, type=int,
                        help="Resolution of final image in dpi.")

    args = parser.parse_args()

    if args.difference and len(args.input) != 2:
        parser.error("When using --difference, exactly two input files are required.")

    if args.difference:
        input_file1, input_file2 = args.input
    else:
        input_file1 = args.input[0]

    if args.difference:
        # Split the filename into base name and extension
        base_name, extension = os.path.splitext(args.output)

        # Modify the output file path for the difference heatmap
        args.output = f"{base_name}_difference{extension}"

        # Process the first input file
        heatmap_generator = HeatmapGenerator(input_file1, args.output, args.lengths, args.adaptors, args.title,
                                             args.resolution, args.difference, args.bidirectional, args.structure)
        heatmap_generator.load_lengths()
        heatmap_generator.load_structure()
        heatmap_generator.import_data()

        noT_labels = [heatmap_generator.TtoU(i) for i in heatmap_generator.clean_labels]

        positions, posterior_prob = heatmap_generator.structure_data()

        # Process the second input file
        second_heatmap_generator = HeatmapGenerator(input_file2, args.output, args.lengths, args.adaptors,
                                                    args.title, args.resolution, args.difference, args.bidirectional,
                                                    args.structure)
        second_heatmap_generator.load_lengths()
        heatmap_generator.load_structure()
        second_heatmap_generator.import_data()

        # Call load_structure() first to ensure the structure is loaded or refreshed
        second_heatmap_generator.load_structure()

        # Get posterior probability matrix for the second file
        _, second_posterior_prob = second_heatmap_generator.structure_data()

        # Calculate the difference matrix
        difference_matrix = second_posterior_prob - posterior_prob

        # Plot the heatmap with the difference matrix and save it with the modified output file path
        heatmap_generator.plot_heatmap(difference_matrix, noT_labels)
    else:
        # Process the input file without calculating the difference
        heatmap_generator = HeatmapGenerator(input_file1, args.output, args.lengths, args.adaptors, args.title,
                                             args.resolution, args.difference, args.bidirectional, args.structure)
        heatmap_generator.load_lengths()
        heatmap_generator.load_structure()
        heatmap_generator.import_data()

        noT_labels = [heatmap_generator.TtoU(i) for i in heatmap_generator.clean_labels]

        positions, posterior_prob = heatmap_generator.structure_data()
        heatmap_generator.plot_heatmap(posterior_prob, noT_labels)


if __name__ == "__main__":
    main()
