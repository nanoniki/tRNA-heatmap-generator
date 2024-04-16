import os
import argparse
import matplotlib
import numpy as np
import matplotlib.pyplot as plt


class HeatmapGenerator:
    def __init__(self, input_file, output_file, lengths_file, adaptors, title, resolution, difference, bidirectional):
        self.input_file = input_file
        self.output_file = output_file
        self.lengths_file = lengths_file
        self.adaptors = [int(adaptors.split(',')[0]), int(adaptors.split(',')[1])]
        self.title = title
        self.difference = difference
        self.bidirectional = bidirectional
        self.resolution = resolution

        self.tRNA_lengths = {}
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
        self.tRNA_max_length = max_length_adjusted

        # If any tRNA length is less than the sum of adaptors, adjust it to 0
        for tRNA, length in self.tRNA_lengths.items():
            if length < sum(self.adaptors):
                self.tRNA_lengths[tRNA] = 0

        # Filter out tRNAs with adjusted length of 0
        self.tRNA_labels = [tRNA for tRNA, length in self.tRNA_lengths.items() if length > 0]

        for name in self.tRNA_labels:
            self.clean_labels.append(name.split("-")[1])

        # Update the tRNA count based on the filtered tRNA labels
        self.tRNA_count = len(self.tRNA_labels)

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
                else:
                    translated_ac = trna[:len(trna)-7] + '-' + anticodon.replace("U", "T") + trna[-4:]
                    
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
            # Check for errors in data length
            if len(posteriors) > self.tRNA_max_length:
                print('Error: The length of', tRNA, 'is longer than max length given. You may not have filtered out '
                                                    'secondary and/or supplementary alignments.')
                print(tRNA, len(posteriors))
            elif len(posteriors) < self.tRNA_max_length:
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
        TICK_FONT_SIZE = 6  # Font size for tick labels
        TITLE_FONT_SIZE = 10  # Font size for the plot title
        LABEL_FONT_SIZE = 8  # Font size for axis labels
        COLOR_BAR_FONT_SIZE = 6  # Font size for color bar labels
        COLOR_BAR_LABEL = 'Reference match probability'  # Label for the color bar

        # Convert centimeters to inches
        cm = 1 / 2.54

        # Create a new figure and axis object
        fig, ax = plt.subplots(figsize=(FIG_WIDTH_CM * cm, FIG_HEIGHT_CM * cm))

        # Define the color map
        if self.difference:
            color_map = matplotlib.colors.LinearSegmentedColormap.from_list("", ["lightgoldenrodyellow", "purple",
                                                                                 "midnightblue"])
        else:
            color_map = 'YlGnBu'

        # Plot the heatmap
        if self.bidirectional:
            im = ax.imshow(pubs, cmap=color_map, vmin=-0.5, vmax=0.5, aspect=1)
        else:
            im = ax.imshow(pubs, cmap=color_map, vmin=0, vmax=1, aspect=1)

        # Set tick marks and labels for the x-axis
        ax.set_xticks([i for i in range(0, (self.tRNA_max_length - 1), TICK_STEP)],
                      labels=[i + 1 for i in range(0, (self.tRNA_max_length - 1), TICK_STEP)])

        # Set tick marks and labels for the y-axis
        ax.set_yticks(np.arange(len(noT_labels)), labels=noT_labels)

        # Set font size and alignment for x-axis tick labels
        plt.setp(ax.get_xticklabels(), fontsize=TICK_FONT_SIZE, ha="center", rotation_mode="anchor")

        # Set font size for y-axis tick labels
        plt.setp(ax.get_yticklabels(), fontsize=TICK_FONT_SIZE)

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
                                                 "generated with marginCaller (see https://github.com/benedictpaten/marginAlign-tRNA) "
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
                                             args.resolution, args.difference, args.bidirectional)
        heatmap_generator.load_lengths()
        heatmap_generator.import_data()

        noT_labels = [heatmap_generator.TtoU(i) for i in heatmap_generator.clean_labels]

        positions, posterior_prob = heatmap_generator.structure_data()

        # Process the second input file
        second_heatmap_generator = HeatmapGenerator(input_file2, args.output, args.lengths, args.adaptors,
                                                    args.title, args.resolution, args.difference, args.bidirectional)
        second_heatmap_generator.load_lengths()
        second_heatmap_generator.import_data()

        # Get posterior probability matrix for the second file
        _, second_posterior_prob = second_heatmap_generator.structure_data()

        # Calculate the difference matrix
        difference_matrix = second_posterior_prob - posterior_prob

        # Plot the heatmap with the difference matrix and save it with the modified output file path
        heatmap_generator.plot_heatmap(difference_matrix, noT_labels)
    else:
        # Process the input file without calculating the difference
        heatmap_generator = HeatmapGenerator(input_file1, args.output, args.lengths, args.adaptors, args.title,
                                             args.resolution, args.difference, args.bidirectional)
        heatmap_generator.load_lengths()
        heatmap_generator.import_data()

        noT_labels = [heatmap_generator.TtoU(i) for i in heatmap_generator.clean_labels]

        positions, posterior_prob = heatmap_generator.structure_data()
        heatmap_generator.plot_heatmap(posterior_prob, noT_labels)


if __name__ == "__main__":
    main()
