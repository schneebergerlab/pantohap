# <editor-fold desc='Step1: Generate haplotype graph'>

# Currently use the graph generated using binning of reference genomes

# </editor-fold>


# <editor-fold desc='Kmer counting'>

# for each 100kb window, get syntenic sequence from query (40 hap) genomes.
get_node_query_sequence()

# Get kmers for each of the fasta file created using get_node_query_sequence()
get_unique_kmer_per_window.sh

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/


# </editor-fold>
