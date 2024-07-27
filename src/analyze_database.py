import os
import sys


def analyze_database(fastq, database, output, rt):
    if rt.lower() != 'illumina' and rt.lower() != 'nanopore':
        sys.exit('Either illumina or nanopore must be given as the read type (rt).')
    # Check if the output directory exists, if not, create it
    if not os.path.exists(output):
        os.makedirs(output)
        print(f"Output directory {output} created.")
    else:
        print(f"Output directory {output} already exists.")

    if rt.lower() == 'illumina':
        fastq.sort()
        for i in range(0, len(fastq), 2):
            if i + 1 < len(fastq):  # Ensure there is a pair
                file_pair = f"{fastq[i]} {fastq[i + 1]}"
                type_stats(file_pair, database, output, rt)
    elif rt.lower() == 'nanopore':
        for file in fastq:
            type_stats(file, database, output, rt)

def type_stats(file, database, output, rt):
    #determine reference species
    if ' ' in file: #assumes a PE string
        single_file = file.split(' ')[0]
        name = os.path.basename(single_file).split('.')[0]
    else:
        name = os.path.basename(file).split('.')[0]
    cmd = f'kma -i {file} -o {output}/{name}_mapping -t_db {database} -mem_mode -Sparse -mf 50000 -ss c -t 4'
    os.system(cmd)

    highest_scoring_template, template_number = highest_scoring_hit(os.path.join(output, f"{name}_mapping.spa"))
    primary_specie = ' '.join(highest_scoring_template.split()[1:3])

    #TBD check these alignments with the stuff from melbourne. Settings could be off.
    #Add nanopore alignment settings
    if rt.lower() == 'illumina':
        cmd = f'kma -i {file} -o {output}/{name}_alignment -t_db {database} -1t1 -mem_mode -Mt1 {template_number} -t 4'
        os.system(cmd)
    elif rt.lower() == 'nanopore':
        cmd = f'kma -i {file} -o {output}/{name}_alignment -t_db {database} -ont -1t1 -mem_mode -Mt1 {template_number} -t 4'
        os.system(cmd)



def highest_scoring_hit(file_path):
    """
    Identifies and returns the highest scoring template from a tab-separated file.

    This function reads through a specified file, assuming it to be tab-separated,
    and identifies the row with the highest score in a designated score column.
    It returns the template corresponding to this highest score.

    Args:
        file_path (str): The path to the file to be read. Assumes a specific format where
                         the score is in the third column and the template in the first.

    Returns:
        str: The identifier of the highest scoring template.
    """

    highest_score = 0
    highest_scoring_template = ""
    template_number = ""

    with open(file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            columns = line.split('\t')
            try:
                # Extract score and compare to find the highest
                score = int(columns[2])  # Score is expected in the 3rd column
                if score > highest_score:
                    highest_score = score
                    highest_scoring_template = columns[0]  # Template is expected in the 1st column
                    template_number = columns[1]
            except ValueError:
                # Skip line if score is not an integer or line is malformed
                continue

    return highest_scoring_template, template_number