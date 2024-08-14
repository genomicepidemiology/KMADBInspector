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

    species_reference_count = 0
    if rt.lower() == 'illumina':
        fastq.sort()
        for i in range(0, len(fastq), 2):
            if i + 1 < len(fastq):  # Ensure there is a pair
                file_pair = f"{fastq[i]} {fastq[i + 1]}"
                count = type_stats(file_pair, database, output, rt)
                species_reference_count += count
        species_reference_count = species_reference_count/(len(fastq)/2)
    elif rt.lower() == 'nanopore':
        for file in fastq:
            count = type_stats(file, database, output, rt)
            species_reference_count += count
        species_reference_count = species_reference_count / len(fastq)

    print(f'Average number of species references: {species_reference_count}')

    res_files = [os.path.join(output, f) for f in os.listdir(output) if f.endswith('.res')]

    # Initialize lists to store the highest identities and coverages
    highest_identities = []
    highest_coverages = []

    # Process each result file
    for result_file in res_files:
        highest_identity, highest_coverage = get_highest_template_identity_and_coverage(result_file)
        highest_identities.append(highest_identity)
        highest_coverages.append(highest_coverage)

    # Calculate and print the average highest template identity
    average_highest_identity = sum(highest_identities) / len(highest_identities) if highest_identities else 0.0
    print(f'Average highest Template Identity: {average_highest_identity}')

    # Calculate and print the average highest query coverage
    average_highest_coverage = sum(highest_coverages) / len(highest_coverages) if highest_coverages else 0.0
    print(f'Average highest Query Coverage: {average_highest_coverage}')

    # K-mer counting
    kmer_count = count_unique_kmers(fastq, rt)
    print(f'Total unique k-mers: {kmer_count}')


def type_stats(file, database, output, rt):
    # Determine reference species
    if ' ' in file:  # Assumes a PE string
        single_file = file.split(' ')[0]
        name = os.path.basename(single_file).split('.')[0]
    else:
        name = os.path.basename(file).split('.')[0]
    cmd = f'kma -i {file} -o {output}/{name}_mapping -t_db {database} -mem_mode -Sparse -mf 50000 -ss c -t 4'
    os.system(cmd)

    highest_scoring_template, template_number = highest_scoring_hit(os.path.join(output, f"{name}_mapping.spa"))
    primary_specie = ' '.join(highest_scoring_template.split()[1:3])

    # TBD: Check these alignments with the settings from Melbourne. Settings could be off.
    # Add nanopore alignment settings. Do we need to perhaps fragment the reads?
    if rt.lower() == 'illumina':
        cmd = f'kma -i {file} -o {output}/{name}_alignment -t_db {database} -1t1 -mem_mode -Mt1 {template_number} -t 4'
        os.system(cmd)
    elif rt.lower() == 'nanopore':
        cmd = f'kma -i {file} -o {output}/{name}_alignment -t_db {database} -ont -1t1 -mem_mode -Mt1 {template_number} -t 4'
        os.system(cmd)

    count = count_species_references(database + '.name', primary_specie)

    return count


def get_highest_template_identity_and_coverage(result_file_path):
    highest_identity = 0.0  # Initialize the highest template identity as zero
    highest_coverage = 0.0  # Initialize the highest query coverage as zero
    with open(result_file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            columns = line.strip().split('\t')
            template_identity = float(columns[4])  # Template_Identity is the 5th column (index 4)
            query_coverage = float(columns[5])  # Query_Coverage is the 6th column (index 5)

            # Update highest identity if the current value is greater
            if template_identity > highest_identity:
                highest_identity = template_identity

            # Update highest coverage if the current value is greater
            if query_coverage > highest_coverage:
                highest_coverage = query_coverage

    return highest_identity, highest_coverage


def highest_scoring_hit(file_path):
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


def count_species_references(file_path, species_name):
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if species_name in line:
                count += 1
    return count


def count_unique_kmers(fastq_files, rt):
    if rt.lower() == 'illumina':
        input_files = ' '.join(fastq_files)  # Join both files if paired-end
    else:
        input_files = fastq_files[0]  # Use the single file for Nanopore

    # Jellyfish is a popular tool for k-mer counting
    cmd = f'jellyfish count -m 21 -s 100M -t 4 -C {input_files} -o mer_counts.jf'
    os.system(cmd)

    # Dump the counts into a text file
    cmd = 'jellyfish dump mer_counts.jf > mer_counts.txt'
    os.system(cmd)

    # Count the number of unique k-mers
    with open('mer_counts.txt', 'r') as file:
        unique_kmers = sum(1 for line in file if line.strip())  # Each line represents a unique k-mer

    # Clean up intermediate files
    os.remove('mer_counts.jf')
    os.remove('mer_counts.txt')

    return unique_kmers


# Example usage:
# fastq_files = ['sample1_R1.fastq', 'sample1_R2.fastq']  # For Illumina
# database = 'database_name'
# output_directory = './output'
# read_type = 'illumina'  # or 'nanopore'

# analyze_database(fastq_files, database, output_directory, read_type)
