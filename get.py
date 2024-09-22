import os 
import zipfile
import json
import subprocess

# Extract the genomic data from an archive file
# Input: file - the archive file containing the genomic data
# Output: a string containing the genomic data in FASTA format
def extract(file: str) -> str:
    if not file.endswith('.zip'):
        raise ValueError('invalid file type')
    
    print(f'Extracting {file}...')
    # Create a directory to store the extracted files
    extract_dir = file.split('.')[0]

    if not os.path.exists(extract_dir):
        os.makedirs(extract_dir, exist_ok=True)

        # Extract the files
        with zipfile.ZipFile(file, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        
        print('Extraction complete.')
    else:
        print('Files already extracted.')

    data_dir = f'{extract_dir}/ncbi_dataset/data'

    # open dir/dataset_catalog.json
    with open(f'{data_dir}/dataset_catalog.json') as f:
        data = json.load(f)
    
    asm = data['assemblies']

    genome = None

    for a in asm:
        # checking for the file containing the nucleotide data
        for f in a['files']:
            if f['fileType'] == 'GENOMIC_NUCLEOTIDE_FASTA':
                # open dir/f['fileType']
                print(f'Reading {f["filePath"]}...')
                with open(f'{data_dir}/{f["filePath"]}') as g:
                    # read the file
                    # remove the first line
                    # remove the newline characters
                    genome = g.readlines()[1:]
                    genome = ''.join(genome).replace('\n', '')
                    # make all characters uppercase
                    genome = genome.upper()

    return genome


# This file is not really a python file, its more like a wrapper for the ncbi CLI.
# Python format for convenience.

def download(taxon):

    file = f'data/{taxon.replace(" ","-")}.zip'

    if os.path.exists(file):
        return file

    cmd = f'ncbi-datasets download genome taxon "{taxon}" --filename {file} --reference'

    os.system(cmd)

    return file


def get_metadata(taxon):
    # Execute the shell command
    result = subprocess.run(
        ['ncbi-datasets', 'summary', 'taxonomy', 'taxon', taxon],
        capture_output=True,
        text=True
    )
    
    # Check if the command was successful
    if result.returncode != 0:
        raise Exception(f"Command failed with error: {result.stderr}")
    
    try :
        # Parse the JSON output
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        return None
    
    return data


def get_fasta(taxon):

    file = download(taxon)

    return extract(file)

