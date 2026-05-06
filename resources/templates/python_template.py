#!/usr/bin/env python3

"""
The script ...

"""

__author__      = "Christopher Janschke"
__version__     = "1.0"
__email__       = "mbxcj2@nottingham.ac.uk"
__date__        = "xx-xx-xxxx"

#imported modules here



#functions here
def read_file_tsv(path):
    """
    Read in a tsv file and return a pandas dataframe.

    Parameters:
    path (string): Path to the tsv file

    Returns:
    pandas dataframe: tsv file as a pandas dataframe
    """
    try:
        tsv_file = pd.read_csv(path, sep='\t') # read in tsv file from param 'path' as pandas dataframe
        print("File successfully read in.\n")
        return tsv_file # return tsv file as pandas dataframe
    except FileNotFoundError:
        print("File not found. Please enter a valid path.")




def main():
    # Your code here










    pass

if __name__ == "__main__":
    main()