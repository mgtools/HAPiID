"""
script used to change all the genome extensions to '.fasta'
"""

import os
import sys

if len(sys.argv) != 3:
   print('to run this script please enter 2 command line argument, exampel to run script i.e. python3 unifyGenomeExtensions.py /dir/to/genome/files .extension')

else:
   in_dir = sys.argv[1]
   extension = sys.argv[2]
   for file in os.listdir(in_dir):
       out_file = file.rsplit('.', 1)[0]
       os.rename(in_dir+file, in_dir+out_file+extension)

