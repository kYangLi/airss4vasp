#!/usr/bin/python3.7
#

import os
import fnmatch
import math

cell_data_file = ''
for filename in os.listdir('.'):
  if fnmatch.fnmatch(filename, '*-out.cell'):
    cell_data_file = filename
    del filename
if '' == cell_data_file:
  print("ERROR: No vaild *.cell file available!!!")
  exit(0)

with open(cell_data_file, 'r') as fpr:
  str_file_lines = fpr.readlines()

while True:
  