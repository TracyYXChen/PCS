#!/usr/bin/python
from pymol.cgo import *    
from pymol import cmd
from glob import glob

cmd.load('D_aligned_20000/1d3z.pdb','a_1d3z')
cmd.load('D_aligned_20000/2_1d3z.pdb','b_1d3z')

curr_dir = os.getcwd()
i = 1
for files in glob("D_aligned_20000/mnodes*.pdb"):
	cmd.load(files,"mnodes")
	cmd.align('a_1d3z','/mnodes//c')
	cmd.align('b_1d3z','/mnodes//d')
	cmd.save('myalign/ali_' + str(i) +'.pdb')
	cmd.delete(all)
	i += 1

