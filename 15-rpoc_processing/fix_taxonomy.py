''' Adds missing fields to list of taxonomy strings. Use: python fix_taxonomy.py my_taxa.txt > new_taxa.txt '''

import sys

with open(sys.argv[1], "r") as f:
	for line in f:
		x = line.rstrip("\n")
		y = x.count(";")
		if y == 0:
			new = (x.split(";")[0] + ";") + ("unclassified;") * 10 + "unclassified"
			print(new)
		elif y < 11:
			add = 11 - y
			new = (";" + x.split(";")[-1]) * add
			print(x + new)
		else:
			print(x)