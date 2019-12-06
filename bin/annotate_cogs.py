#!/usr/bin/env python

import sys
from collections import Counter

#func_categories_file = '../data/functional_categories.txt'

cog_file = sys.argv[1] # file with OGs which fulfilled a condition
output_file = sys.argv[2] #file to write results to
func_categories_file = sys.argv[3] #file containing functional categories descriptions

func_file = sys.stdin #file with annotations for each OG

# Initialize variables
category2description = {} #dict
cog2category = {} #dict
category2count = Counter() #counter
my_cogs = set() #set


# Step 1: Read description of functional categories
with open(func_categories_file) as fin:
    for line in fin:
        line = line.strip()
        if line == "":
            continue
        if line[0] != "[":
            continue
        line = line.split(' ')
        #get category letter
        cat = line[0][1]
        
        x = len(line)
        #get description of category
        descr = " ".join(line[1:x])
        
        #save in dictionary
        category2description[cat] = descr
          
    

# Step 2: Read file with OGs of interest
with open(cog_file) as fin:
    for line in fin:
        line = line.strip()
        if line.startswith("#"):
            continue
        line = line.split("\t")
        my_cogs.add(line[0])


# Step 3: Read file with functional annotations
#with open(func_file) as fin:
for line in func_file:
    line = line.strip().split("\t")
    
    cog = line[1]
    cat = line[2]
    
    cog2category[cog] = cat
    

# Step 4: Count and output the categories for OGs of interest

for cog in my_cogs:
    category2count.update(cog2category[cog])

# Step 5: Write counted categories to file

with open(output_file, 'w') as fout:
    for cat in category2count:
        output = "\t".join([str(category2count[cat]), category2description[cat]])
        print(output, file=fout)