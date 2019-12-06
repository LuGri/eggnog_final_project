###### Use the eggNOG5.0 database to identify conserved single-copy genes.

#create directories if they don't exist for whatever reason
#shouldn't really happen because the whole project structure is intended to be cloned from github
mkdir -p bin data doc results

#set permissions
cd ./bin
chmod +x annotate_cogs.py
chmod +x read_members_file.py
cd ..

#get data from FH server
cd ./data

DATADIR=/mirror/eggnog/eggnog_5.0/per_tax_level/

#members files symbolic links
ln -s $DATADIR/2/2_members.tsv.gz 2_members.tsv.gz
ln -s $DATADIR/2157/2157_members.tsv.gz 2157_members.tsv.gz
ln -s $DATADIR/2759/2759_members.tsv.gz 2759_members.tsv.gz

#annotations file symbolic link
ln -s $DATADIR/2157/2157_annotations.tsv.gz 2157_annotations.tsv.gz

cd ..

### Set variables

# members files
bacteria=./data/2_members.tsv.gz
archaea=./data/2157_members.tsv.gz
eukaryota=./data/2759_members.tsv.gz
# annotations files
archaea_fct=./data/2157_annotations.tsv.gz
# or download files from http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/

#directories
results_dir=./results/
bin_dir=./bin/

### 1. Which genes (OGs) occur in at least 99% of all genomes in the eggNOG5 database in each domain of life, respectively?

zcat $bacteria | $bin_dir/read_members_file.py -min_occurence 99 > $results_dir/cogs_bacteria_o99.txt
zcat $archaea | $bin_dir/read_members_file.py -min_occurence 99 > $results_dir/cogs_archaea_o99.txt
zcat $eukaryota | $bin_dir/read_members_file.py -min_occurence 99 > $results_dir/cogs_eukaryota_o99.txt


### 2. Which bacterial genes occur in at least 50% of all bacterial genomes, and in at least 99% thereof as single-copy?

zcat $bacteria | $bin_dir/read_members_file.py -min_occurence 50 -min_uniqueness 99 > $results_dir/cogs_bacteria_o50_u99.txt

# How many of these OGs were also identified as universal bacterial OGs (previous question)?
#comm -12 cogs_bacteria_o99.txt cogs_bacteria_o50_u99.txt

zcat $bacteria | $bin_dir/read_members_file.py -min_occurence 99 -min_uniqueness 99 > $results_dir/cogs_bacteria_o99_u99.txt

### 3. Identify all OGs that occur as single-copy in at least 97% of all archaea

zcat $archaea | $bin_dir/read_members_file.py -min_occurence_as_singlecopy 97 > $results_dir/cogs_archaea_os97.txt

# Are there archaea which lack 4 or more of those universal OGs?
zcat $archaea | $bin_dir/read_members_file.py -min_occurence_as_singlecopy 97 -missing 4


### 4. Compile an overview of the functional categories of these 121 archaeal OGs
zcat $archaea_fct | $bin_dir/annotate_cogs.py $results_dir/cogs_archaea_os97.txt archaea_counted_categories.txt ./data/functional_categories.txt
sort -nr archaea_counted_categories.txt > $results_dir/archaea_counted_categories_sorted.txt
rm archaea_counted_categories.txt