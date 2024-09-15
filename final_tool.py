#setting up the logger first and foremost
#importing 'logging', prompting the user to install it if not available.
try:
    import logging
except ModuleNotFoundError:
    print(f"Hello, the logging module was not found, please make sure it is installed, and run the script again.")
    raise SystemExit(1)

#setting up two loggers - one to print serious errors in the terminal, the other to generate the log file.
file_logger = logging.getLogger('file_logger') #getting the logger. the name is passed because multiple loggers are being used.
file_logger.setLevel(logging.INFO) #setting the level to info as we want to write information into the log file.
fhandler = logging.FileHandler('log.log') #setting up filehandler to produce an output log file called 'log.log'.
fhandler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - {%(pathname)s:%(lineno)d} - %(message)s')) #formatting the logger messages
file_logger.addHandler(fhandler) #adding the handler to the logger

#repeating the same for a terminal logger
stream_logger = logging.getLogger('stream_logger')
stream_logger.setLevel(logging.ERROR)
shandler = logging.StreamHandler()
shandler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - {%(pathname)s:%(lineno)d} - %(message)s'))
stream_logger.addHandler(shandler)

#importing all the necessary modules, showing an error and exiting the program if they're not installed.
try:
    import matplotlib.pyplot as plt
    import os
    import gffutils
    import seaborn as sns
    from operator import attrgetter
    from statistics import mean 
    import pandas as pd
    import argparse
    import numpy as np
    import copy
except ModuleNotFoundError as e:
    stream_logger.error(f"Hello, the following module was not found, please make sure it is installed and run the script again.\n{e}")
    raise SystemExit(1)

#first create a large data structure - should be easy
#second, calculate positions of introns in every gene.
#third, match. 

#creating a command line
#parser is created, containing the program name, description, epilog, and formatter.
parser = argparse.ArgumentParser(prog='summer project',description='this program creates takes junctions data files, a genome sequence file, and a genome annotations file, and assesses accuracy of gene models.',epilog='thanks and regards',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#arguments are added
#here we create flags to take files/folder at the command line
parser.add_argument('--junctions', required=True, help='use this argument to pass the folder containing junctions data file(s). these may be organized in any way, it will recursively extract all files.')
parser.add_argument('--gff_whole', required=True, help='use this argument to pass the GFF file containing all the feature annotations for the whole genome (including intron features).')
parser.add_argument('--gff_introns_only', required=True, help=f"use this argument to pass the GFF file containing only 'intron' feature annotations.")
parser.add_argument('--genome_sequence', required=True, help='use this argument to pass the whole genome sequence fasta file.')
#adding the optional arguments with defaults
parser.add_argument('--AED_threshold', required=False, type= float, default= 0.5,help='use this argument to pass the desired AED (Annotation Edit Distance) cutoff. this is optional.')
parser.add_argument('--intron_cutoff', required=False, type= float, default= 0.2,help='use this argument to pass a cutoff for filtering introns. this is optional.')
#a list containing the arguments passed in the command line is created.
args = parser.parse_args()

#defining variables
folder = str(args.junctions)
gff_whole = str(args.gff_whole)
gff_introns_only = str(args.gff_introns_only)
genome_fasta = str(args.genome_sequence)
AED_threshold = float(args.AED_threshold)
cutoff = float(args.intron_cutoff)

#error handling to check that the paths supplied exist
if os.path.exists(folder) and os.path.exists(gff_whole) and os.path.exists(gff_introns_only) and os.path.exists(genome_fasta):
    pass
else:
    stream_logger.error(f"Please check that all the files/folders supplied exist. Not all were found. Exiting program.")
    raise SystemExit(1)

#now the program is starting.

##### class for introns ######
#this class was partly written by Dr. Kathryn Crouch based on a list structure in an earlier draft of the script written by me 
# Note: we don't actually care about non-unique counts - chosen not to capture this information
class Intron:
    _gene_list = False
    _in_annotation = False

    # constructor creates
    def __init__(self, contig, start, end, strand, unique, NU):
        self._key = f'{contig}:{start}-{end}:{strand}'
        self._contig = contig

        # error handling to check start and end are ints
        try:
            assert type(start) == int
            assert type(end) == int
        except AssertionError:
            stream_logger.error("During initialization of this intron object, start and end positions were not found to be integers. Please check and try again. Raising SystemExit.")
            raise SystemExit(1)
        
        self._start = (start)
        self._end = (end)

        # error handling to check + or -
        try:
            assert str(strand) == "+" or str(strand) == "-"
        except AssertionError:
            stream_logger.error("During initialization of this intron object, strand was not '+' or '-'. Please check and try again. Raising SystemExit.")
            raise SystemExit(1)
        
        self._strand = strand

        # error handling to check this is an int
        try:
            assert type(unique) == int
        except AssertionError:
            stream_logger.error("During initialization of this intron object, unique read depth was not found to be an integer. Please check and try again. Raising SystemExit.")
            raise SystemExit(1)
        
        self._unique = unique

    # getters

    @property
    def key(self):
        return self._key

    @property
    def contig(self):
        return self._contig

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    @property
    def unique(self):
        return self._unique

    @property
    def in_annotation(self):
        return self._in_annotation

    @property
    #if we have any records of this intron existing in any gene(s), this will be listed. else, it shows a result of the intron being in zero genes.
    def gene_list(self):
        if self._gene_list:
            return self._gene_list
        else:
            return 0

    # setters - only for attributes that we want to update

    # unique setter should have some checking for ints
    @unique.setter
    def unique(self, count):
        try:
            self._unique = int(count)
        except ValueError:
            stream_logger.error(f"The unique could not be set since it is not an integer. Please check, raising SystemExit.")
            raise SystemExit(1)

    @in_annotation.setter
    def in_annotation(self, bool):
        self._in_annotation = bool

    @gene_list.setter
    def gene_list(self, genes):
        self._gene_list = genes
 
    # method to increment unique counts
    def add_unique(self, count):
        # check that count is an int
        try:
            self.unique += int(count)
        except ValueError:
            stream_logger.error(f"The unique count could not be incremented since the supplied count is not an integer. Please check, raising SystemExit.")
            raise SystemExit(1)

    # make a list of genes that our intron is in
    def add_genes(self, *genes):
        if self.gene_list:
            for gene in genes:
                self.gene_list.append(gene)
        else:
            self.gene_list = list(genes)

    # what to show when we use print
    def __repr__(self):
        gene_ids = ', '.join([gene.id for gene in self.gene_list]) if self.gene_list else 'None'
        return f'{self.contig}:{self.start}-{self.end}:{self.strand}\nUnique counts:\t{self.unique}\nGenes:\t{gene_ids}\n'

#defining function to check clashes
def clash_checker(intronA, intronB):
    # KC: would not intronB.start <= intronA.end be a simpler test here? #MG: no, it would not cover clash case II (consult clash_cases diagram in appendix folder)

    if int(intronB.start) - int(intronA.end) <= 1: #i.e. B's start pos - A's end pos
        return True
    else:
        return False
    
#comparing all our predicted gene models (post filtering) to annotation models using AED 
def AED_calculator(ann_model, predicted_model):    
    #ann model will be the query, and predicted models will be the reference (created from experimental evidence)
    
    #both should be ordered by the start pos. 
    sorted_query = sorted(ann_model, key=lambda x: x.start)
    sorted_ref = sorted(predicted_model, key=lambda x: x.start)
   
    #mathematical concept taken from maker2 paper
    #citation for the same: Holt, C., & Yandell, M. (2011). MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. BMC bioinformatics, 12, 491. https://doi.org/10.1186/1471-2105-12-491

    #introns will be compared in pairs. for each pair of unknown position and length (with the only assumption being directionality and order)
    #we'll go intron-wise in the query 

    #there are 13 possible cases (check diagram). writing statements for these, and creating a variable to store intersection (base pair wise)

    intersection = 0
    
    for query_intron in sorted_query:
        #print("\nQuery Intron")
        #print(query_intron)
        for ref_intron in sorted_ref:
            #print("\nRef Intron")
            #print(ref_intron)
            query_end = int(query_intron.end)
            query_start = int(query_intron.start)
            ref_end = int(ref_intron.end)
            ref_start = int(ref_intron.start)

    #ann is query, junc is ref
            #case I - perfect match
            if query_end - ref_end == 0 and query_start - ref_start == 0:
                common_bp = (query_end - query_start) + 1
                intersection = intersection + common_bp

            #case II - ref_intron within query_intron
            elif query_end - ref_end > 0 and query_start - ref_start < 0:
                common_bp = (ref_end - ref_start) + 1
                intersection = intersection + common_bp

            #case III - query_intron within ref_intron
            elif query_end - ref_end < 0 and query_start - ref_start > 0:
                common_bp = (query_end - query_start) + 1
                intersection = intersection + common_bp

            #case IV - ref_intron ahead, query lagging
            elif query_end - ref_end > 0 and query_start - ref_start > 0 and query_start - ref_end < 0:
                common_bp = (ref_end - query_start) + 1
                intersection = intersection + common_bp

            #case V - query ahead, ref lagging
            elif query_end - ref_end < 0 and query_start - ref_start < 0 and query_end - ref_start > 0:
                common_bp = (query_end - ref_start) + 1
                intersection = intersection + common_bp

            #case VI - ref outside and after query
            elif query_end - ref_end < 0 and query_start - ref_start < 0 and query_end - ref_start < 0:
                #there are no more introns that have any common sequence with the query intron, breaking the loop
                break

            #case VII - same start, query is longer
            elif query_end - ref_end > 0 and query_start - ref_start == 0:
                common_bp = (ref_end - ref_start) + 1
                intersection = intersection + common_bp

            #case VIII - same start, reference is longer
            elif query_end - ref_end < 0 and query_start - ref_start == 0:
                common_bp = (query_end - query_start) + 1
                intersection = intersection + common_bp

            #case IX - same end, query is longer
            elif query_end - ref_end == 0 and query_start - ref_start < 0:
                common_bp = (ref_end - ref_start) + 1
                intersection = intersection + common_bp

            #case X - same end, ref is longer
            elif query_end - ref_end == 0 and query_start - ref_start > 0:
                common_bp = (query_end - query_start) + 1
                intersection = intersection + common_bp

            #case XI - exact trail, query first, then ref
            elif query_end - ref_end < 0 and query_start - ref_start < 0 and query_end - ref_start == 0:
                common_bp = 1
                intersection = intersection + common_bp

            #case XII - exact trail, ref first, then query
            elif query_end - ref_end > 0 and query_start - ref_start > 0 and query_start - ref_end == 0:
                common_bp = 1
                intersection = intersection + common_bp

            #case XIII - ref and query apart - ref first, then query
            elif query_end - ref_end > 0 and query_start - ref_start > 0 and query_start - ref_end > 0:
                pass #there is no commonality, hence no intersection to be counted

    #now calculating the number of bp in introns in the query
    query_bp = 0
    for query_intron in sorted_query:
        bp = (query_intron.end) - (query_intron.start) + 1
        query_bp = query_bp + bp

    #similarly, calculating for ref
    ref_bp = 0
    for ref_intron in sorted_ref:
        bp = (ref_intron.end) - (ref_intron.start) + 1
        ref_bp = ref_bp + bp

    #calculating sensitivity (SN)
    try:
        SN = intersection / ref_bp
    except ZeroDivisionError:
        SN = intersection / 0.000001

    #calculating specificity (SP)
    try:
        SP = intersection / query_bp
    except ZeroDivisionError: #there may be a case where query_bp is zero - the ann model containing zero introns. here, we can tweak it to be a very small value to avoid a Zero Division Error
        SP = intersection / 0.000001

    #calculating congruency (C)
    C = (SN + SP) / 2

    #calculating incongruency (D), also known as AED
    AED = 1 - C

    return AED

#creating a function to predict gene models from a list of Intron class objects
#this function contains some modifications done by Dr. Kathryn Crouch with the purpose of removing bugs, each of these is highlighted with a comment starting with #KC.
def gene_models_predictor(introns):
    
    #sorting directionally
    introns = sorted(introns, key=lambda x: x.start)  

    gene_models_list = []
    list1 = [] #first gene model

    for intron in introns:

        #we start with model number 1
        if len(list1) == 0: #this means we're adding the first intron in, this can be added as is, without any testing
            list1.append(intron)
            gene_models_list.append(list1.copy()) #KC: added .copy()

        else: #we have an existing set of model(s)
            
            #we need to check if there is a clash between current 'intron' and last intron in every gene model present in gene_models_list
            #how does this work? - if all are clashing, new model(s) needs to be prepared and added in, otherwise we can just add the intron into whichever gene model it is not clashing with

            #testing the former
            # KC changed this to any - you end up with more models to assess later, but you were missing some with all
            if any(clash_checker(model[-1], intron) for model in gene_models_list): 
                
                #what now? - in that case, we create a new list containing all models with the last intron removed. then duplicates are removed (since there might be models that are same with only the last intron differing), and this will be repeated until either a model is found which doesn't clash with the intron-of-interest, or until all models are emptied out, in which case a new model will be created that contains only the intron-of-interest.
                new_models = copy.deepcopy(gene_models_list) #KC: changed copy to deepcopy
                
                while len(new_models) != 0 and all(clash_checker(model[-1], intron) for model in new_models): #note - an empty model will throw an error with clash_checker, and the model[-1] bit. this won't be an issue since if the first condition in an 'and' statement evaluates false, the second isn't evaluated at all!

                    #we need to pop last item from each list. if after popping the model is empty, the entire model is removed.
                    models_to_remove = []
                    for model in new_models:
                        try:
                            model.pop()
                        except IndexError: #popping from empty list
                            pass
                        if len(model) == 0:
                            models_to_remove.append(model.copy()) #KC: added copy

                    for model in models_to_remove:
                        if model in new_models:
                            new_models.remove(model)

                    del models_to_remove  

                    #removing duplicates
                    new_models_without_duplicates = []
                    for model in new_models:
                        if model not in new_models_without_duplicates:
                            new_models_without_duplicates.append(model)
                    new_models = new_models_without_duplicates

                #once the while loop has stopped running, either the list is empty (there are no more models left), or there is at least one model which is not clashing with the intron of interest.
                #checking former
                if len(new_models) == 0:
                    new_gene_model = [intron]
                    gene_models_list.append(copy.deepcopy(new_gene_model)) #KC: added deepcopy

                #checking latter
                elif not all(clash_checker(model[-1], intron) for model in new_models):
                    #manually going through each potential new model. if there is still a clash, it can be skipped, if there is no clash, the new intron can be appended and the new model can be added to the master list
                    for model in new_models:
                        if clash_checker(model[-1], intron):
                            pass
                        else:
                            model.append(intron)
                            gene_models_list.append(copy.deepcopy(model)) #KC: added deepcopy

            #testing the latter
            else: #if this is true, intron may be safely added to at least one of the existing models

                for model in gene_models_list:
                    if clash_checker(model[-1], intron):
                        pass
                    else:
                        model.append(intron)

    return gene_models_list

#defining a function to translate DNA sequence
#function adapted from - https://www.geeksforgeeks.org/dna-protein-python-3/
def translate(DNA_seq): 
	
	codon_table = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				 
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
		'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
	} 
	protein_seq ="" 
	for i in range(0, len(DNA_seq), 3): 
		codon = DNA_seq[i:i + 3] 
		protein_seq += codon_table[codon] 
	return protein_seq 

#extracting all junctions files found in given folder (including all sub-folders)   
files = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(folder)) for f in fn]

dict_all = dict() #this is a dictionary that will contain all introns found in the junctions files.

#iterating through junction files
for file in files:
    with open(file) as in_file:

        #skipping header
        next(in_file)

        #iterating
        for line in in_file:

            # unpacking everything in one go             
            location, strand, unique, NU = line.rstrip().rsplit('\t')

            #creating key for dictionary
            key = str(location) + ':' + str(strand)

            #creating value for the key created
            contig = location.split(':')[0]
            start_pos = ((location.split(':')[1]).split('-')[0])
            end_pos = ((location.split(':')[1]).split('-')[1])
            # creating your value here isn't necessary - only create it below where you need it

            # we don't care about introns with no unique reads
            # by excluding these, we can cut down processing time for the GFF part
            if int(unique) >= 1:

                if key in dict_all:
                    # I made the intron data into a class - we can use the class method to increment these counts
                    # this is safer and more efficient than partially overwriting the list like you were doing
                    # note that although I wrote a lot of boiler plate above, the code here is now reduced from 10 lines to 1 line, and it's safer!
                    dict_all[key].add_unique(int(unique))

                else:
                    # This is now a dict where each value is an Intron object, with its identifier key
                    dict_all[key] = Intron(contig, int(start_pos), int(end_pos), strand, int(unique), int(NU))

file_logger.info(f"these are all the introns in the all read junctions files - {len(dict_all)}")

#writing a gff file listing all introns in the junction files.
with open('all_introns.gff','w') as out_file:
    for intron in dict_all.values():
        out_file.write(f"{intron.contig}\tVEuPathDB\tintron\t{intron.start}\t{intron.end}\t.\t{intron.strand}\t.\tID={intron.key}\n")
        #works as expected

file_logger.info(f"the file 'all_introns.gff' was successfully written. this contains gff records of all junctions supplied.")

dict_ann_only = dict() #this is a dictionary that contains introns found only in annotations, that were not found in the junctions files.

#now opening gff file - introns only
with open(gff_introns_only) as in_file:

    for line in in_file:

        line = line.rstrip().rsplit('\t')

        #reproducing the key
        key = str(line[0]) + ':' + str(line[3]) + '-' + str(line[4]) + ':' + str(line[6])

        if key in dict_all:
            dict_all[key].in_annotation= True

        else:
            dict_ann_only[key] = line
        
file_logger.info(f"these are all the introns found in annotations only, unsupported by reads - {len(dict_ann_only)}")

dict_read_only = dict() #this is a dictionary containing introns that we found only in reads, not annotated. we expect these to be high in number, the majority in fact, as there will be many instances of noise 'introns' with very low unique read depth

dict_both = dict() #this is a dictionary that contains introns found in both the reads (junctions files) and annotations

for key, intron in dict_all.items():
    
    if not intron.in_annotation:
        dict_read_only[key] = intron

    else:
        dict_both[key] = intron

file_logger.info(f"these are all the introns found in the reads only, i.e. not annotated - {len(dict_read_only)}")

file_logger.info(f"these are all the introns found in reads that were annotated - {len(dict_both)}") #for SRR6493537 (used to test a bit of this script), this is 30402 (total junctions in reads) - 3250 (junctions in reads not annotated) #30402 - 3250 = 27152

#creating histograms to show the distribution of unique read depth in reads_only_introns and reads_that_were_annotated.
read_only_unique = []
for intron in dict_read_only.values():
    read_only_unique.append(intron.unique)

both_unique = []
for intron in dict_both.values():
    both_unique.append(intron.unique)

#creating distribution plots

plt.hist(read_only_unique, color='lightgreen', ec='black', bins=15)
#the plot is titled
plt.title('Read Depth (Unique) in Introns Found in Reads Only')
#the plot is saved
plt.savefig('readsonly_depth.png')
#the plot object is erased.
plt.clf()

plt.hist(both_unique, color='lightgreen', ec='black', bins=15)
#the plot is titled
plt.title('Read Depth (Unique) in Introns Found in Reads That Were Annotated')
#the plot is saved
plt.savefig('readsandanns_depth.png')
#the plot object is erased.
plt.clf()

file_logger.info(f"two histograms were successfully created - named 'readsonly_depth.png' and 'readsandanns_depth.png'. \nboth saved in the working directory.")

#removing extra lists from memory
del both_unique, read_only_unique

#creating the genome annotation database from the gff file
db = gff_whole.replace('.gff','.db')
#if the database already exists, this is loaded. else, a new database is created.
if not os.path.isfile(db):
    db = gffutils.create_db(gff_whole,dbfn=db,force=True,keep_order=True,merge_strategy='merge',id_spec=['ID', 'Name'],sort_attribute_values=True)
else:
    db = gffutils.FeatureDB(db, keep_order=True)

#this section of the code (lines 567 - 622) was written by Dr. Kathryn Crouch as an improvement upon a previous code written by me that had similar functionality but less-than-ideal efficiency.

# it is more efficient to run one query to get all the genes, then iterate through it gene-wise
# also, since introns cannot cross chromosomes, I've split the intron dict into chromosomes so you only consider the introns that are on the correct chromosome each time
# this reduces the computational burden further

# make a new dict where the chromosomes are the keys
# the value is a list of Intron objects in that chromosome
dict_by_chr = {}
for key, intron in dict_all.items():
    if intron.contig in dict_by_chr:
        dict_by_chr[intron.contig].append(intron)
    else:
        dict_by_chr[intron.contig] = [intron]

# now we consider only one chromosome at a time
#while doing this, also creating a master list of ALL genes
genes_master = []

for chr, intronList in dict_by_chr.items():
    # sort the list of introns in this chromosome by start
    # this will make it easier for us to avoid iterating through every intron for every gene
    intronList.sort(key=lambda x: x.start)
    # I'm going to query just once per chromosome and iterate through the list of genes
    for gene in db.region(seqid=chr,featuretype='protein_coding_gene', completely_within=True):
        #adding gene to master list
        #creating the key
        key = str(gene.id) + ':' + str(gene.start) + '-' + str(gene.end) + ':' + str(gene.strand)
        #adding to list
        genes_master.append(key)
        
        # I believe this list of genes is ordered already
        
        # this list will help us to avoid upstream introns in future iterations
        introns_to_remove = 0

        # Note that we are only looking at introns that we know are already on the same contig as the current gene
        for intron in intronList:
            # one end of the intron is in a gene
            # we only want to add the gene to the intron if they're both on the same strand.
            if ((gene.end >= intron.start and intron.start >= gene.start and str(gene.strand) == str(intron.strand)) or (gene.end >= intron.end and intron.end >= gene.start and str(gene.strand) == str(intron.strand))):
                # note that we add the whole gene feature object here - this contains the location and strand information
                # you can use this later to determine if the intron is incompletely within a gene and if the strand matches
                dict_all[intron.key].add_genes(gene)
                
            # if the end of this intron is before the end of the current gene, the intron cannot span to the next gene            
            # we will remove it so that we do not waste time considering introns that cannot be in genes
            if intron.end <= gene.end:
                introns_to_remove += 1

            # if the intron starts AFTER this gene, there is no point in going through the rest of the list, so stop this loop
            if intron.start >= gene.end:
                break

        # we remove introns that we  know must be upstream of the next gene before the next loop so that we don't consider them
        # this must be done outside of the loop that iterates this list
        del intronList[0:introns_to_remove]

# a list of the genes is now stored in each Intron object, but if you wanted to access by gene
#creating a dict to store information on genes
gene_dict = {}
for intron in dict_all.values():
    if intron.gene_list:
        for gene in intron.gene_list:
            #creating the key
            key = str(gene.id) + ':' + str(gene.start) + '-' + str(gene.end) + ':' + str(gene.strand)
            #print(key) #example - TGME49_318880:5521-10311:+
            if key in gene_dict:
                gene_dict[key].append(intron)
            else:
                gene_dict[key] = [intron]

#setting a cutoff value
#calculating int/intronX over whole genome
int_by_intX = []

for gene in gene_dict.keys():
    #assigning 'Intron X' to a variable
    IntronX = max(gene_dict[gene], key=attrgetter('unique'))
    #print(IntronX)

    for intron in gene_dict[gene]:
            int_by_intX.append(int(intron.unique)/int(IntronX.unique))

    #print(int_by_intX)

#plotting this in a univariate plot
sns.violinplot(y=int_by_intX).figure.savefig('cutoff.png') #doing violin instead of swarmplot because the latter was taking hours and hours to run due to too many data points
#the plot object is erased.
plt.clf()
#based on this plot, setting default cutoff at 0.2 until 5000, after which we'll count all introns above 1000 as valid. #this will be added as a default to this program's command line.

#report the creation of plot in the logger
file_logger.info(f"a violin plot of distribution of cutoff values across the whole genome was created. this is named 'cutoff.png' and may be found in the working directory.")

#deleting extra stuff from memory
del int_by_intX, dict_ann_only, dict_both, dict_read_only

#creating counts for each flag
likely_valid_count = 0
likely_invalid_count = 0
unsure_count = 0

#now going gene by gene and writing both the gene output file, and appending predictions in gff format to the main supplied gff file.
with open(gff_whole,'a') as gff_file:
    with open('gene_info.csv','w') as out_file:

        #writing the header in the output file
        out_file.write(f"Gene\tFlag\tCase\tRef Transcript\tAED Scores\tTotal Introns Found\tNumber of Valid Introns\tMean Read Depth Before Filtering\tMean Read Depth After Filtering\tMaximum Read Depth\tCutoff Used\tEvidence of Validity of Ref Transcript?\tEvidence of Existence of Alternative Transcript?\tBest Alternative Model\n")

        #here, adding genes to gene dict which didn't coincide with *any* intron from junction files. because of the way the gene dict has been written, so far only those genes will be included for which at least one intron was found in junctions. we also want to report information on genes for which no introns were found
        for gene in genes_master:
            if gene in gene_dict:
                pass
            else:
                gene_dict[gene] = [] #the value is an empty list, which represents that no introns were found

        #iterating through genes. writing output accordingly.
        for gene, intron_list in gene_dict.items():

            #defining some variables
            gene_ID = gene.split(':')[0]
            gene_start = int(gene.split(':')[1].split('-')[0])
            gene_end = int(gene.split(':')[1].split('-')[1])
            gene_strand = gene.split(':')[2]
            total_introns = len(intron_list)

            #checking if any introns were found in junction files for this gene. if not, we can write this in with the 'unsure' flag
            if len(intron_list) == 0:
                case_ID = "C"
                flag = "Unsure - Intron Evidence not Found in Supplied Dataset, No Alternatives"
                out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t-\t-\t{total_introns}\t0\t-\t-\t-\t-\t-\tFalse\t-\n")
                unsure_count = unsure_count + 1
                continue

            #assigning 'Intron X' 
            IntronX = max(intron_list, key=attrgetter('unique'))
            
            #sorting the intron_list by start
            intron_list.sort(key=lambda x: x.start)

            #defining some more variables
            for intron in intron_list:
                contig = intron.contig
                break
            read_depth_values = [o.unique for o in intron_list]
            mean_read_depth_all = mean(read_depth_values) 
            max_read_depth = IntronX.unique

            #here we need to consider the special case of exon only genes - since IntronX will not be a true representation of the expression of that gene, we need to impose an absolute metric. for our purposes, an intron with more than 1000 reads is good enough.
            exon_only = False
            if len(list(db.children(id=gene_ID, completely_within= True, featuretype= 'intron'))) == 0:
                exon_only = True

            #filtering out invalid introns
            filtered_intron_list = []
            for intron in intron_list:
                if IntronX.unique > 5000 or exon_only == True:
                    cutoff_used = "More than 1000 (absolute)"
                    if intron.unique > 1000:
                        filtered_intron_list.append(intron)

                else:
                    cutoff_used = f"More than or equal to {cutoff} * max read depth"
                    if int(intron.unique)/int(IntronX.unique) >= cutoff:
                        filtered_intron_list.append(intron)
            
            del intron_list

            #defining some more variables
            total_introns_after_filtering = len(filtered_intron_list)
            read_depth_values = [o.unique for o in filtered_intron_list]
            if total_introns_after_filtering == 0:
                mean_read_depth_filtered = 0
            else:
                mean_read_depth_filtered = mean(read_depth_values)

            #for any exon-only gene for which no introns remained after filtering, we can say that they're valid
            if exon_only == True and len(filtered_intron_list) == 0:
                flag = "Likely Valid - Exon Only, No Alternatives"
                case_ID = "A"
                out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t-\t-\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\t-\t-\n")
                likely_valid_count = likely_valid_count + 1
                continue

            #if exon-only is False, meaning gene has an intron, but filtered_intron_list is still zero: that means that the record contains a gene with an intron that has not been found in junction files. 
            # KC: it likely means that the intron has been added based on other RNA-seq experiments that you are not considering here. It does not necessarily mean it is manually added
            if exon_only == False and len(filtered_intron_list) == 0:
                flag = "Unsure - Intron Evidence not Found in Supplied Dataset, No Alternatives"
                case_ID = "C"
                out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t-\t-\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\t-\t-\n")
                unsure_count = unsure_count + 1
                continue

            #iterating through introns to find and check flanking introns. the moment the first valid one is found, the loop is stopped and the information is recorded.
            #if there is a flanking intron, then either the start of the first intron will be out of bounds of the gene, or the end of one of the introns in the list will be out of bounds
            
            #checking if start is out of bounds
            if filtered_intron_list[0].start <= gene_start: 
                flag = "Likely Invalid - Flanking Intron, No Alternatives"
                case_ID = "D"
                out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t-\t-\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tFalse\t-\t-\n")
                likely_invalid_count = likely_invalid_count + 1
                continue

            #checking if end is out of bounds
            # KC: When you checked the start above, you just checked the first intron in the list # why not just check the last intron in the list here? Faster than iterating....
            #MG: The reason for iterating is that the list is arranged by start position (ascending), not end position. Further, this is a list of all (filtered) introns, and not a gene model, hence there may well be clashes. In this case, there may be introns in the middle of this list, which are longer than subsequent introns and go out of bounds, while shorter introns towards the end may be within bounds. For this reason, I iterate and break the loop at the first instance of a flanking intron. 
            #MG: The other option was to rearrange this list by end position in descending order, but I didn't want to mess with an order that is uniform and works.
            for intron in filtered_intron_list:
                continue_to_next_gene = False
                if intron.end >= gene_end:
                    flag = "Likely Invalid - Flanking Intron, No Alternatives"
                    case_ID = "E"
                    out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t-\t-\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tFalse\t-\t-\n")
                    continue_to_next_gene = True #once info has been written, we can move onto the next gene
                    likely_invalid_count = likely_invalid_count + 1
                    break

            if continue_to_next_gene:
                continue
            
            #if exon-only is true, for which at least one intron remained after filtering, we can say that it is probably not valid as it is not accounting for what the evidence suggests
            if exon_only == True and len(filtered_intron_list) != 0:
                flag = "Likely Invalid - Alternative(s) Available"
                case_ID = "B"
                #creating predictions
                predicted_gene_models = gene_models_predictor(filtered_intron_list)
                #since AED will be 1 (since the annotation is completely different from evidence, it contains no introns), we will not calculate it at all at all
                #instead, we'll just write these predictions in the gff file and calculate the best alternative model (the bit of code for this is repeated, this was a last minute addition so not much time for refactoring)
                
                #writing the gff file, and finding mean read depth
                ID_counter = 0
                mean_read_depth_list = []
                
                for prediction in predicted_gene_models:
                    ID_counter = ID_counter + 1
                    gff_file.write(f"{contig}\tVEuPathDB_Prediction\tPredicted_Transcript\t{gene_start}\t{gene_end}\t.\t{gene_strand}\t.\tID={gene_ID}_p{ID_counter};Parent={gene_ID}\n")
                    #writing introns first, it is straightforward
                    for counter, intron in enumerate(prediction):
                        #the ID for each intron will be {prediction.id}-I{counter+1} (because counter is zero indexed but these IDs in gff files are 1 indexed)
                        gff_file.write(f"{contig}\tVEuPathDB_Prediction\tintron\t{intron.start}\t{intron.end}\t.\t{gene_strand}\t.\tID={gene_ID}_p{ID_counter}-I{counter+1};Parent={gene_ID}_p{ID_counter};gene_id={gene_ID}\n")

                    #writing exons, now, we need a counter, we need to write the first and last exon manually, and the rest within a loop over pairs of introns
                    #creating this counter 
                    exon_counter = 1 #starting at one because this is for IDs and IDs are 1-indexed
                    #writing first exon
                    gff_file.write(f"{contig}\tVEuPathDB_Prediction\texon\t{gene_start}\t{(prediction[0].start)-1}\t.\t{gene_strand}\t.\tID={gene_ID}_p{ID_counter}-E{exon_counter};Parent={gene_ID}_p{ID_counter};gene_id={gene_ID}\n")
                    exon_counter = 2
                    #writing all middling exons
                    for intronA, intronB in zip(prediction, prediction[1:]):
                        gff_file.write(f"{contig}\tVEuPathDB_Prediction\texon\t{(intronA.end)+1}\t{(intronB.start)-1}\t.\t{gene_strand}\t.\tID={gene_ID}_p{ID_counter}-E{exon_counter};Parent={gene_ID}_p{ID_counter};gene_id={gene_ID}\n")
                        exon_counter = exon_counter + 1
                    #writing the last exon
                    gff_file.write(f"{contig}\tVEuPathDB_Prediction\texon\t{(prediction[-1].end)+1}\t{(gene_end)}\t.\t{gene_strand}\t.\tID={gene_ID}_p{ID_counter}-E{exon_counter};Parent={gene_ID}_p{ID_counter};gene_id={gene_ID}\n")
                    #gff written

                    read_depth = []
                    for intron in prediction:
                        read_depth.append(intron.unique)
                    mean_read_depth_list.append(mean(read_depth))

                #best model suggestion
                #finding which index has the highest number in mean_read_depth_list
                index_max = np.argmax(mean_read_depth_list)
                #this is the index of the model in predicted_gene_models_list that is the best alt model
                best_alt_model_ID = f"{gene_ID}_p{index_max+1}" #+1 added because index will be zero indexed but IDs are 1 indexed

                #writing in gene_info
                out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t-\t-\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\tTrue\t{best_alt_model_ID}\n")
                likely_invalid_count = likely_invalid_count + 1
                continue
            
            #here on out, we only have genes for which ALL the filtered introns are found within bounds
            # for all of these genes, we want to go transcript-wise.
            for transcript in (list(db.children(id=gene_ID, completely_within= True, featuretype= 'mRNA'))):
                
                ref_transcript_ID = transcript.id

                #for this transcript, finding the beginning position of the first CDS and the end position of the last CDS to find beginning and end of coding sequence
                #strand does not matter right now - we're just interested in the beginning and end. once the DNA sequence has been found
                CDS_start = int((list(db.children(id=transcript.id, completely_within= True, featuretype= 'CDS'))[0]).start)
                CDS_end = int((list(db.children(id=transcript.id, completely_within= True, featuretype= 'CDS'))[-1]).end)
                
                #creating predictions
                predicted_gene_models = gene_models_predictor(filtered_intron_list)
                
                #checking these predictions to see if there are any intermittent stop codons, if there are, these models are removed. 
                for prediction in predicted_gene_models:
                    
                    predictions_to_remove = []
                    #how would this work??? we'd have to declare the sequence invalid if the CDS start OR end is within the bounds of any intron. checking this first
                    for intron in prediction:
                        try: 
                            assert CDS_start < intron.start or CDS_start > intron.end
                            assert CDS_end < intron.start or CDS_end > intron.end
                        except AssertionError:
                            predictions_to_remove.append(prediction)
                            break

                for prediction in predictions_to_remove:
                    if prediction in predicted_gene_models:
                        predicted_gene_models.remove(prediction)
                
                #now for the actual sequence extraction, this is discussed in report's future work section
                
                #formatting this annotated transcript 
                ann_gene_model = []
                for ann_intron in list(db.children(id=transcript.id, completely_within= True, featuretype= 'intron')):
                    ann_intron = Intron(contig, int(ann_intron.start), int(ann_intron.end), ann_intron.strand, 0, 0) #unique and non-unique values set to zero because these are annotations, so this information does not apply
                    ann_gene_model.append(ann_intron)

                #checking if ann_gene_model is valid
                ann_seq = ""
                if gene_strand == "+":
                    CDSes = db.children(id=transcript.id,featuretype='CDS',order_by='start') #all the CDS children of the transcript are found, and are ordered by their start position.
                elif gene_strand == "-":
                    CDSes = db.children(id=transcript.id,featuretype='CDS',order_by='start',reverse=True) #all the CDS children of the transcript are found, and are ordered by their start position, and then reversed because the strand is negative.
                for CDS in CDSes:
                    ann_seq = ann_seq + str(CDS.sequence(genome_fasta,use_strand=True)) #use_strand = True ensures that the coding sequence is read with respect to the strand directionality.

                #our sequence should be divisible by 3. if not, it is not valid.
                if len(ann_seq) % 3 != 0:
                    ann_seq_valid = False
                else: #it can be translated and checked for stop codons 
                    try:
                        assert (translate(ann_seq)).count("*") == 1 and (translate(ann_seq)).endswith("*")
                        ann_seq_valid = True
                    except AssertionError:
                        ann_seq_valid = False

                #filtering out based on AED
                predictions_to_remove = []
                for prediction in predicted_gene_models:
                    if float(AED_calculator(ann_gene_model,prediction)) > AED_threshold:
                        predictions_to_remove.append(prediction)
                for prediction in predictions_to_remove:
                    if prediction in predicted_gene_models:
                        predicted_gene_models.remove(prediction)

                #if AED based filtering filtered out *all* the models, this could mean two things - 1) either there wasn't enough evidence to validate if annotation was correct or not, or 2) there was plenty of evidence but the annotation was *so* incorrect that it can't be used as even a rough yardstick to compare predictions to, and so the script is unable to recommend predictions. 
                #in the former case, the annotation may be correct or incorrect, in the latter, the annotation is wrong. but there is no way to know which it is (at least programmatically within the paradigm of this tool)
                #for this reason, the flag unsure, no alternatives is allotted. 
                if len(predicted_gene_models) == 0:
                    flag = "Unsure - Alternative(s) Not Available"
                    case_ID = "O"
                    unsure_count = unsure_count + 1
                    out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t-\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\t-\tFalse\t-\n")
                    continue

                #out of the models we have left, creating pandas dataframe to store information
                #class objects can't be properly stored in a pandas dataframe, so storing indexes to the predicted_gene_models list instead
                indices = []
                for i in range(0,len(predicted_gene_models),):
                    indices.append(i)
                
                #initiating dataframe
                df = pd.DataFrame(indices, columns=['Gene_Model'])
                
                #adding AED and mean read depth info
                AED_list = []
                mean_read_depth_list = []
                for prediction in predicted_gene_models:
                    AED_list.append(float(AED_calculator(ann_gene_model,prediction)))
                    read_depth = []
                    for intron in prediction:
                        read_depth.append(intron.unique)
                    mean_read_depth_list.append(mean(read_depth))
                df.insert(0,'AED',AED_list,True)
                df.insert(0,'Mean_Read_Depth',mean_read_depth_list,True)
                
                #sorting models based on AED
                df.sort_values('AED')
                
                #adding a column for IDs - the ID for each prediction will be {transcript.id}+'_p'+{index} (this may also be considered AED-based ranking)
                prediction_IDs = []
                for i in range(1,len(df.index)+1,):
                    prediction_IDs.append(f"{transcript.id}_p{i}") 
                df.insert(0,'IDs',prediction_IDs,True)

                #writing this df into the gff file 
                for index, row in df.iterrows():
                    #each row is one predicted transcript
                    gff_file.write(f"{contig}\tVEuPathDB_Prediction\tPredicted_Transcript\t{gene_start}\t{gene_end}\t.\t{gene_strand}\t.\tID={row['IDs']};Parent={gene_ID}\n")
                    model = predicted_gene_models[int(row['Gene_Model'])]
                    #writing introns first, it is straightforward
                    for counter, intron in enumerate(model):
                        #the ID for each intron will be {prediction.id}-I{counter+1} (because counter is zero indexed but these IDs in gff files are 1 indexed)
                        gff_file.write(f"{contig}\tVEuPathDB_Prediction\tintron\t{intron.start}\t{intron.end}\t.\t{gene_strand}\t.\tID={row['IDs']}-I{counter+1};Parent={row['IDs']};gene_id={gene_ID}\n")

                    #writing exons, now, we need a counter, we need to write the first and last exon manually, and the rest within a loop over pairs of introns
                    #creating this counter 
                    exon_counter = 1 #starting at one because this is for IDs and IDs are 1-indexed
                    #writing first exon
                    gff_file.write(f"{contig}\tVEuPathDB_Prediction\texon\t{gene_start}\t{(model[0].start)-1}\t.\t{gene_strand}\t.\tID={row['IDs']}-E{exon_counter};Parent={row['IDs']};gene_id={gene_ID}\n")
                    exon_counter = 2
                    #writing all middling exons
                    for intronA, intronB in zip(model, model[1:]):
                        gff_file.write(f"{contig}\tVEuPathDB_Prediction\texon\t{(intronA.end)+1}\t{(intronB.start)-1}\t.\t{gene_strand}\t.\tID={row['IDs']}-E{exon_counter};Parent={row['IDs']};gene_id={gene_ID}\n")
                        exon_counter = exon_counter + 1
                    #writing the last exon
                    gff_file.write(f"{contig}\tVEuPathDB_Prediction\texon\t{(model[-1].end)+1}\t{(gene_end)}\t.\t{gene_strand}\t.\tID={row['IDs']}-E{exon_counter};Parent={row['IDs']};gene_id={gene_ID}\n")
                #gff written

                #finally writing in the gene_info

                #creating a string to input as AED result
                AED_output = ""
                for index, row in df.loc[:, ['AED', 'IDs']].iterrows():
                    AED_result = str(row['IDs']) + " : " + str(row['AED']) + "; "
                    AED_output = AED_output + AED_result

                #finding the best alternative model using (highest) mean read depth
                try:
                    best_alt_model_ID = df.iloc[(df['Mean_Read_Depth'].idxmax())]['IDs']
                except ValueError: #if pandas is throwing a ValueError, that means that there were no predictions. in this case, there is no ID. this should be extremely rare.
                    best_alt_model_ID = "-"

                #let's consider cases to allot an appropriate flag, we have the following bits of information - 1) validity of ann_sequence, 2) AEDs, 3) Mean Read Depth of each predicted transcript, 4) Number of Predictions
                # KC: These make sense - but I would be a little cautious about saying "Valid" or "Invalid" - the language is a bit strong given the caveats. It's more like "this is the most likely" or "there is an alternative that is more likely" - you need to find a concise way to phrase that.
                #MG: adding a unique case ID to output and a table in the report for additional context, instead of making flags too long/complicated/diverse outside of three major categories - valid, invalid, unsure.
                
                #Case I - (Validity and AED == 0.0 found) + only one prediction - this means the annotation is highly likely to be valid, all evidence points it
                #Case II - (Validity and AED == 0.0 found) + more than one prediction - if the highest mean read depth is also p1, then same as above, else it is likely valid with alternatives available
                #Case III - Validity is true, but AED == 0.0 was not found, but there was only one prediction - this means that the prediction and ann are different, but ann is valid. in this case, we're not sure
                #Case IV - validity is false, but AED == 0.0 and there was only one prediction - this means that our prediction was the same as the (invalid) annotation so we don't know the real answer. hence, we're not sure
                #Case V - (validity and AED == 0.0) both are false, and there is only one prediction - in this case the program has found decisively that the ann is probably not valid and instead there is a prediction available
                #Case VI - validity is true but AED == 0.0 is false (i.e. the ann was not in predictions), and there are multiple predictions - unsure about ann's validity since it is not found in evidence-based predictions, but it is valid, hence this is unsure
                #Case VII - validity is false, AED == 0.0 was found, there are multiple predictions - likely invalid with alternatives available
                #Case VIII - (validity and AED == 0.0) both false, multiple predictions found, same as above - likely invalid with alternatives available

                #Case I
                if ann_seq_valid == True and df.iloc[0]['AED'] == 0.0 and len(df.index) == 1:
                    flag = "Likely Valid - No Alternatives"
                    case_ID = "F"
                    likely_valid_count = likely_valid_count + 1
                    out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t{AED_output}\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\tFalse\t-\n")

                #Case II
                elif ann_seq_valid == True and df.iloc[0]['AED'] == 0.0 and len(df.index) > 1:
                    if best_alt_model_ID == df.iloc[0]['IDs']:
                        flag = "Likely Valid - No Alternatives"
                        case_ID = "G"
                        likely_valid_count = likely_valid_count + 1
                        out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t{AED_output}\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\tFalse\t-\n")
                    else:
                        flag = "Likely Valid - Alternative(s) Available"
                        case_ID = "H"
                        likely_valid_count = likely_valid_count + 1
                        out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t{AED_output}\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\tTrue\t{best_alt_model_ID}\n")

                #Case III
                elif ann_seq_valid == True and df.iloc[0]['AED'] != 0.0 and len(df.index) == 1:
                    flag = "Unsure - Alternative(s) Available"
                    case_ID = "I"
                    unsure_count = unsure_count + 1
                    out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t{AED_output}\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\tTrue\t{best_alt_model_ID}\n")

                #Case IV
                elif ann_seq_valid == False and df.iloc[0]['AED'] == 0.0 and len(df.index) == 1: 
                    flag = "Unsure - Alternative(s) Not Available"
                    case_ID = "J"
                    unsure_count = unsure_count + 1
                    out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t{AED_output}\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\tFalse\t-\n")

                #Case V
                elif ann_seq_valid == False and df.iloc[0]['AED'] != 0.0 and len(df.index) == 1:
                    flag = "Likely Invalid - Alternative(s) Available"
                    case_ID = "K"
                    likely_invalid_count = likely_invalid_count + 1
                    out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t{AED_output}\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tFalse\tTrue\t{best_alt_model_ID}\n")

                #Case VI
                elif ann_seq_valid == True and df.iloc[0]['AED'] != 0.0 and len(df.index) > 1:
                    flag = "Unsure - Alternative(s) Available"
                    case_ID = "L"
                    unsure_count = unsure_count + 1
                    out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t{AED_output}\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\tTrue\t{best_alt_model_ID}\n")

                #Case VII
                elif ann_seq_valid == False and df.iloc[0]['AED'] == 0.0 and len(df.index) > 1:
                    flag = "Likely Invalid - Alternative(s) Available"
                    case_ID = "M"
                    likely_invalid_count = likely_invalid_count + 1
                    out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t{AED_output}\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tTrue\tTrue\t{best_alt_model_ID}\n")

                #Case VIII
                elif ann_seq_valid == False and df.iloc[0]['AED'] != 0.0 and len(df.index) > 1:
                    flag = "Likely Invalid - Alternative(s) Available"
                    case_ID = "N"
                    likely_invalid_count = likely_invalid_count + 1
                    out_file.write(f"{gene_ID}\t{flag}\t{case_ID}\t{ref_transcript_ID}\t{AED_output}\t{total_introns}\t{total_introns_after_filtering}\t{mean_read_depth_all}\t{mean_read_depth_filtered}\t{max_read_depth}\t{cutoff_used}\tFalse\tTrue\t{best_alt_model_ID}\n")

file_logger.info(f"the output file 'gene_info.csv' has been written, this is found in the working directory. further, predicted gene models have been added to the supplied whole genome gff file.\n")

#adding valid, invalid, and unsure counts to final log file
file_logger.info(f"here are the cumulative counts of the flags:\nlikely valid: {likely_valid_count}\nlikely invalid: {likely_invalid_count}\nunsure: {unsure_count}\n")