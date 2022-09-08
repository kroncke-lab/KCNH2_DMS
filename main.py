import models.file_structure
import models.labelvariant

# Define the folder where the fastq files are and generate a list of all fastq files
folder = "\\Users\\KRONCKE\\OneDrive\\Desktop\\KCNH2-sequencing-tests\\"  # input file

suffix, file_stem = models.file_structure.file_grab(folder)

#  second input expects a list, if given a string of characters it will not work, i.e., use '[file_stem][0]' not
#  'file_stem' with one input
var_file_list = models.file_structure.file_clean(folder, list(file_stem), suffix)

out_list = models.labelvariant.convert_to_mutation(list(var_file_list), folder)


'''
################### FIND BARCODE COUNTS FROM SORTED CELLS ###################
#  Define the folder where the fastq files are and generate a list of all fastq files
folder = "\\Users\\KRONCKE\\Desktop\\KCNH2-sequencing-tests\\"  # input file

suffix, file_stem = models.file_structure.file_grab(folder)

bc_file_list = models.file_structure.sorted_file_clean(file_stem, suffix)

#  count and combine
models.labelvariant.count_sorted(bc_file_list)
'''



