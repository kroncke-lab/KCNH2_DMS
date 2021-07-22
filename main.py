import models.labelvariant
import models.file_structure



# Define the folder where the fastq files are and generate a list of all fastq files
folder = "D:\\Tile4-Library4-LV369\\barcode-key\\20210716-6605-RU\\"  # input file

suffix, file_stem = models.file_structure.file_grab(folder)

#  second input expects a list, if given a string of characters it will not work, i.e., use '[file_stem][0]' not
#  'file_stem' with one input
var_file_list = models.file_structure.file_clean(folder, list(file_stem), suffix)

models.labelvariant.convert_to_mutation(list(var_file_list), folder)






