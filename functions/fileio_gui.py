import os
import sys
import glob

class fileio():
    """Functions to process files and directories"""

    def __init__(self):
        pass

    def create_new_folder(self, directory, folder):
        folder_path = os.path.join(directory, folder)
        if not os.path.exists(os.path.join(directory, folder)):
            os.mkdir(folder_path)

    def change_file_name(self, directory, folder, oldsuffix, newsuffix):
        var = len(oldsuffix)
        for filename in glob.iglob(os.path.join(directory, folder, '*' + oldsuffix)):
            os.rename(filename, filename[:-var] + newsuffix)

    def get_file_list(self, directory, folder, suffix):
        fileList = os.listdir(os.path.join(directory, folder))
        returnList = []
        for file in fileList:
            if os.path.splitext(file)[1] == suffix:
                returnList.append(file)
        return(returnList)

    def get_exec_path(self, executable):
        if os.path.exists('/usr/local/bin/' + executable):
            return True
        if os.path.exists('/usr/bin/' + executable):
            return True
        return False

    def get_sam_filelist(self, directory, infolder):
        file_list = []
        if os.path.exists(os.path.join(directory, infolder)):
            file_list = self.get_file_list(directory, infolder, '.sam')
        return file_list

    def get_chromosomes_list(self, directory, tag, printio_handle):
        chromList = []
        for chr_name in printio_handle.get_text_block_as_array(tag):
            chromList.append(chr_name.rstrip())
        return chromList

    def remove_file(self, directory, folder, file_list):
        for fi in file_list:
            os.remove(os.path.join(directory, folder, fi))
            sys.stdout.write(">>> Cleaned up file %s in folder %s\n" % (fi, folder))

    def make_FASTA(self, junctionFile, outputFile):
        sys.stdout.write('>>> Writing %s as FASTA...\n' % os.path.split(junctionFile)[-1])
        sys.stdout.flush()
        counter = 0
        inFile = open(junctionFile, 'r')
        outFile = open(outputFile, 'w')
        for line in inFile:
            line.strip()
            split = line.split()
            try:
                outFile.write(">%s\n%s\n" % (str(split[0]), str(split[5])))
                counter += 1
            except Exception:
                continue
        inFile.close()
        outFile.close()
        sys.stdout.write('>>> Converted  %d junctions in %s to a FASTA file.\n' % (counter, os.path.split(
                junctionFile)[-1]))
        sys.stdout.flush()
