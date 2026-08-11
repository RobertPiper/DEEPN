import sys, os

class printio():
    """Functions to print comments"""
    def __init__(self):
        pass

    def print_centered(self, string):
        print()
        if sys.stdout.isatty():
            if not os.popen('stty size', 'r').read():
                r, c = (100, 100)
            else:
                r, c = os.popen('stty size', 'r').read().split()
        else:
            r, c = (100, 100)
        for n in range(((int(c) + 1) - len(string)) // 2):
            sys.stdout.write(' ')
        print(string, "\n")

    def get_text_block(self, tag):
        found = False
        block = ''
        handle = open(os.path.join(os.curdir, "Y2Hreadme.txt"), 'r')
        for line in handle.readlines():
            if found:
                if line.strip() == ">>>END":
                    break
                block += line
            else:
                if line.strip() == ">>>" + tag:
                    found = True
        handle.close()
        return block

    def get_text_block_as_array(self, tag):
        found = False
        block = []
        handle = open(os.path.join(os.curdir, "Y2Hreadme.txt"), 'r')
        for line in handle.readlines():
            if found:
                if line.strip() == ">>>END":
                    break
                block.append(line.strip())
            else:
                if line.strip() == ">>>" + tag:
                    found = True
        handle.close()
        return block

    def print_progress(self, infile, x, y, z, a):
        if a == 1:
            sys.stdout.write('>>> Processing {}:\n'.format(infile))
            sys.stdout.flush()
        if a == 2:
            sys.stdout.write('>>> Counting. {} sequences with blast hits in file: {}'.format(x, infile))
            sys.stdout.flush()
        if a == 3:
            sys.stdout.write('>>> Finished with file {}.\n{} sequences with blast hits.\n{} total hits in {} '
                             'minutes'.format(infile, x, y, z))
            sys.stdout.flush()
        if a == 4:
            sys.stdout.write('\t\n   1st sequence:  {}\n   2nd sequence:  {}\n   3rd sequence:  {}\n\n'.format(x, y, z))
            sys.stdout.flush()
