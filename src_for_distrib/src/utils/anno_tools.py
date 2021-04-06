#!usr/bin/python

# Tools for managing gff files
# These are really just files that contain a set of 9 tab delimited fields:
#  Chrom_name Data_origin Site_type start end . +/-(direction) . comments

# additional functions for splitting the last field of a gff entry
def newSplit(value):
        import shlex
        lex = shlex.shlex(value)
        lex.quotes = '"'
        lex.whitespace_split = True
        lex.commenters = ''
        return list(lex)

def make_comment_dict(gff_entry):
        gff_entry.comment_dict = {}
        keyvalue = gff_entry.comments.split(";")
        for pair in keyvalue:
                key_value= newSplit(pair)
                key = key_value[0]
                if len(key_value) < 2:
                        value = ""
                else:
                        value = "".join(key_value[1:])
                gff_entry.comment_dict[key] = value


class BEDGraphEntry:
    '''A simple container class for the equivalent of a bedgraph file line.
    '''
    # contig_name start end sitename score
    FORMAT_STRING = "{}\t{}\t{}\t{}"

    def __init__(self, line=None):

        if line is not None:
            self.parse_bedgraph_line(line)
        else:
            self.chrom_name = ""
            self.start = 0
            self.end = 0
            self.score = 0

    def parse_bedgraph_line(self, line):
        """
        Set this entry's values to those of a line from a gff file
        """

        datarray = line.rstrip().split("\t")
        self.chrom_name = datarray[0]
        self.start = datarray[1]
        self.end = datarray[2]
        self.score = datarray[4]

    def __repr__(self):
        """
        Return a formatted bedgraph line, which can be used to reconstitute the object or be written directly to a bedgraph file
        """
        return BEDGraphEntry.FORMAT_STRING.format(
            self.chrom_name,
            self.start,
            self.end,
            self.score
        )


class WigEntry:
    '''A simple container class for the equivalent of a wiggle file line.
    '''
    # contig_name start end sitename score
    FORMAT_STRING = "{}\t{}"

    def __init__(self, line=None):

        if line is not None:
            self.parse_wig_line(line)
        else:
            self.chrom_name = ""
            self.start = 0
            self.end = 0
            self.length = self.end - self.start
            self.score = 0

    #def parse_wig_line(self, line):
    #    """
    #    Set this entry's values to those of a line from a gff file
    #    """
    #
    #    datarray = line.rstrip().split("\t")
    #    self.chrom_name = datarray[0]
    #    self.start = datarray[1]
    #    self.end = datarray[2]
    #    self.score = datarray[4]

    def __repr__(self):
        """
        Return a formatted wig line, which can be used to reconstitute the object or be written directly to a bedgraph file
        """
        return WigEntry.FORMAT_STRING.format(
            self.start,
            self.score
        )


class AnnotationData:
    '''Class with some utilities for manipulating genome annotation files
    '''

    def __init__(self):
        self.data = []
        self.index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == len(self.data):
            self.index = 0
            raise StopIteration
        self.index = self.index + 1
        return self.data[self.index-1]

    def clear_db(self):
        self.data = []

    def cleanup(self):
        """
        sort entries based on starting position
        """

        self.data = list(set(self.data))
        # test this
        self.data = sorted(self.data, key=lambda x: (x.chrom_name, x.start))
  
    def find_entry(self, findfunc, findall=False):
        """
        Return the line or lines for which findfunc is true when given a gff item

        If findall is false, only the first such entry is returned

        Findfunc should take a gffline object as its only argument
        """
        matches = filter(findfunc, self.data)

        if len(matches) == 0:
            return []

        if (findall):
            return matches
        else:
            return matches[0]

    def add_entry(self, new_entry):
        # Add an externally constructed Entry object
        self.data.append(new_entry)

    def write_file(self, filename):
        """
        Write the current contents of my data to a file
        """
        with open(filename, "w") as ostr:
            for line in self:
                ostr.write("{}\n".format(line))


class WigData(AnnotationData):
    """
    Class for storing and manipulating bedgraph data
    """

    # super().__init__() keeps all parent attributes and methods
    def __init__(self):
        super().__init__()

    def addline(self, chrom_name, start, end, score):
        """
        Add a line with the given data
        """

        newobj = WigEntry()
        newobj.chrom_name = chrom_name
        newobj.start = start
        newobj.end = end
        newobj.score = score
        newobj.length = end - start

        self.data.append(newobj)

    def write_file(self, filename):
        self.cleanup()
        lengths = [rec.length for rec in self.data]
        length_count = len(set(lengths))
        if length_count == 1:
            fixed_step = True
        else:
            fixed_step = False

        with open(filename, "w") as ostr:
            ostr.write("track type=wiggle_0\n")
            prior_chrom = ''
            for line in self:
                if line.chrom_name != prior_chrom:
                    if fixed_step:
                        ostr.write("fixedStep chrom={} start=0 step={}\n".format(line.chrom_name,lengths[0]))
                    else:
                        ################## implement span in future
                        ostr.write("variableStep chrom={}\n".format(line.chrom_name))
                if fixed_step:
                    ostr.write("{}\n".format(line.score))
                else:
                    ostr.write("{}\t{}\n".format(line.start,line.score))
                prior_chrom = line.chrom_name



class BEDGraphData(AnnotationData):
    """
    Class for storing and manipulating bedgraph data
    """

    # super().__init__() keeps all parent attributes and methods
    def __init__(self):
        super().__init__()

    def parse_bedgraph_file(self, filename, clear=True):
        """
        Parse a bed file and store the lines in self.data

        The current contents of this object are overwritten iff clear 
        """

        if clear:
            self.clear_db()

        instr = open(filename, "r")
        for line in instr:
            if line.startswith("#"):
                continue

            newline = BEDGraphEntry(line)
            self.data.append(newline)

    def addline(self, chrom_name, start, end, score):
        """
        Add a line with the given data
        """

        newobj = BEDGraphEntry()
        newobj.chrom_name = chrom_name
        newobj.start = start
        newobj.end = end
        newobj.score = score

        self.data.append(newobj)

class GffEntry:
    """
    A simple container class for the equivalent of a gff file line
    """

    FORMAT_STRING = "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}"

    def __init__(self, line=None):

        if line is not None:
            self.parse_gff_line(line)
        else:
            self.chrom_name = ""
            self.data_origin = ""
            self.site_type = ""
            self.start = 0
            self.end = 0
            self.direction = "+"
            self.comments = ""


    def parse_gff_line(self, line):
        """
        Set this entry's values to those of a line from a gff file
        """

        datarray = line.rstrip().split("\t")
        self.chrom_name = datarray[0]
        self.data_origin = datarray[1]
        self.site_type = datarray[2]
        self.start = int(datarray[3])
        self.end = int(datarray[4])
        self.direction = datarray[6]
        self.comments = " ".join(datarray[8:])
        self.comments = self.comments.replace("\t", " ")
        self.comments = self.comments.replace("\n", "")

    def __repr__(self):
        """
        Return a formatted gff line, which can be used to reconstitute the object or be written directly to a gff file
        """
        return GffEntry.FORMAT_STRING.format(self.chrom_name, self.data_origin, self.site_type, self.start, self.end, self.direction, self.comments)


class GffData:
    """
    Class for storing and manipulating gff data
    """
    # super().__init__() keeps all parent attributes and methods
    def __init__(self):
        super().__init__()

    def parse_gff_file(self,filename, clear=True):
        """
        Parse a gff file and store the lines in self.data

        The current contents of this object are overwritten iff clear 
        """

        if (clear):
            self.clear_db()

        instr = open(filename, "r")
        for line in instr:
            if line[0] == "#":
                    continue

            newline = GffEntry(line)
            self.data.append(newline)

    def addline(self, chrom_name, data_origin, site_type,
                start, end, direction, comments):
        """
        Add a line with the given data
        """

        newobj = GffEntry()
        newobj.chrom_name = chrom_name
        newobj.data_origin = data_origin
        newobj.site_type = site_type
        newobj.start = start
        newobj.end = end
        #newobj.score = score
        newobj.direction = direction
        newobj.comments = comments

        self.data.append(newobj)


