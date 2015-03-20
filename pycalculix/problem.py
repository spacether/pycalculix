"""This module stores the Problem class which is used to solve different types
of analysis.
"""

import subprocess # used to launch ccx solver
import os # used to see if there is a results file

from . import environment
from . import base_classes
from . import results_file

class Problem(base_classes.Idobj):
    """Makes a problem which can be analyzed with Calculix ccx.

    Args:
        feamodel (FeaModel): the parent FeaModel
        problem_type (str): model type, options:

            - 'struct': structural
        fname (str): file prefix for the problem .inp and results files
            If value is '' it will default to the project name of the FeaModel
    Attributes:
        fea (FeaModel): parent FeaModel
        __ptype (str): problem type, options:

            - 'struct': structural
        solved (bool): boolean storing whether or not the analysis was solved
        fname (str): the problem database and results file prefix
        rfile (Results_File):
            Results_File is loaded in after the model has been solved
    """
    def __init__(self, feamodel, problem_type, fname=''):
        self.fea = feamodel
        self.__ptype = problem_type
        self.solved = False
        if fname == '':
            fname = self.fea.fname
        self.fname = fname
        self.rfile = results_file.ResultsFile(self)
        base_classes.Idobj.__init__(self)
        self.fea.problems.append(self)

    @staticmethod
    def __get_ntxt(nodes):
        """Returns list of strings defining all nodes.

        Args:
            nodes (list): list of all nodes
        """
        res = []
        res.append('*NODE, NSET=nodes')
        for node in nodes:
            res.append(node.ccx())
        print('INFO: %i nodes' % len(nodes))
        return res

    def __get_etxt(self, elements):
        """Returns list of strings defining all elements.

        Args:
            elements (list): list of all elements
        """
        res = []
        ccxtypes = set([e.ccxtype for e in elements])
        eall_written = False
        eall = []
        for ccxtype in ccxtypes:
            setname = ccxtype
            if len(ccxtypes) == 1:
                setname = 'EAll'
                eall_written = True
            res.append('*ELEMENT, TYPE=%s, ELSET=%s' % (ccxtype, setname))
            elist = [e for e in elements if e.ccxtype == ccxtype]
            print('INFO: %i %s elements' % (len(elist), ccxtype))
            for element in elist:
                res.append(element.ccx())
            eall += elist
            print('INFO: %i total elements' % len(eall))
        if eall_written == False:
            tmp = self.__get_eset('EALL', eall)
            res += tmp
        return res

    @staticmethod
    def __get_ctxt(components):
        """Returns list of strings defining all components.

        Args:
            components (list): list of all components
        """
        res = []
        for comp in components:
            res += comp.ccx()
        return res

    @staticmethod
    def __get_eset(name, elements):
        """Returns list of strings defining components of elements.

        Args:
            name (str): component name
            elements (list): list of component elements
        """
        res = []
        items_per_line = 6
        res.append('*ELSET,ELSET='+name)
        grouped_els = base_classes.chunk_list(elements, items_per_line)
        for group in grouped_els:
            item_ids = [str(e.id) for e in group]
            line = ', '.join(item_ids)
            if group != grouped_els[-1]:
                line += ','
            res.append(line)
        return res

    @staticmethod
    def __fix_line(line, fstr):
        """Fixes the line and returns a fixed line.

        Args:
            line (str): line to fix
            fstr (str): format string to use
        """
        items = fstr.split(',')
        node_num_size = int(items[2][1:])
        pre_len = 3 + node_num_size
        # we are missing + prefix on some numbers, some are too long
        start_str = line[:pre_len]
        end_str = line[pre_len:]
        fields = end_str.count('E')
        low_ind = 0
        values = []
        while len(values) < fields:
            high_ind = end_str.find('E', low_ind)+5
            value = float(end_str[low_ind:high_ind])
            values.append(value)
            low_ind = high_ind
        values = ['%12.5e' % val for val in values]
        end_str = ''.join(values)
        new_str = start_str + end_str + '\n'
        return new_str

    def __fix_frd(self):
        """Fixes the frd file on win32 systems. Text formatting fixed.

        On win32 Calculix, results file formatting is  not written correctly.
        Nodal results and stresses are not written as fixed length fields.
        """
        if 'win32' in environment.CCX:
            frd_file = self.fname+'.frd'
            lines = []
            try:
                with open(frd_file, "r") as infile:
                    lines = infile.readlines()
                numlines = len(lines)
                ind = 0
                fix = False
                fstr = ''
                while ind < numlines:
                    line = lines[ind]
                    if '1PSTEP' in line:
                        # we are in a results block
                        ind += 1
                        line = lines[ind]
                        format_ = int(line.split()[-1])
                        if format_ == 0:
                            fstr = "1X,I2,I5,6E12.5"
                        elif format_ == 1:
                            fstr = "1X,I2,I10,6E12.5"
                        ind += 1
                        line = lines[ind]
                        if 'DISPR' in line or 'FORCR' in line:
                            ind += 5
                        else:
                            ind += 7
                        fix = True
                        line = lines[ind]
                    if line[:3] == ' -3':
                        fix = False
                    if fix:
                        lines[ind] = self.__fix_line(line, fstr)
                    ind += 1
                with open(frd_file, "w") as outfile:
                    outfile.writelines(lines)
                print('File %s had its formatting fixed!' % frd_file)
            except IOError:
                print('Error reading .frd file!')

    def solve(self):
        """Solves the model in Calculix ccx."""
        inp = []

        # store what results we'll be outputting for each type of analysis
        out_el = {}
        out_el['struct'] = 'E,S' # strain, stress
        out_node = {}
        out_node['struct'] = 'RF,U' # reaction forces, displacement

        if self.__ptype == 'struct':
            # store nodes, elements, components
            box = {}
            box['nodes'] = []
            box['elements'] = []
            box['components'] = set()

            # store all loads in the parts in this model
            load_dict = self.fea.loads

            # store all nodes, elements, and part element sets
            box['nodes'] = self.fea.view.nodes
            box['elements'] = self.fea.view.elements

            # store all node and element components
            for time in load_dict:
                for load in load_dict[time]:
                    if load.ltype not in ['press', 'press_fluid']:
                        box['components'].add(load.comp)

            box['nodes'] = self.__get_ntxt(box['nodes'])
            box['elements'] = self.__get_etxt(box['elements'])
            box['components'] = self.__get_ctxt(box['components'])

            # add text definition for nodes, elelents, components
            inp += box['nodes']+['']
            inp += box['elements']+['']
            inp += box['components']+['']

            # read in all materials
            for matl in self.fea.matls:
                inp += matl.ccx()

            # write all steps and loads
            for time in load_dict:
                if time == 0:
                    # this is for thicknesses and materials
                    for load in load_dict[time]:
                        inp += load.ccx()
                else:
                    # only write times >= 1
                    inp.append('*STEP')
                    inp.append('*STATIC')

                    for load in load_dict[time]:
                        inp += load.ccx()

                    # make output frd file for cgx
                    inp.append('*EL FILE')
                    inp.append(out_el[self.__ptype])
                    inp.append('*NODE FILE')
                    inp.append(out_node[self.__ptype])

                    # make output dat file for integration point results
                    inp.append('*EL PRINT,ELSET=EALL')
                    inp.append('S')

                    # end step
                    inp.append('*END STEP')

            # write CCX inp file to the local directory
            fname = self.fname+'.inp'
            with open(fname, 'w') as outfile:
                for line in inp:
                    #print (line)
                    outfile.write(line+'\n')
            print('File: %s was written' % (fname))

            # run file
            runstr = "%s %s" % (environment.CCX, self.fname)
            print(runstr)
            subprocess.call(runstr, shell=True)
            print('Solving done!')

            # select the probem's parts and load the results file
            frd_file = self.fname+'.frd'
            if os.path.isfile(frd_file):
                self.solved = True
                self.__fix_frd()
                self.rfile.load()
            else:
                print('ERROR: results .frd file was not written!')
