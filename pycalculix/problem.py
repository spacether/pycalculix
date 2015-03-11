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
            subprocess.call("%s %s" % (environment.CCX, self.fname),
                            shell=True)
            print('Solving done!')

            # select the probem's parts and load the results file
            if os.path.isfile(fname):
                self.solved = True
                self.rfile.load()
