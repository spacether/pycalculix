import subprocess # used to launch ccx solver

from . import environment
from . import base_classes
from . import results_file

class Model(base_classes.Idobj):
    """Makes a model which can be analyzed with Calculix ccx.
    
    Args:
      parent (FeaModel): the parent FeaModel
      parts (Part or list of Part): stores the parts to run the analysis on
      mtype (str): model type, options:
        'struct': structural
    
    Attributes:
      p (FeaModel): parent FeaModel
      parts (list of Part): stores the parts to run the analysis on
      mtype (str): model type, options:
        'struct': structural
      rfile (None or Results_File): None by default
        Results_File is loaded in after the model has been solved
      
    """
    def __init__(self, parent, parts, mtype):
        self.p = parent
        if not isinstance(parts, list):
            parts = [parts]
        self.parts = parts
        self.mtype = mtype
        self.rfile = None
        base_classes.Idobj.__init__(self)

    def get_ntxt(self, nodes):
        """Returns list of strings defining all nodes.

        Args:
          nodes (list): list of all nodes
        """
        res = []
        res.append('*NODE, NSET=nodes')
        for n in nodes:
            res.append(n.ccx())
        return res

    def get_etxt(self, elements):
        """Returns list of strings defining all elements.

        Args:
          elements (list): list of all elements
        """
        res = []
        types = set([e.ccxtype for e in elements])
        eall_written = False
        es = []
        for t in types:
            tname = t
            if len(types) == 1:
                tname = 'Eall'
                eall_written = True
            res.append('*ELEMENT, TYPE='+t+', ELSET='+tname)            
            eset = [e for e in elements if e.ccxtype == t]
            es += eset
            for e in eset:
                res.append(e.ccx())
        if eall_written == False:
            tmp = self.get_eset('EALL', es)
            res += tmp
        return res

    def get_ctxt(self, components):
        """Returns list of strings defining all components.

        Args:
          components (list): list of all components
        """
        res = []
        for c in components:
            res += c.ccx()
        return res

    def get_eset(self, name, elements):
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
        
        if self.mtype == 'struct':
            # store nodes and elements
            N = []
            E = []
            C = []
            P = []

            # store all loads in the parts in this model
            these_loads = self.p.loads
            
            # store all nodes, elements, and part element sets
            for part in self.parts:
                E += part.elements
                N += part.nodes
                #P += self.get_eset(part.get_name(), part.elements)
            
            # store all nodal components
            for time in these_loads:
                for load in these_loads[time]:
                    if load.ltype not in ['press', 'press_fluid']:
                        C.append(load.comp)
            
            N = self.get_ntxt(N)
            E = self.get_etxt(E)
            C = self.get_ctxt(C)
            
            # nodes
            inp += N
            # elements
            inp += E
            # part components
            inp += P
            # load components
            inp += C
            
            # read in all materials
            for matl in self.p.matls:
                inp += matl.ccx()
            
            # write all steps and loads
            for time in these_loads:
                if time == 0:
                    # this is for thicknesses and materials
                    for load in these_loads[time]:
                        inp += load.ccx()
                else:
                    # only write times >= 1
                    inp.append('*STEP')
                    inp.append('*STATIC')
                    
                    for load in these_loads[time]:
                        inp += load.ccx()                

                    # make output frd file for cgx
                    inp.append('*EL FILE')
                    inp.append(out_el[self.mtype])
                    inp.append('*NODE FILE')
                    inp.append(out_node[self.mtype])
        
                    # make output dat file for integration point results 
                    inp.append('*EL PRINT,ELSET=EALL')
                    inp.append('S')

                    # end step
                    inp.append('*END STEP')
                        
            # write CCX i np file to the local directory
            fname = self.p.fname+'.inp'
            f = open(fname,'w')
            for line in inp:
                #print (line)
                f.write(line+'\n')        
            print ('File: '+ fname + ' was written')
            f.close()
            
            # run file
            subprocess.call("%s %s" % (environment.CCX, self.p.fname),
                                shell=True)
            print('Solving done!')
            
            # read the results file in
            self.rfile = results_file.Results_File(self, self.p.fname)
            self.p.select(list(self.parts))
