import re # used to get info from frd file
from numpy import roots, linspace # need to find S1-S3, and to make contours
import os #need to check if results file exists
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

from . import environment # needed for DPI
from . import base_classes # needed for RESFIELDS

class Results_File(object):
    """Makes a results file.

    Args:
    solved_model (Model): parent model
    fname (str): project/results file prefix
    
    Attributes:
    p (Model): parent model
    fname (str): results file prefix
    steps (list): a list of float time steps
    results (dict): a dict sotring the results file.
      results[step]['node'][nnum][field] --> value
      results[step]['element'][enum][field] --> value
      field = 'ux' or 'Seqv' or 'ey' etc.
    time (float): current time we are looking at, defaults to -1
      When a file is loaded in, the first time step is loaded.
    """

    def __init__(self, solved_model, fname):
        self.p = solved_model
        self.fname = fname
        self.steps = [] # this stores a list of time steps
        self.results = {} # stores results, nested dicts
        self.time = -1
        self.read_frd() # read nodal results
        self.read_dat() # read element integration pt results

    def nplot(self, field, fname='', display=True, levels=21, gradient=False, gmult=1.0):
        # plot the results in the given field for the selected object

        # store the selected nodes and elements
        sel = {}
        sel['nodes'] = self.p.p.sel['nodes']
        sel['elements'] = self.p.p.sel['elements']
        sel['faces'] = self.p.p.sel['faces']
        
        # sort nodes low to high so index is correct
        # we have index to id below so showing subsets works
        sel['nodes'] = list(sel['nodes'])
        sel['nodes'] = sorted(sel['nodes'], key=lambda k: k.id)
        
        '''
        faceplot = False
        if len(sel['elements']) == 0 and len(sel['faces']) > 0
        NEED TO ADD CODE HERE TO MAKE INTERPOLATION TRIANGLES FOR ELEMENT FACES
        '''
        
        # store results at nodes
        axials = []
        radials = []
        zs = []
        id_to_ind = {}
        for n in sel['nodes']:
            id_to_ind[n.id] = len(axials)
            ax = n.y + gmult*self.results[self.time]['node'][n.id]['uy']
            rad = n.x + gmult*self.results[self.time]['node'][n.id]['ux']
            axials.append(ax)
            radials.append(rad)
            zs.append(self.results[self.time]['node'][n.id][field])

        # make a list of triangles, given by indices, looping anticlockwise
        triangles = []
        mylist = []
        if len(sel['elements']) > 0:
            mylist = sel['elements']
        elif len(sel['faces']) > 0:
            mylist = sel['faces']
        for e in mylist:
            tris = e.get_tris()     # list of triangle nodes defined by node id
            for t in tris:
                for ind, nid in enumerate(t):
                    t[ind] = id_to_ind[nid]     # convert id to index
            triangles += tris
        
        # check to see if selected nodes and elements are
        # in the parent model's nodes and elements
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # need to set tick list here
        vmin = min(zs)
        vmax = max(zs)
        tick_list = [vmin]
        if vmax != vmin:
            # we have a range of values we're plotting
            tick_list = linspace(vmin,vmax,levels)        
        
        # plot using a gradient(shaded) or levels
        if gradient:
            # This one is shaded
            plt.tripcolor(axials, radials, triangles, zs, shading='gouraud')
        else:
            # this one is not shaded
            plt.tricontourf(axials, radials, triangles, zs, levels=tick_list)

        # code required for the colorbar
        cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
        if vmax != vmin:
            # we have a range of values we're plotting
            if gradient:
                cmap = plt.cm.jet
            else:
                cmap = plt.get_cmap('jet', levels)
        scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=cmap)            
        scalarMap._A = [] # need to set this for it to work
        plt.colorbar(scalarMap, orientation='vertical', ticks=tick_list)
            
        # set the horizontal and vertical axes
        vert = max(radials) - min(radials)
        horiz = max(axials) - min(axials)
        vadder = (vert)/5
        hadder = (horiz)/5    
        (vmax,vmin) = (max(radials)+vadder, min(radials)-vadder)
        (hmax,hmin) = (max(axials)+hadder, min(axials)-hadder)
        plt.xlim(hmin, hmax)
        plt.ylim(vmin, vmax)
        
        # set units
        [f_unit, d_unit, t_unit] = self.p.p.get_units(field, 'dist', 'time')

        # set plot axes
        plt.title('Node %s%s\nTime=%f%s' % (field, f_unit, self.time, t_unit))
        plt.xlabel('axial, y'+d_unit)
        plt.ylabel('radial, x'+d_unit)
        ax.set_aspect('equal')
        if gmult != 1:
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])

        if fname != '':
            # save the image
            fname += '.png'
            if environment.DPI != None:
                plt.savefig(fname, dpi=environment.DPI, bbox_inches='tight')
            else:
                plt.savefig(fname, bbox_inches='tight')
            print('File %s was saved' % (fname))
        
        if display:
            plt.tight_layout()
            plt.show()
            
        # remove all figures
        plt.close()

    def eplot(self, field, fname='', display=True, levels=21, gradient=False, gmult=1.0):
        # plot the results in the given field for the selected object

        # store the selected nodes and elements
        sel = {}
        sel['nodes'] = self.p.p.sel['nodes']
        sel['elements'] = self.p.p.sel['elements']
        sel['faces'] = self.p.p.sel['faces']
        
        # sort nodes low to high so index is correct
        # we have index to id below so showing subsets works
        sel['nodes'] = list(sel['nodes'])
        sel['nodes'] = sorted(sel['nodes'], key=lambda k: k.id)
        
        '''
        faceplot = False
        if len(sel['elements']) == 0 and len(sel['faces']) > 0
        NEED TO ADD CODE HERE TO MAKE INTERPOLATION TRIANGLES FOR ELEMENT FACES
        '''
        
        # store results at nodes
        axials = []
        radials = []
        zs = []
        id_to_ind = {}
        for n in sel['nodes']:
            id_to_ind[n.id] = len(axials)
            ax = n.y + gmult*self.results[self.time]['node'][n.id]['uy']
            rad = n.x + gmult*self.results[self.time]['node'][n.id]['ux']
            axials.append(ax)
            radials.append(rad)

        # make a list of triangles, given by indices, looping anticlockwise
        triangles = []
        mylist = []
        if len(sel['elements']) > 0:
            mylist = sel['elements']
        elif len(sel['faces']) > 0:
            mylist = sel['faces']
        for e in mylist:
            val = self.results[self.time]['element'][e.id][field]
            tris = e.get_tris()     # list of triangle nodes defined by node id
            for t in tris:
                zs.append(val)
                for ind, nid in enumerate(t):
                    t[ind] = id_to_ind[nid]     # convert id to index
            triangles += tris
        
        # check to see if selected nodes and elements are
        # in the parent model's nodes and elements
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        # need to set tick list here
        vmin = min(zs)
        vmax = max(zs)
        tick_list = [vmin]
        if vmax != vmin:
            # we have a range of values we're plotting
            tick_list = linspace(vmin,vmax,levels)        
        
        # plot using levels
        plt.tripcolor(axials, radials, triangles, zs, shading='flat')

        # code required for the colorbar
        cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
        if vmax != vmin:
            # we have a range of values we're plotting
            if gradient:
                cmap = plt.cm.jet
            else:
                cmap = plt.get_cmap('jet', levels)
        scalarMap = cmx.ScalarMappable(norm=cNorm,cmap=cmap)            
        scalarMap._A = [] # need to set this for it to work
        plt.colorbar(scalarMap, orientation='vertical', ticks=tick_list)
            
        # set the horizontal and vertical axes
        vert = max(radials) - min(radials)
        horiz = max(axials) - min(axials)
        vadder = (vert)/5
        hadder = (horiz)/5    
        (vmax,vmin) = (max(radials)+vadder, min(radials)-vadder)
        (hmax,hmin) = (max(axials)+hadder, min(axials)-hadder)
        plt.xlim(hmin, hmax)
        plt.ylim(vmin, vmax)
        
        # set units
        [f_unit, d_unit, t_unit] = self.p.p.get_units(field, 'dist', 'time')

        # set plot axes
        plt.title('Element %s%s\nTime=%f%s' %
                  (field, f_unit, self.time, t_unit))
        plt.xlabel('axial, y'+d_unit)
        plt.ylabel('radial, x'+d_unit)
        ax.set_aspect('equal')
        if gmult != 1:
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])

        if fname != '':
            # save the image
            fname += '.png'
            if environment.DPI != None:
                plt.savefig(fname, dpi=environment.DPI, bbox_inches='tight')
            else:
                plt.savefig(fname, bbox_inches='tight')
            print('File %s was saved' % (fname))
        
        if display:
            plt.tight_layout()
            plt.show()
            
        # remove all figures
        plt.close()
        
    def set_time(self, val):
        # sets the current time to time val
        self.time = val
        print('Results file time set to: %f' % (self.time))

    def utot(self, vals):
        # computes sum of the squares
        res = [a**2 for a in vals]
        res = ( sum(res) )**0.5
        return res
    
    def seqv(self, vals):
        # calculates the equivalent stress
        [s11,s22,s33,s12,s13,s23] = vals
        a = s11 - s22
        b = s22 - s33
        c = s33 - s11
        d = s12**2 + s23**2 +s13**2
        res = (0.5*(a**2 + b**2 + c**2 +6*d))**0.5
        return res

    def principals(self, vals):
        # calculates and returns principal stresses, S1, S2, S3
        [s11,s22,s33,s12,s13,s23] = vals
        a = 1
        b = (s11 + s22 + s33)*-1.0
        c = (s11*s22 + s11*s33 + s22*s33 - s12**2 - s13**2 - s23**2)
        d = (s11*s22*s33 + 2*s12*s13*s23 - s11*(s23**2) - s22*(s13**2) - s33*(s12**2))*-1.0
        res = list(roots([a,b,c,d]))
        res = sorted(res, reverse=True)        
        return res
    
    def get_nmax(self, field):
        # returns the max value of a given results field
        res = [ ndict[field] for ndict in self.results[self.time]['node'].values() ]
        return max(res)

    def get_nmin(self, field):
        # returns the min value of a given results field
        res = [ ndict[field] for ndict in self.results[self.time]['node'].values() ]
        return min(res)
    
    def get_fsum(self, item):
        # returns fsum on nodes under given item (point or line)
        (fx,fy,fz) = ( [],[],[] )
        nodes = item.nodes
        nodes = [n.id for n in nodes]
        for n in nodes:
            x = self.results[self.time]['node'][n]['fx']
            y = self.results[self.time]['node'][n]['fy']
            z = self.results[self.time]['node'][n]['fz']
            if x != 0 or y != 0 or z != 0:
                #print('Node %i, (fx, fy, fz) = (%12.11f,%12.11f,%12.11f)' % (n, x, y, z))
                fx.append(x)
                fy.append(y)
                fz.append(z)
        fx = sum(fx)
        fy = sum(fy)
        fz = sum(fz)
        return [fx, fy, fz]
    
    def get_emax(self, stype):
        # returns the max stress of given type
        res = 0.0
        for step in self.steps:
            for (enum, edict) in self.results[step]['element'].items():
                for (ipnum, ipdict) in edict.items():
                    if ipdict[stype] > res:
                        res = ipdict[stype]
        return res

    def get_vals(self, fstr, line):
        # this returns a list of items based on an input format string
        res = []
        fstr = fstr.split(',')
        for item in fstr:
            if item[0] == "'":
                # strip off the char quaotes
                item = item[1:-1]
                # this is a string entry, grab the val out of the line
                ind = len(item)
                fwd = line[:ind]
                line = line[ind:]
                res.append(fwd)
            else:
                # format is: 1X, A66, 5E12.5, I12
                # 1X is number of spaces                
                (m,c) = (1, None)
                m_pat = re.compile(r'^\d+') # find multiplier
                c_pat = re.compile(r'[XIEA]') # find character
                if m_pat.findall(item) != []:
                    m = int(m_pat.findall(item)[0])
                c = c_pat.findall(item)[0]
                if c == 'X':
                    # we are dealing with spaces, just reduce the line size
                    line = line[m:]
                elif c == 'A':
                    # character string only, add it to results
                    fwd = line[:m].strip()
                    line = line[m:]
                    res.append(fwd)
                else:
                    # IE, split line into m pieces
                    w_pat = re.compile(r'[IE](\d+)') # find the num after char
                    w = int(w_pat.findall(item)[0])
                    for i in range(m):
                        # only add items if we have enough line to look at
                        if w <= len(line):
                            substr = line[:w]
                            line = line[w:]
                            substr = substr.strip() # remove space padding
                            if c == 'I':
                                substr = int(substr)
                            elif c == 'E':
                                substr = float(substr)
                            res.append(substr)
        return res
                
    def read_frd(self):
        # this reads a ccx frd results file which contains nodal results
        fname = self.fname+'.frd'
        if os.path.isfile(fname):
            
            print('Reading results file: '+fname)
            f = open(fname)
            mode = 'off'
            time = 0.0
            FORMAT = None
            rfstr = ''
            node_ids = []            
            while True:
                line = f.readline()
                if not line:
                    break
                
                #-------------
                # set the modes
                #--------------
                if '2C' in line:
                    # we are in nodal definition
                    fstr = "1X,'   2','C',18X,I12,37X,I1"
                    t = self.get_vals(fstr, line)
                    [KEY,CODE,NUMNOD,FORMAT] = t

                    # set results format to short, long or binary
                    # only short and long are parsed so far
                    if FORMAT == 0:
                        # short format
                        rfstr = "1X,'-1',I5,3E12.5"
                    elif FORMAT == 1:
                        # long format
                        rfstr = "1X,'-1',I10,3E12.5"
                    elif format == 2:
                        # binary
                        pass
                    
                    line = f.readline()
                    mode = 'nodes'
                    print('Reading '+mode)

                elif '1PSTEP' in line:
                    # we are in a results block
                    
                    # next line has the time in it
                    line = f.readline()
                    fstr = "1X,' 100','C',6A1,E12.5,I12,20A1,I2,I5,10A1,I2"
                    t = self.get_vals(fstr, line)
                    [KEY,CODE,SETNAME,VALUE,NUMNOD,TEXT,ICTYPE,NUMSTP,ANALYS,FORMAT] = t

                    # set results format to short, long or binary
                    # only short and long are parsed so far
                    if FORMAT == 0:
                        # short format
                        rfstr = "1X,I2,I5,6E12.5"
                    elif FORMAT == 1:
                        # long format
                        rfstr = "1X,I2,I10,6E12.5"
                    elif format == 2:
                        # binary
                        pass

                    # set the time
                    time = VALUE
                    if time not in self.steps:
                        self.steps.append(time)
                    if time not in self.results:
                        self.results[time] = {'node':{},'element':{}}
                        for nid in node_ids:
                            # make dict for each node
                            self.results[time]['node'][nid] = {}
                        
                    # get the name to determine if stress or displ
                    line = f.readline()
                    fstr = "1X,I2,2X,8A1,2I5"
                    t = self.get_vals(fstr, line)
                    [KEY, NAME,NCOMPS,IRTYPE] = t
                    
                    if NAME == 'DISPR':
                        mode = 'displ'
                        # read 4 lines to get to data
                        for i in range(4):
                            f.readline()
                    elif NAME == 'STRESSR':
                        mode = 'stress'
                        # read 6 lines to get to data
                        for i in range(6):
                            f.readline()
                    elif NAME == 'TOSTRAIR':
                        mode = 'strain'
                        # read 6 lines to get to data
                        for i in range(6):
                            f.readline()
                    elif NAME == 'FORCR':
                        mode = 'force'
                        # read 4 lines to get to data
                        for i in range(4):
                            f.readline()
                            
                    print('Reading '+mode+' storing: '+
                          ','.join(base_classes.RESFIELDS[mode]))
                    line = f.readline()
                
                #----------------------------
                # Store results
                #----------------------------
                
                # reset the read mode if we hit an end of record
                if line[:3] == ' -3':
                    mode = 'off'
                
                if mode == 'nodes':
                    # node definition, store node numbers only
                    t = self.get_vals(rfstr, line)
                    [KEY, NODE, x, y, z] = t
                    node_ids.append(NODE)
                    
                elif mode == 'displ':  
                    # displacements
                    t = self.get_vals(rfstr, line)
                    [KEY, NODE, ux, uy, uz] = t                    
                    labs = base_classes.RESFIELDS[mode]
                    vals = [ux, uy, uz]
                    utot = self.utot(vals)
                    vals.append(utot)
                    #print(vals)
                    d = self.results[time]['node'][NODE]
                    for (k, val) in zip(labs, vals):
                        d[k] = val

                elif mode == 'stress':
                    # stresses
                    t = self.get_vals(rfstr, line)
                    [KEY, NODE, sx, sy, sz, sxy, syz, szx] = t
                    labs = base_classes.RESFIELDS[mode]
                    vals = [sx, sy, sz, sxy, syz, szx]
                    seqv = self.seqv(vals)
                    [s1, s2, s3] = self.principals(vals)
                    vals.append(seqv)
                    vals += [s1, s2, s3]
                    #print(vals)
                    d = self.results[time]['node'][NODE]
                    for (k, val) in zip(labs, vals):
                        d[k] = val

                elif mode == 'strain':
                    # strains
                    t = self.get_vals(rfstr, line)
                    [KEY, NODE, ex, ey, ez, exy, eyz, ezx] = t
                    labs = base_classes.RESFIELDS[mode]
                    vals = [ex, ey, ez, exy, eyz, ezx]
                    eeqv = self.seqv(vals)
                    [e1, e2, e3] = self.principals(vals)
                    vals.append(eeqv)
                    vals += [e1, e2, e3]
                    #print(vals)
                    d = self.results[time]['node'][NODE]
                    for (k, val) in zip(labs, vals):
                        d[k] = val
                
                elif mode == 'force':
                    # reaction forces
                    t = self.get_vals(rfstr, line)
                    [KEY, NODE, fx, fy, fz] = t
                    labs = base_classes.RESFIELDS[mode]
                    vals = [fx, fy, fz]
                    d = self.results[time]['node'][NODE]
                    for (k, val) in zip(labs, vals):
                        d[k] = val                    

            f.close()
            print('The following times have been read: ',self.steps)
            print('Done reading file: %s' % (fname))
            self.set_time(self.steps[0])
            
        else:
            # Show an error
            print("Error: %s file not found" % fname)        

    def read_dat(self):
        # reads the element results file
        fname = self.fname+'.dat'
        if os.path.isfile(fname):
            
            f = open(fname)
            mode = 'off'
            time = 0.0
            while True:
                line = f.readline()
                if not line:
                    break
                line = line.strip()
                
                # check for read flags for displ or stress
                if 'stress' in line:
                    words = line.split()
                    # add time if not present
                    time = float(words[-1])
                    if time not in self.steps:
                        self.steps.append(time)
                    if time not in self.results:
                        self.results[time] = {'node':{},'element':{}}
                    # set mode
                    if 'stress' in line:
                        mode = 'stress'
                    f.readline()
                    line = f.readline().strip()
                
                # reset the read type if we hit a blank line
                if line == '':
                    mode = 'off'
                
                if mode == 'stress':
                    # store stress results
                    w = line.split()
                    labs = ['Sx', 'Sy', 'Sz', 'Sxy', 'Sxz', 'Syz', 'Seqv']
                    vals = [float(a) for a in w[2:]]
                    seqv = self.seqv(vals)
                    [s1, s2, s3] = self.principals(vals)
                    vals.append(seqv)
                    vals += [s1, s2, s3]
                    #print(vals)
                    enum = int(w[0])
                    ipnum = int(w[1]) # integration point number
                    d = {}
                    for (key, val) in zip(labs, vals):
                        d[key] = val
                    if enum not in self.results[time]['element']:
                        self.results[time]['element'][enum] = {}
                    # each line is an integration point result
                    self.results[time]['element'][enum][ipnum] = d

            f.close()
            
            # loop over all element results, calculating avg element result
            # by averaging integration point vals
            for t in self.steps:
                for (enum, edict) in self.results[t]['element'].items():
                    ipnums = [a for a in edict.keys() if isinstance(a, int)]
                    ipts = [edict[a] for a in ipnums]
                    
                    labs = ['Sx', 'Sy', 'Sz', 'Sxy', 'Sxz', 'Syz']
                    vals = []
                    for (ind, field) in enumerate(labs):
                        vlist = [a[field] for a in ipts]
                        val = sum(vlist)/len(vlist)
                        edict[field] = val
                        vals.append(val)
                    # for each element, caclulate Seqv, S1, S2, S3
                    seqv = self.seqv(vals)
                    [s1, s2, s3] = self.principals(vals)
                    edict['Seqv'] = seqv
                    edict['S1'] = s1
                    edict['S2'] = s2
                    edict['S3'] = s3
            
            print('The following times have been read: ',self.steps)
            print('Results from file: %s have been read.' % (fname))
            
        else:
            # Show an error
            print("Error: %s file not found" % fname)        
