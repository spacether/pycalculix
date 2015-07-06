"""This module stores the Results_File class."""

import re # used to get info from frd file
import math # used for metric number conversion
from numpy.lib.polynomial import roots # need to find S1-S3
from numpy.core.function_base import linspace # need to make contours
import os #need to check if results file exists
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np

from . import base_classes # needed for RESFIELDS
from . import mesh

CMAP = 'jet'

class ResultsFile(object):
    """Makes a results file.

    Args:
        problem (Problem): problem that was solved

    Attributes:
        __problem (Problem): parent problem
        __steps (list): a list of float time steps
        __results (dict): a dict storing the results data
            results[step]['node'][nnum][field] --> value
            results[step]['element'][enum]['avg'][field] --> value
            results[step]['element'][enum]['max'][field] --> value
            results[step]['element'][enum]['min'][field] --> value
            results[step]['element'][enum]['ipoints'][ipnum][field] --> value
            field = 'ux' or 'Seqv' or 'ey' etc.
        __time (float): current time we are looking at, defaults to -1
            When a file is loaded in, the first time step is loaded.
    """

    def __init__(self, problem):
        self.__problem = problem
        self.__steps = [] # this stores a list of time steps
        self.__results = {} # stores results, nested dicts
        self.__time = -1
        if self.__problem.solved:
            self.load()

    @property
    def steps(self):
        """Returns a list of loaded time steps.

        Note: this is read only, you can not assign a value to it.
        """
        return self.__steps

    @property
    def time(self):
        """Returns the current time (float) in the results file.

        Note: this is read only, you can not assign a value to it.
        """
        return self.__time

    @staticmethod
    def __metric_num(number, sig_figs=3, sci=False):
        """Returns string of number, with only 10**3 suffixes.

        If 0 <= number < 1000 no suffix is added.
        This is useful for quantities usually given in metric: stress, displacement.

        Args:
            number (float or int): the number we want converted to __metric_number
            sig_figs (int): number of significant figures to use, right of decimal
            sci (bool): True means use scientific formatting, False use metric
        """
        if sci:
            format_str = "%.{}e".format(sig_figs)
            my_str = format_str % number
        else:
            format_str = "%.{}f".format(sig_figs)
            my_str = format_str % number
            if number != 0:
                # get the scientific exponent
                exp = math.floor(math.log10(abs(number)))
                metric_exp = exp - (exp % 3)
                new_float = number/(10**metric_exp)
                if metric_exp != 0:
                    format_str = "%.{}fe%i".format(sig_figs)
                    my_str = format_str % (new_float, metric_exp)
        return my_str

    def load(self):
        """Loads the results file with problem.fname prefix."""
        self.__read_frd() # read nodal results
        self.__read_dat() # read element integration pt results

    def nplot(self, field, fname='', display=True, levels=21, gradient=False,
              gmult=1.0, max_val=None, min_val=None, title=''):
        """Plots nodal results.

        Args:
            field (str): results item to plot, examples: 'ux', 'ey', 'Seqv'
            fname (str): prefix of png file name, if writing an image
            display (bool): True = interactively show the plot
            levels (int): number of levels to use in the colorbar
            gradient (bool): True = results plotted with gradient
                False = results plotted with filled areas
            gmult (int): geometric multiplier on displacement of nodes
                displayed_node_loc = model_node_loc + gmult*node_displacement
            max_val (float or None): max value in the colorbar

                - None: max from selected data used
                - float: use the passed float
            min_val (float): min value in the colorbar

                - None: min from selected data used
                - float: use the passed float
            title (str): third line in the plot title
        """
        # store the selected nodes and elements
        sel = {}
        sel['nodes'] = self.__problem.fea.view.nodes
        sel['elements'] = self.__problem.fea.view.elements
        sel['faces'] = self.__problem.fea.view.faces

        # sort nodes low to high so index is correct
        # we have index to id below so showing subsets works
        sel['nodes'] = list(sel['nodes'])
        sel['nodes'] = sorted(sel['nodes'], key=lambda k: k.id)

        # store results at nodes
        axials = []
        radials = []
        zvals = []
        id_to_ind = {}
        for node in sel['nodes']:
            id_to_ind[node.id] = len(axials)
            axi = node.y + gmult*self.__results[self.__time]['node'][node.id]['uy']
            rad = node.x + gmult*self.__results[self.__time]['node'][node.id]['ux']
            axials.append(axi)
            radials.append(rad)
            zvals.append(self.__results[self.__time]['node'][node.id][field])

        # make a list of triangles, given by indices, looping anticlockwise
        triangles = []
        mylist = []
        if len(sel['elements']) > 0:
            mylist = sel['elements']
        elif len(sel['faces']) > 0:
            mylist = sel['faces']
        for element in mylist:
            tris = element.get_tris() # list of triangle nodes
            for tri in tris:
                for ind, nid in enumerate(tri):
                    tri[ind] = id_to_ind[nid]     # convert id to index
            triangles += tris

        # check to see if selected nodes and elements are
        # in the parent model's nodes and elements
        fig = plt.figure()
        ax_ = fig.add_subplot(111)

        # need to set tick list here
        vmin = min(zvals)
        vmax = max(zvals)
        stop_plot = False
        if max_val != None and min_val == None:
            if max_val < vmin:
                stop_plot = True
                print('Error:')
                print(' Only max was passed but it is < the data min!')
                print(' Pass a max_val that is > the data min of %f' % vmin)
            else:
                vmax = max_val
        elif min_val != None and max_val == None:
            if min_val > vmax:
                stop_plot = True
                print('Error:')
                print(' Only min was passed but it is > the data max!')
                print(' Pass a min_val that is < the data max of %f' % vmax)
            else:
                vmin = min_val
        elif max_val != None and min_val != None:
            if max_val < min_val:
                stop_plot = True
                print('Error:')
                print(' Min and max passed, but min > max!')
                print(' Pass a min_val that is < max_val')
            else:
                vmax = max_val
                vmin = min_val
        # exit if stop plot flag is on
        if stop_plot:
            return None

        tick_list = [vmin]
        if vmax != vmin:
            # we have a range of values we're plotting
            tick_list = linspace(vmin, vmax, levels+1)

        # plot using a gradient(shaded) or levels
        # code required for the colorbar, needs to go before plotting for colormap
        cnorm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
        if vmax != vmin:
            # we have a range of values we're plotting
            if gradient:
                cmap = plt.get_cmap(CMAP)
            else:
                cmap = plt.get_cmap('jet', levels)
        cmap.set_under('0.3', 0.8)
        cmap.set_over('0.7', 0.8)
        if gradient or len(tick_list) == 1:
            # This one is shaded
            plt.tripcolor(axials, radials, triangles, zvals, shading='gouraud',
                          cmap=cmap, norm=cnorm)
        else:
            # this one is not shaded
            plt.tricontourf(axials, radials, triangles, zvals, levels=tick_list,
                            cmap=cmap, norm=cnorm, extend='both')

        scalarmap = cmx.ScalarMappable(norm=cnorm, cmap=cmap)
        scalarmap.set_array([])
        cbar = plt.colorbar(scalarmap, orientation='vertical', ticks=tick_list)

        scibool = False
        if field[0] == 'e':
            # strain plotting, use scientific numbering
            scibool = True
        met_max = self.__metric_num(max(zvals), sci=scibool)
        met_min = self.__metric_num(min(zvals), sci=scibool)
        label = 'Max: %s\nMin: %s' % (met_max, met_min)
        tick_list = [self.__metric_num(tick, sci=scibool) for tick in tick_list]
        cbar.ax.set_yticklabels(tick_list)
        cbar.ax.set_xlabel(label, labelpad=10, x=0, ha='left')
        cbar.ax.xaxis.set_label_position('top')

        # set the horizontal and vertical axes
        base_classes.plot_set_bounds(plt, axials, radials)

        # set units
        alist = self.__problem.fea.get_units(field, 'dist', 'time')
        [f_unit, d_unit, t_unit] = alist

        # set plot axes
        plot_title = ('Node %s%s\nTime=%f%s' %
                      (field, f_unit, self.__time, t_unit))
        if title != '':
            plot_title += '\n%s' % title
        plt.title(plot_title)
        plt.xlabel('axial, y'+d_unit)
        plt.ylabel('radial, x'+d_unit)
        ax_.set_aspect('equal')
        if gmult != 1:
            ax_.xaxis.set_ticklabels([])
            ax_.yaxis.set_ticklabels([])
        base_classes.plot_finish(plt, fname, display)

    def eplot(self, field, fname='', display=True, levels=21,
              gmult=1.0, mode='avg', max_val=None, min_val=None, title=''):
        """Plots element results.

        Args:
            field (str): results item to plot. Only stresses supported.
                Examples: 'Sx', 'Sxy', 'S1', 'Seqv' etc.
            fname (str): prefix of png file name, if writing an image
            display (bool): True = interactively show the plot
            levels (int): number of levels to use in the colorbar
            gmult (int): geometric multiplier on displacement of nodes
                displayed_node_loc = model_node_loc + gmult*node_displacement
            mode (str): the type of element result to plot

                - 'avg': integration points averaged to avg element result
                - 'max': max value of field in the integration points plotted
                - 'min': min value of field in the integration points plotted
            max_val (float or None): max value in the colorbar

                - None: max from selected data used
                - float: use the passed float
            min_val (float): min value in the colorbar

                - None: min from selected data used
                - float: use the passed float
            title (str): third line in the plot title
        """
        # store the selected nodes and elements
        sel = {}
        sel['nodes'] = self.__problem.fea.view.nodes
        sel['elements'] = self.__problem.fea.view.elements
        sel['faces'] = self.__problem.fea.view.faces

        # sort nodes low to high so index is correct
        # we have index to id below so showing subsets works
        sel['nodes'] = list(sel['nodes'])
        sel['nodes'] = sorted(sel['nodes'], key=lambda k: k.id)

        # store results at nodes
        axials = []
        radials = []
        zvals = []
        id_to_ind = {}
        for node in sel['nodes']:
            id_to_ind[node.id] = len(axials)
            axi = node.y + gmult*self.__results[self.__time]['node'][node.id]['uy']
            rad = node.x + gmult*self.__results[self.__time]['node'][node.id]['ux']
            axials.append(axi)
            radials.append(rad)

        # make a list of triangles, given by indices, looping anticlockwise
        triangles = []
        mylist = []
        if len(sel['elements']) > 0:
            mylist = sel['elements']
        elif len(sel['faces']) > 0:
            mylist = sel['faces']
        for ele in mylist:
            val = self.__results[self.__time]['element'][ele.id][mode][field]
            tris = ele.get_tris() # list of triangle nodes defined by node id
            for tri in tris:
                zvals.append(val)
                for ind, nid in enumerate(tri):
                    tri[ind] = id_to_ind[nid]     # convert id to index
            triangles += tris

        # check to see if selected nodes and elements are
        # in the parent model's nodes and elements

        fig = plt.figure()
        ax_ = fig.add_subplot(111)

        # need to set tick list here
        vmin = min(zvals)
        vmax = max(zvals)
        stop_plot = False
        if max_val != None and min_val == None:
            if max_val < vmin:
                stop_plot = True
                print('Error:')
                print(' Only max was passed but it is < the data min!')
                print(' Pass a max_val that is > the data min of %f' % vmin)
            else:
                vmax = max_val
        elif min_val != None and max_val == None:
            if min_val > vmax:
                stop_plot = True
                print('Error:')
                print(' Only min was passed but it is > the data max!')
                print(' Pass a min_val that is < the data max of %f' % vmax)
            else:
                vmin = min_val
        elif max_val != None and min_val != None:
            if max_val < min_val:
                stop_plot = True
                print('Error:')
                print(' Min and max passed, but min > max!')
                print(' Pass a min_val that is < max_val')
            else:
                vmax = max_val
                vmin = min_val
        # exit if stop plot flag is on
        if stop_plot:
            return None

        tick_list = [vmin]
        if vmax != vmin:
            # we have a range of values we're plotting
            tick_list = linspace(vmin, vmax, levels+1)

        # code required for the colorbar, needs to go before plotting for cmap
        cnorm = colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = colors.ListedColormap(['b', 'b']) # default to plot one val
        if vmax != vmin:
            # we have a range of values we're plotting
            cmap = plt.get_cmap(CMAP, levels)
        cmap.set_under('0.3', 0.8)
        cmap.set_over('0.7', 0.8)

        # plot using levels
        plt.tripcolor(axials, radials, triangles, zvals,
                      shading='flat', cmap=cmap, norm=cnorm)

        scalarmap = cmx.ScalarMappable(norm=cnorm, cmap=cmap)
        scalarmap.set_array([])
        cbar = plt.colorbar(scalarmap, orientation='vertical', ticks=tick_list)
        scibool = False
        if field[0] == 'e':
            # strain plotting, use scientific numbering
            scibool = True
        met_max = self.__metric_num(max(zvals), sci=scibool)
        met_min = self.__metric_num(min(zvals), sci=scibool)
        label = 'Max: %s\nMin: %s' % (met_max, met_min)
        tick_list = [self.__metric_num(tick, sci=scibool) for tick in tick_list]
        cbar.ax.set_yticklabels(tick_list)
        cbar.ax.set_xlabel(label, labelpad=10, x=0, ha='left')
        cbar.ax.xaxis.set_label_position('top')

        # set the horizontal and vertical axes
        base_classes.plot_set_bounds(plt, axials, radials)

        # set units
        alist = self.__problem.fea.get_units(field, 'dist', 'time')
        [f_unit, d_unit, t_unit] = alist

        # set plot axes
        plot_title = ('Element %s %s%s\nTime=%f%s' %
                      (mode, field, f_unit, self.__time, t_unit))
        if title != '':
            plot_title += '\n%s' % title
        plt.title(plot_title)
        plt.xlabel('axial, y'+d_unit)
        plt.ylabel('radial, x'+d_unit)
        ax_.set_aspect('equal')
        if gmult != 1:
            ax_.xaxis.set_ticklabels([])
            ax_.yaxis.set_ticklabels([])
        base_classes.plot_finish(plt, fname, display)

    def set_time(self, time):
        """Sets the time point we're looking at in the results file.

        Args:
            time (float): time we are setting
        """
        if time in self.steps:
            self.__time = time
            print('Results file time set to: %f' % (self.__time))
        else:
            print('Time %f is not in the loaded times. Valid times are:')
            print(self.steps)
    
    
    def plot_gradient(self, start_point, end_point, field, fname='', display=True, title='', max_val=None, min_val=None, curve_fitting=True, n_poly=3, n_subpoints=500, legend=True):
        """Create diagram with data projected onto line on the undeformed geometry.

        Args:
            start_point [(float), (float)]: starting point of line. [x, y]
            end_point [(float), (float)]: end point of line. Example: [x, y]
            field (str): results item to plot, examples: 'ux', 'ey', 'Seqv'
            
        Kargs:
            fname (str): prefix of png file name, if writing an image
            display (bool): True = interactively show the plot
            title (str): third line in the plot title
            max_val (float or None): max value in the y-axis
                - None: max from selected data used
                - float: use the passed float
            min_val (float or None): min value in the y-axis
                - None: min from selected data used
                - float: use the passed float
            curve_fitting (bool): True = a curve is fitted to the gradient
            n_poly (int): numbers of polygons for fitting
            n_subpoints (int): numbers of points the line is subdivided into
            legend (bool): True = legend with fitted equation is shown
        """
        
        # store the selected nodes and elements
        sel = {}
        sel['nodes'] = self.__problem.fea.view.nodes

        # sort nodes low to high so index is correct
        # we have index to id below so showing subsets works
        sel['nodes'] = list(sel['nodes'])
        sel['nodes'] = sorted(sel['nodes'], key=lambda k: k.id)

        # store results at nodes
        node_position = np.zeros((len(sel['nodes']),2))
        field_values = np.zeros(len(sel['nodes']))
        
        for idx, node in enumerate(sel['nodes']):
            
            node_position[idx] = [node.x, node.y]
            field_values[idx] = self.__results[self.__time]['node'][node.id][field]
            
        
        #create subpoints on line
        subpoints = np.zeros((n_subpoints, 3))  #[x, y, line position]
        
        subpoints[:,0] = np.linspace(start_point[0], end_point[0], n_subpoints)
        subpoints[:,1] = np.linspace(start_point[1], end_point[1], n_subpoints)
        subpoints[:,2] = np.arange(n_subpoints) / n_subpoints * np.sqrt(np.sum( (np.array(start_point) - np.array(end_point))**2))
        
        #calculate weighted field value for every subpoint
        wfield = np.zeros(n_subpoints)
        
        for idx in range(n_subpoints):
            
            #calculate inverse of distance from nodes to subpoints
            dist = np.sqrt(np.sum((node_position-subpoints[idx,0:2])**2,axis=1))
            
            #calculte weighted field value
            #dist[dist < 1E-10] = 1E-10
            #inv_dist = 1. / dist**3
            #wfield[idx] = np.average(field_values, weights=inv_dist)
            
            #use nearest value
            wfield[idx] = field_values[min(range(len(dist)),key=dist.__getitem__)]
            
            
        #plot diagram
        
        fig = plt.figure(figsize=(10,6))
        ax_ = fig.add_subplot(111)
        
        plt.plot(subpoints[:,2], wfield, '-r', linewidth=2.5, label=field)
        
        if curve_fitting==True:
            #execute curve fitting if needed
            poly = np.polyfit(subpoints[:,2], wfield, n_poly)
            
            #string for equation of fitted function
            funcstring = [str(np.round(poly[i]))+u'*x^'+str(np.arange(n_poly,0,-1)[i]) for i in range(n_poly)]
            funcstring.append(str(np.round(poly[-1])))
            funcstring = '+'.join(funcstring)
            
            func = np.poly1d(poly)
            
            plt.plot(subpoints[:,2], func(subpoints[:,2]), '--k', linewidth=1.5, label=funcstring)
        
        
        # set units
        alist = self.__problem.fea.get_units(field, 'dist', 'time')
        [f_unit, d_unit, t_unit] = alist

        # set plot axes
        plot_title = ('Gradient %s%s\nTime=%f%s' %(field, f_unit, self.__time, t_unit))
        if title != '':
            plot_title += '\n%s' % title
        plt.title(plot_title)
        plt.xlabel('path position'+d_unit)
        plt.ylabel(field + ' ' +f_unit)
        
        #show legend if needed
        if legend == True:
            plt.legend()
            
        #set limits on y-axis
        if min_val!=None:
            plt.gca().set_ylim(bottom=min_val)
        if max_val!=None:
            plt.gca().set_ylim(top=max_val)
            
        plt.grid()
        base_classes.plot_finish(plt, fname, display)
    
    def get_relative_gradient(self, start_point, end_point, field, n_poly=3, n_subpoints=500):
        """Calculte relative stress gradient (gradient/start_value)

        Args:
            start_point [(float), (float)]: starting point of line. [x, y]
            end_point [(float), (float)]: end point of line. Example: [x, y]
            field (str): results item to plot, examples: 'ux', 'ey', 'Seqv'
            
        Kargs:
            n_poly (int): numbers of polygons for fitting, min=2
            n_subpoints (int): numbers of points the line is subdivided into
        """
        
        # store the selected nodes and elements
        sel = {}
        sel['nodes'] = self.__problem.fea.view.nodes

        # sort nodes low to high so index is correct
        # we have index to id below so showing subsets works
        sel['nodes'] = list(sel['nodes'])
        sel['nodes'] = sorted(sel['nodes'], key=lambda k: k.id)

        # store results at nodes
        node_position = np.zeros((len(sel['nodes']),2))
        field_values = np.zeros(len(sel['nodes']))
        
        for idx, node in enumerate(sel['nodes']):
            
            node_position[idx] = [node.x, node.y]
            field_values[idx] = self.__results[self.__time]['node'][node.id][field]
            
        
        #create subpoints on line
        subpoints = np.zeros((n_subpoints, 3))  #[x, y, line position]
        
        subpoints[:,0] = np.linspace(start_point[0], end_point[0], n_subpoints)
        subpoints[:,1] = np.linspace(start_point[1], end_point[1], n_subpoints)
        subpoints[:,2] = np.arange(n_subpoints) / n_subpoints * np.sqrt(np.sum( (np.array(start_point) - np.array(end_point))**2))
        
        #calculate weighted field value for every subpoint
        wfield = np.zeros(n_subpoints)
        
        for idx in range(n_subpoints):
            
            #calculate inverse of distance from nodes to subpoints
            dist = np.sqrt(np.sum((node_position-subpoints[idx,0:2])**2,axis=1))
            
            #use nearest value
            wfield[idx] = field_values[min(range(len(dist)),key=dist.__getitem__)]
            
            
        #curve fitting
        poly = np.polyfit(subpoints[:,2], wfield, n_poly)
        
        rel_grad = abs(poly[-2])/abs(poly[-1])
        
        return rel_grad
        
        
        
        

    @staticmethod
    def __utot(vals):
        """Returns the total displacement distance, given [dx,dy,dz].

        Args:
            vals (list): [dx, dy, dz] list of displacements in x, y, and z axes

        Returns:
            res (float): displacement
        """
        # computes sum of the squares
        res = [a**2 for a in vals]
        res = (sum(res))**0.5
        return res

    @staticmethod
    def __seqv(vals):
        """Returns the Von Mises stress, which will be stored as 'Seqv'.

        Args:
            vals (list): list of six stresses [s11,s22,s33,s12,s13,s23]

        Returns:
            res (float): Von Mises stress
        """
        [s11, s22, s33, s12, s13, s23] = vals
        aval = s11 - s22
        bval = s22 - s33
        cval = s33 - s11
        dval = s12**2 + s23**2 +s13**2
        res = (0.5*(aval**2 + bval**2 + cval**2 +6*dval))**0.5
        return res

    @staticmethod
    def __principals(vals):
        """Returns principal stresses [S1,S2,S3].

        Args:
            vals (list): six stresses [s11,s22,s33,s12,s13,s23]

        Returns:
            res (list): principal stresses [S1,S2,S3] stresses are high-to-low
        """
        # calculates and returns principal stresses, S1, S2, S3
        [s11, s22, s33, s12, s13, s23] = vals
        aval = 1
        bval = (s11 + s22 + s33)*-1.0
        cval = (s11*s22 + s11*s33 + s22*s33 - s12**2 - s13**2 - s23**2)
        dval = (s11*s22*s33 + 2*s12*s13*s23 - s11*(s23**2) - s22*(s13**2)
                - s33*(s12**2))*-1.0
        res = list(roots([aval, bval, cval, dval]))
        res = sorted(res, reverse=True)
        return res

    def __get_data_dict(self, time, type_str):
        """Returns the data dict at the correct time for element or node.

        Args:
            time (float): None or the time we want, if None use current time
            type_str: 'element' or 'node'

        Returns:
            res (dict or None): dictionary with field values in it
                None if the time was invalid
        """
        res = self.__results[self.__time][type_str]
        if time != None:
            if time not in self.steps:
                print('Error: passed time is not in steps!')
                print(' Pass a time in the steps:')
                print(self.steps)
                return None
            else:
                res = self.__results[time][type_str]
        return res

    def get_nmax(self, field, time=None):
        """Returns the max value of node results field in selected nodes.

        Reports results for the current time.

        Args:
            field (str): results field, for example 'ux', 'ey', 'S1', 'Seqv'
            time (None or float): the time to query

                - None: uses the current time
                - float: uses the passed float time
        Returns:
            res (float): max value
        """
        nodes = self.__problem.fea.view.nodes
        node_ids = [node.id for node in nodes]
        data_dict = self.__get_data_dict(time, 'node')
        if data_dict == None:
            return None
        ndicts = [data_dict[nid] for nid in node_ids]
        res = [ndict[field] for ndict in ndicts]
        res = max(res)
        return res

    def get_nmin(self, field, time=None):
        """Returns the min value of node results field in selected nodes.

        Reports results for the current time.

        Args:
            field (str): results field, for example 'ux', 'ey', 'S1', 'Seqv'
            time (None or float): the time to query

                - None: uses the current time
                - float: uses the passed float time
        Returns:
            res (float): min value
        """
        nodes = self.__problem.fea.view.nodes
        node_ids = [node.id for node in nodes]
        data_dict = self.__get_data_dict(time, 'node')
        if data_dict == None:
            return None
        ndicts = [data_dict[nid] for nid in node_ids]
        res = [ndict[field] for ndict in ndicts]
        res = min(res)
        return res

    def get_nval(self, node, field, time=None):
        """Returns the field result value under node.

        Result will be returned whether or not passed node is selected.

        Args:
            node (str or Node): node we are asking about
            field (str): the results item we want: 'ux', 'Sy', 'Seqv', 'fx'
            time (None or float): the time to query

                - None: uses the current time
                - float: uses the passed float time
        Returns:
            res (float or None): float value if field exists, None otherwise
        """
        items = self.__problem.fea.get_item(node)
        if len(items) == 1:
            if isinstance(items[0], mesh.Node):
                nnum = items[0].id
                data_dict = self.__get_data_dict(time, 'node')
                if data_dict == None:
                    return None
                ndict = data_dict[nnum]
                if field in ndict:
                    res = ndict[field]
                    return res
                else:
                    print('Passed field is not in the results!')
                    return None
            else:
                print('You did not pass in a node!')
                print('A single node or string node name must be passed in!')
                return None
        else:
            print('A single node or string node name must be passed in!')
            return None

    def get_fsum(self, item):
        """Returns the force sum on nodes under a given point or line.

        Reports results for the current time.

        Args:
            item (Point or SignLine): item that has reaction forces on its nodes

        Returns:
            list: [fx, fy, fz] reaction forces in each axis, force units
        """
        (fxx, fyy, fzz) = ([], [], [])
        nodes = item.nodes
        nodes = [n.id for n in nodes]
        for node in nodes:
            f_x = self.__results[self.__time]['node'][node]['fx']
            f_y = self.__results[self.__time]['node'][node]['fy']
            f_z = self.__results[self.__time]['node'][node]['fz']
            if f_x != 0 or f_y != 0 or f_z != 0:
                fxx.append(f_x)
                fyy.append(f_y)
                fzz.append(f_z)
        fxx = sum(fxx)
        fyy = sum(fyy)
        fzz = sum(fzz)
        return [fxx, fyy, fzz]

    def get_emax(self, field, time=None, mode='avg'):
        """Returns the max results field value of selected elements at curent time.

        Args:
            field (str): results field, stresses supported 'S1', 'Sx', etc.
            time (None or float): the time to query

                - None: uses the current time
                - float: uses the passed float time
            mode (str): type of element result to give back

                - 'max': for each element only use the max value of field over
                    all of its integration points
                - 'min': for each element only use the min value of field over
                    all of its integration points
                - 'avg': for each element only use an average of all integration
                    points in the eleemnt. Principal streses and Seqv are
                    calculated after averaging 6 stress components.
        Returns:
            res (float): max value
        """
        res = []
        elements = self.__problem.fea.view.elements
        data_dict = self.__get_data_dict(time, 'element')
        if data_dict == None:
            return None
        for element in elements:
            enum = element.id
            edict = data_dict[enum][mode]
            res.append(edict[field])
        res = max(res)
        return res

    def get_emin(self, field, time=None, mode='avg'):
        """Returns the min results field value of selected elements at curent time.

        Args:
            field (str): results field, stresses supported 'S1', 'Sx', etc.
            time (None or float): the time to query

                - None: uses the current time
                - float: uses the passed float time
            mode (str): type of element result to give back

                - 'max': for each element only use the max value of field over
                    all of its integration points
                - 'min': for each element only use the min value of field over
                    all of its integration points
                - 'avg': for each element only use an average of all integration
                    points in the eleemnt. Principal streses and Seqv are
                    calculated after averaging 6 stress components.
        Returns:
            res (float): min value
        """
        res = []
        elements = self.__problem.fea.view.elements
        data_dict = self.__get_data_dict(time, 'element')
        if data_dict == None:
            return None
        for element in elements:
            enum = element.id
            edict = data_dict[enum][mode]
            res.append(edict[field])
        res = min(res)
        return res

    def get_eval(self, element, field, time=None, mode='avg'):
        """Returns the field result value under element.

        Result will be returned whether or not passed element is selected.

        Args:
            element (str or Element): element we are asking about
            field (str): the results item we want: 'Sy', 'Seqv'
            mode (str): the type of element result to get

                - 'avg': integration points averaged to avg element result
                - 'max': max value of field in the integration points plotted
                - 'min': min value of field in the integration points plotted
        Returns:
            res (float or None): float value if field exists, None otherwise
        """
        items = self.__problem.fea.get_item(element)
        if len(items) == 1:
            if isinstance(items[0], mesh.Element):
                enum = items[0].id
                data_dict = self.__get_data_dict(time, 'element')
                if data_dict == None:
                    return None
                edict = data_dict[enum][mode]
                if field in edict:
                    res = edict[field]
                    return res
                else:
                    print('Passed field is not in the results!')
                    return None
            else:
                print('You did not pass in a element!')
                print('A single element or string element name must be given!')
                return None
        else:
            print('A single element or string element name must be given!')
            return None

    @staticmethod
    def __get_vals(fstr, line):
        """Returns a list of typed items based on an input format string.

        Args:
            fst (str): C format string, commas separate fields
            line (str): line string to parse

        Returns:
            res (list): list of typed items extracted from the line
        """
        res = []
        fstr = fstr.split(',')
        thestr = str(line)
        for item in fstr:
            if item[0] == "'":
                # strip off the char quaotes
                item = item[1:-1]
                # this is a string entry, grab the val out of the line
                ind = len(item)
                fwd = thestr[:ind]
                thestr = thestr[ind:]
                res.append(fwd)
            else:
                # format is: 1X, A66, 5E12.5, I12
                # 1X is number of spaces
                (mult, ctype) = (1, None)
                m_pat = re.compile(r'^\d+') # find multiplier
                c_pat = re.compile(r'[XIEA]') # find character
                if m_pat.findall(item) != []:
                    mult = int(m_pat.findall(item)[0])
                ctype = c_pat.findall(item)[0]
                if ctype == 'X':
                    # we are dealing with spaces, just reduce the line size
                    thestr = thestr[mult:]
                elif ctype == 'A':
                    # character string only, add it to results
                    fwd = thestr[:mult].strip()
                    thestr = thestr[mult:]
                    res.append(fwd)
                else:
                    # IE, split line into m pieces
                    w_pat = re.compile(r'[IE](\d+)') # find the num after char
                    width = int(w_pat.findall(item)[0])
                    while mult > 0:
                        # only add items if we have enough line to look at
                        if width <= len(thestr):
                            substr = thestr[:width]
                            thestr = thestr[width:]
                            substr = substr.strip() # remove space padding
                            if ctype == 'I':
                                substr = int(substr)
                            elif ctype == 'E':
                                substr = float(substr)
                            res.append(substr)
                        mult -= 1
        return res

    @staticmethod
    def __skip_lines(infile, numlines):
        """Reads forward numlines from infile"""
        lines = int(numlines)
        while lines > 0:
            infile.readline()
            lines -= 1

    def __read_frd(self):
        """Reads a ccx frd results file which contains nodal results."""
        fname = self.__problem.fname+'.frd'
        if os.path.isfile(fname):

            print('Reading results file: '+fname)
            infile = open(fname, 'r')
            mode = 'off'
            time = 0.0
            format_ = None
            rfstr = ''
            node_ids = []
            while True:
                line = infile.readline()
                if not line:
                    break

                #-------------
                # set the modes
                #--------------
                if '2C' in line:
                    # we are in nodal definition
                    fstr = "1X,'   2','C',18X,I12,37X,I1"
                    # [key, code, numnod, format_]
                    format_ = self.__get_vals(fstr, line)[3]

                    # set results format to short, long or binary
                    # only short and long are parsed so far
                    if format_ == 0:
                        rfstr = "1X,'-1',I5,3E12.5"
                    elif format_ == 1:
                        rfstr = "1X,'-1',I10,3E12.5"
                    elif format_ == 2:
                        # binary
                        pass

                    line = infile.readline()
                    mode = 'nodes'
                    print('Reading '+mode)

                elif '1PSTEP' in line:
                    # we are in a results block

                    # next line has the time in it
                    line = infile.readline()
                    fstr = "1X,' 100','C',6A1,E12.5,I12,20A1,I2,I5,10A1,I2"
                    tmp = self.__get_vals(fstr, line)
                    #[key, code, setname, value, numnod, text, ictype, numstp, analys, format_]
                    value, format_ = tmp[3], tmp[9]

                    # set results format to short, long or binary
                    # only short and long are parsed so far
                    if format_ == 0:
                        rfstr = "1X,I2,I5,6E12.5"
                    elif format_ == 1:
                        rfstr = "1X,I2,I10,6E12.5"
                    elif format_ == 2:
                        # binary
                        pass

                    # set the time
                    time = value
                    if time not in self.__steps:
                        self.__steps.append(time)
                    if time not in self.__results:
                        self.__results[time] = {'node':{}, 'element':{}}
                        for nid in node_ids:
                            # make dict for each node
                            self.__results[time]['node'][nid] = {}

                    # get the name to determine if stress or displ
                    line = infile.readline()
                    fstr = "1X,I2,2X,8A1,2I5"
                    # [key, name, ncomps, irtype]
                    name = self.__get_vals(fstr, line)[1]

                    if name == 'DISPR':
                        mode = 'displ'
                        self.__skip_lines(infile, 4)
                    elif name == 'STRESSR':
                        mode = 'stress'
                        self.__skip_lines(infile, 6)
                    elif name == 'TOSTRAIR':
                        mode = 'strain'
                        self.__skip_lines(infile, 6)
                    elif name == 'FORCR':
                        mode = 'force'
                        self.__skip_lines(infile, 4)

                    print('Reading '+mode+' storing: '+
                          ','.join(base_classes.RESFIELDS[mode]))
                    line = infile.readline()

                #----------------------------
                # Store results
                #----------------------------

                # reset the read mode if we hit an end of record
                if line[:3] == ' -3':
                    mode = 'off'

                if mode == 'nodes':
                    # node definition, store node numbers only
                    # [key, node, x, y, z]
                    node = self.__get_vals(rfstr, line)[1]
                    node_ids.append(node)

                elif mode == 'displ':
                    # displacements
                    # [key, node, ux, uy, uz]
                    [node, ux_, uy_, uz_] = self.__get_vals(rfstr, line)[1:]
                    labs = base_classes.RESFIELDS[mode]
                    vals = [ux_, uy_, uz_]
                    utot = self.__utot(vals)
                    vals.append(utot)
                    adict = self.__results[time]['node'][node]
                    for (label, val) in zip(labs, vals):
                        adict[label] = val

                elif mode == 'stress':
                    # stresses
                    tmp = self.__get_vals(rfstr, line)
                    # [key, node, sx, sy, sz, sxy, syz, szx]
                    [node, sxx, syy, szz, sxy, syz, szx] = tmp[1:]
                    labs = base_classes.RESFIELDS[mode]
                    vals = [sxx, syy, szz, sxy, syz, szx]
                    seqv = self.__seqv(vals)
                    [s_1, s_2, s_3] = self.__principals(vals)
                    vals.append(seqv)
                    vals += [s_1, s_2, s_3]
                    adict = self.__results[time]['node'][node]
                    for (label, val) in zip(labs, vals):
                        adict[label] = val

                elif mode == 'strain':
                    # strains
                    tmp = self.__get_vals(rfstr, line)
                    # [key, node, ex, ey, ez, exy, eyz, ezx]
                    [node, exx, eyy, ezz, exy, eyz, ezx] = tmp[1:]
                    labs = base_classes.RESFIELDS[mode]
                    vals = [exx, eyy, ezz, exy, eyz, ezx]
                    eeqv = self.__seqv(vals)
                    [e_1, e_2, e_3] = self.__principals(vals)
                    vals.append(eeqv)
                    vals += [e_1, e_2, e_3]
                    adict = self.__results[time]['node'][node]
                    for (label, val) in zip(labs, vals):
                        adict[label] = val

                elif mode == 'force':
                    # reaction forces
                    # [key, node, fx, fy, fz]
                    [node, f_x, f_y, f_z] = self.__get_vals(rfstr, line)[1:]
                    labs = base_classes.RESFIELDS[mode]
                    vals = [f_x, f_y, f_z]
                    adict = self.__results[time]['node'][node]
                    for (label, val) in zip(labs, vals):
                        adict[label] = val

            infile.close()
            print('The following times have been read:')
            print(self.__steps)
            print('Done reading file: %s' % fname)
            self.set_time(self.__steps[0])

        else:
            print("Error: %s file not found" % fname)

    def __read_dat(self):
        """Reads ccx dat results file. It has element integration point results.
        """
        fname = self.__problem.fname+'.dat'
        if os.path.isfile(fname):

            infile = open(fname, 'r')
            mode = 'off'
            time = 0.0
            while True:
                line = infile.readline()
                if not line:
                    break
                line = line.strip()

                # check for read flags for displ or stress
                if 'stress' in line:
                    words = line.split()
                    # add time if not present
                    time = float(words[-1])
                    if time not in self.__steps:
                        self.__steps.append(time)
                    if time not in self.__results:
                        self.__results[time] = {'node':{}, 'element':{}}
                    # set mode
                    if 'stress' in line:
                        mode = 'stress'
                    infile.readline()
                    line = infile.readline().strip()

                # reset the read type if we hit a blank line
                if line == '':
                    mode = 'off'

                if mode == 'stress':
                    # store stress results
                    tmp = line.split()
                    labs = ['Sx', 'Sy', 'Sz', 'Sxy', 'Sxz', 'Syz', 'Seqv']
                    vals = [float(a) for a in tmp[2:]]
                    seqv = self.__seqv(vals)
                    [s_1, s_2, s_3] = self.__principals(vals)
                    vals.append(seqv)
                    vals += [s_1, s_2, s_3]
                    #print(vals)
                    enum = int(tmp[0])
                    ipnum = int(tmp[1]) # integration point number
                    adict = {}
                    for (label, val) in zip(labs, vals):
                        adict[label] = val
                    if enum not in self.__results[time]['element']:
                        start_val = {'ipoints': {},
                                     'avg': {},
                                     'min': {},
                                     'max': {}}
                        self.__results[time]['element'][enum] = start_val
                    # each line is an integration point result
                    self.__results[time]['element'][enum]['ipoints'][ipnum] = adict

            infile.close()

            # loop over all element results, calculating avg element result
            # by averaging integration point vals
            for time in self.__steps:
                for edict in self.__results[time]['element'].values():
                    ipoints = edict['ipoints'].values()
                    labs = ['Sx', 'Sy', 'Sz', 'Sxy', 'Sxz', 'Syz']
                    field_dict = {}
                    # set field values in max, min, avg locations
                    for field in labs:
                        field_vals = [ipt[field] for ipt in ipoints]
                        field_avg = sum(field_vals)/len(field_vals)
                        field_max, field_min = max(field_vals), min(field_vals)
                        edict['avg'][field] = field_avg
                        edict['max'][field] = field_max
                        edict['min'][field] = field_min
                        field_dict[field] = field_vals
                    # for each element, caclulate Seqv, S1, S2, S3 at average
                    locs = ['avg', 'max', 'min']
                    for loc in locs:
                        seqv, s_1, s_2, s_3 = None, None, None, None
                        if loc == 'avg':
                            vals = [edict[loc][lab] for lab in labs]
                            seqv = self.__seqv(vals)
                            [s_1, s_2, s_3] = self.__principals(vals)
                        elif mode == 'max':
                            seqv = max(field_dict['Seqv'])
                            s_1 = max(field_dict['S1'])
                            s_2 = max(field_dict['S2'])
                            s_3 = max(field_dict['S3'])
                        elif mode == 'min':
                            seqv = min(field_dict['Seqv'])
                            s_1 = min(field_dict['S1'])
                            s_2 = min(field_dict['S2'])
                            s_3 = min(field_dict['S3'])
                        edict[loc]['Seqv'] = seqv
                        edict[loc]['S1'] = s_1
                        edict[loc]['S2'] = s_2
                        edict[loc]['S3'] = s_3

            print('The following times have been read:')
            print(self.__steps)
            print('Results from file: %s have been read.' % fname)

        else:
            print('Error: %s file not found' % fname)
