import numpy as np #numerical python module
import pylab as plt #plotting module
import matplotlib as mplt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import math

################################################################################
#    Class from AJ Medford (?)          #
################################################################################

from scipy.interpolate import InterpolatedUnivariateSpline as spline #used for barrier curves. Don't worry about it
class Plotter:
    def __init__(self,energies,barriers=[],labels=[]):
        self.energies = energies
        self.barriers = barriers
        self.labels = labels
        self.energy_line_args = {'color':'k','lw':2}
        self.barrier_line_args = {'color':'k','lw':2}
        self.label_args = {'color':'k','size':8,'rotation':0}
        self.label_positions= None
        self.initial_energy = 0
        self.initial_stepnumber = 0
        self.energy_mode ='relative' #absolute
        self.energy_line_widths = 0.5

    def draw(self,ax=None):

        def attr_to_list(attrname,required_length=len(self.energies)):
            try: 
                getattr(self,attrname)[0] #Ensure that it is a list
                iter(getattr(self,attrname)) #Ensure that it is a list...
                if len(getattr(self,attrname)) == required_length:
                    pass
                else:
                    raise ValueError(attrname + ' list is of length '+ \
                            str(len(getattr(self,attrname)))+ \
                            ', but needs to be of length ' + \
                            str(required_length))
                return getattr(self,attrname)
            except:
                return [getattr(self,attrname)]*required_length
        
        barrier_line_args = attr_to_list('barrier_line_args',
                len(self.energies)-1)
        energy_line_widths = attr_to_list('energy_line_widths')
        energy_line_args = attr_to_list('energy_line_args')
        label_args =attr_to_list('label_args')
        label_positions=attr_to_list('label_positions')

        #plot energy lines
        energy_list = np.array(self.energies)
        energy_list = (energy_list - energy_list[0])
        energy_list = list(energy_list)
        if self.energy_mode == 'relative':
            cum_energy = [energy_list[0]]
            for i,e in enumerate(energy_list[1:]):
                last = cum_energy[i]+e
                cum_energy.append(last)
            energy_list = cum_energy
        energy_list = np.array(energy_list) + self.initial_energy
        energy_list = list(energy_list)
        energy_lines = [
                [[i+self.initial_stepnumber,i+width+self.initial_stepnumber],
                    [energy_list[i]]*2] 
                for i,width in enumerate(energy_line_widths)]
        for i,line in enumerate(energy_lines):
            if i!=0:
                energy_line_args[i]['label']=None
            ax.plot(*line,**energy_line_args[i])

        #create barrier lines
        barrier_lines = []
        if not self.barriers: self.barriers = [0]*(len(self.energies)-1)
        for i,barrier in enumerate(self.barriers):
            xi = energy_lines[i][0][1]
            xf = energy_lines[i+1][0][0]
            yi = energy_lines[i][1][0]
            yf = energy_lines[i+1][1][0]
            if barrier == 0 or barrier <= yf-yi:
                line = [[xi,xf],[yi,yf]]
            else:
                yts = yi+barrier
                barrier_rev = barrier + (yi-yf)
                ratio = np.sqrt(barrier)/(np.sqrt(barrier)+np.sqrt(barrier_rev))
                xts = xi + ratio*(xf-xi)
                xs = [xi,xts,xf]
                ys = [yi,yts,yf]
                f = spline(xs,ys,k=2)
                newxs = np.linspace(xi,xf,20)
                newys = f(newxs)
                line = [newxs,newys]
            barrier_lines.append(line)
        #plot barrier lines
        for i,line in enumerate(barrier_lines):
            ax.plot(*line,**barrier_line_args[i])
        #add labels
        trans = ax.get_xaxis_transform()
        for i,label in enumerate(self.labels):
            xpos = sum(energy_lines[i][0])/len(energy_lines[i][0])-0.1
            label_position = label_positions[i]
            args = label_args[i]
            if label_position in ['top','ymax']:
                ypos = 1
                args['transform'] = trans
                ax.text(xpos,ypos,label,**args)
            elif label_position in ['bot','bottom','ymin']:
                ypos = -0.1 
                ax.xaxis.set_ticks([float(sum(line[0])/len(line[0])) 
                    for line in energy_lines])
                ax.set_xticklabels(self.labels)
                for attr in args.keys():
                    try:
                        [getattr(t,'set_'+attr)(args[attr]) 
                                for t in ax.xaxis.get_ticklabels()]
                    except:
                        pass
            else:
                ypos = energy_lines[i][1][0]
                ax.annotate(label,[xpos,ypos],**args)

################################################################################
#   End of energy diagram class. Actual problem starts below.                  #
################################################################################


