from pylab import *
import matplotlib as mpl 
from matplotlib.font_manager import fontManager, FontProperties
import pylab
from pylab import arange,pi,sin,cos,sqrt

def plot(x, y):    
    markers =     ['-k', '-b',  '-c', '-y',  '-r', '-m', '-b', '-.^g', '-vm', '--Dr', '--Dy']
    markersface = [ 'w',  'w',   'w',  'w',   'w',  'w',  'w', 'g',   'w',    'w',   'y']
    markersedge = [ 'w',  'w',   'w',  'w',   'w',  'w',  'w', 'g',   'm',    'r',   'y']

    scale = 2.0

    # Linear regression of points in 
    # logarithmic scale using x and y vector.
    xp  = log(x)
    yp  = log(y);
    m, b = polyfit(xp, yp, 1)


    x_half_order = []
    y_half_order = []
    vx = scale * x[0]
    vy = sqrt(scale) * y[0]
    for i in xrange(len(x)+2):
        x_half_order.append(vx)
        y_half_order.append(vy)
        vx /= scale
        vy /= sqrt(scale)


    x_first_order = []
    y_first_order = []
    vx = scale * x[0]
    vy = scale * y[0]
    for i in xrange(len(x)+2):
        x_first_order.append(vx)
        y_first_order.append(vy)
        vx /= scale
        vy /= scale


    x_second_order = []
    y_second_order = []
    vx = scale * x[0]
    vy = scale * scale * y[0]
    for i in xrange(len(x)+2):
        x_second_order.append(vx)
        y_second_order.append(vy)
        vx /= scale
        vy /= scale * scale


    x_third_order = []
    y_third_order = []
    vx = scale * x[0]
    vy = scale * scale * scale * y[0]
    for i in xrange(len(x)+2):
        x_third_order.append(vx)
        y_third_order.append(vy)
        vx /= scale
        vy /= scale * scale * scale


    x_fourth_order = []
    y_fourth_order = []
    vx = scale * x[0]
    vy = scale * scale * scale * scale * y[0]
    for i in xrange(len(x)+2):
        x_fourth_order.append(vx)
        y_fourth_order.append(vy)
        vx /= scale
        vy /= scale * scale * scale * scale

    x_fifth_order = []
    y_fifth_order = []
    vx = scale * x[0]
    vy = scale * scale * scale * scale * scale * y[0]
    for i in xrange(len(x)+2):
        x_fifth_order.append(vx)
        y_fifth_order.append(vy)
        vx /= scale
        vy /= scale * scale * scale * scale * scale

    x_sixth_order = []
    y_sixth_order = []
    vx = scale * x[0]
    vy = scale * scale * scale * scale * scale * scale * y[0]
    for i in xrange(len(x)+2):
        x_sixth_order.append(vx)
        y_sixth_order.append(vy)
        vx /= scale
        vy /= scale * scale * scale * scale * scale * scale

    #First order convergence
    methods = [['1/2th order slope', y_half_order, x_half_order]]
    methods += [['1st order slope', y_first_order, x_first_order]]
    methods += [['2nd order slope', y_second_order, x_second_order]]
    methods += [['3rd order slope', y_third_order, x_third_order]]
    methods += [['4th order slope', y_fourth_order, x_fourth_order]]
    methods += [['5th order slope', y_fifth_order, x_fifth_order]]
    methods += [['6th order slope', y_sixth_order, x_sixth_order]]
    
    #Double precision:
    methods += [['Implementation ($p = %2.2f$)' % m, y, x]]


    # Plotting
    fig = gcf()
    fig.set_facecolor('white')

    legends = []
    k = 0
    for i in methods:
        plot = loglog(i[2], i[1], markers[k], linewidth=1.5,           \
            markerfacecolor=markersface[k], markeredgecolor=markersedge[k], \
               markersize=8, antialiased = True, basex=10)    
        k += 1
        legends += [i[0]]    

    axes = plot[0].get_axes()

    mpl.rcParams['grid.color'] = '1.0'

    # Setting graph parameters
    params = {'xtick.labelsize' : 20,\
                  'ytick.labelsize' : 20,\
                  'grid.color' : (0.5, 0.5, 0.5, 1.0)
              }
    pylab.rcParams.update(params)
    #pylab.xlim(4.5 * 10**(-3), 1.0)

    # Plotting settings
    grid(True)
    xlabel(r"$\log(N)$", size=25)
    ylabel('$\log(L_\infty$ norm) ', size=25)

    font= FontProperties(size=15);
    legend(legends, loc='lower right', prop=font)

    show()
