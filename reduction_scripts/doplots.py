def some_plots(xx,yy,dolegend=False,dohalpha=False,donh3=False):
    c13 = h213cocube.get_spectrum(xx,yy)
    c12 = h2cocube.get_spectrum(xx,yy)
    c12.baseline(exclude=[-225,200],order=5)
    c13.baseline(exclude=[-225,200],order=5)
    c13.plotter(label='H$_{2}$$^{13}$CO',axis=gca())
    c12.plotter(axis=c13.plotter.axis,clear=False,color='b',label="H$_{2}$CO")
    (c13*6).plotter(label='6$\\times$H$_{2}$$^{13}$CO',axis=gca(),color='r',clear=False)
    if dolegend:
        legend(loc='best')
    if dohalpha:
        halpha = h110acube.get_spectrum(xx,yy)
        halpha.baseline(exclude=[-225,200],order=5)
        (halpha*5).plotter(axis=c13.plotter.axis,clear=False,color='darkgreen')
    if donh3:
        nh3 = nh3cube.get_spectrum(xx,yy)
        nh3.baseline(exclude=[-225,200],order=5)
        (nh3*5).plotter(axis=c13.plotter.axis,clear=False,color='orange')
    gca().set_xlim(-100,100)
    draw()

