# Tk-Plots of scattering functions


from Tkinter import *

from Scientific import N
from Scientific.IO.NetCDF import NetCDFFile
from Scientific.TkWidgets.TkPlotCanvas import PlotCanvas, PolyLine, VerticalLine, PlotGraphics
from Scientific.TkWidgets import FloatEntry, IntEntry, StatusBar, ModalDialog

class ScatteringFunctionPlot(Frame):

    def __init__(self, master, filename):
        Frame.__init__(self, master)

        file = NetCDFFile(filename)
        self.q = file.variables['q'][1:]
        self.t = file.variables['time'][:]
        self.fn = file.variables['sf'][1:, :]

        Label(self, text='S(t) for various q').pack(side=TOP, fill=X)
        self.plot1 = PlotCanvas(self, 600, 250, zoom=1)
        self.plot1.pack(side=TOP, fill=BOTH, expand=YES)
        Label(self, text='S(q) for various t').pack(side=TOP, fill=X)
        self.plot2 = PlotCanvas(self, 600, 250, zoom=1)
        self.plot2.pack(side=TOP, fill=BOTH, expand=YES)
        frame = Frame(self)
        frame.pack(side=TOP, fill=X)
        self.first_q = IntEntry(frame, "q range:  from ", 0, 0, len(self.q))
        self.first_q.grid(row=0, column=0)
        self.first_q.bind('<Return>', self.draw)
        self.last_q = IntEntry(frame, " to ", len(self.q), 0, len(self.q))
        self.last_q.grid(row=0, column=1)
        self.skip_q = IntEntry(frame, " skip ", (len(self.q)+10)/10,
                               1, len(self.q))
        self.skip_q.grid(row=0, column=2)
        self.first_t = IntEntry(frame, "t range:  from ", 0, 0, len(self.t))
        self.first_t.grid(row=1, column=0)
        self.last_t = IntEntry(frame, " to ", len(self.t), 0, len(self.t))
        self.last_t.grid(row=1, column=1)
        self.skip_t = IntEntry(frame, " skip ", (len(self.t)+10)/10,
                               1, len(self.t))
        self.skip_t.grid(row=1, column=2)

        self.first_q.bind('<Return>', self.draw)
        self.last_q.bind('<Return>', self.draw)
        self.skip_q.bind('<Return>', self.draw)
        self.first_t.bind('<Return>', self.draw)
        self.last_t.bind('<Return>', self.draw)
        self.skip_t.bind('<Return>', self.draw)

        self.draw()

    def draw(self, event=None):
        try:
            first_q = self.first_q.get()
            last_q = self.last_q.get()
            skip_q = self.skip_q.get()
            first_t = self.first_t.get()
            last_t = self.last_t.get()
            skip_t = self.skip_t.get()
        except ValueError:
            return
        graphics1, graphics2 = self.plot_objects(first_q, last_q, skip_q,
                                                 first_t, last_t, skip_t)
        self.plot1.clear()
        self.plot1.draw(PlotGraphics(graphics1), 'automatic', 'automatic')
        self.plot2.clear()
        self.plot2.draw(PlotGraphics(graphics2), 'automatic', 'automatic')

    def plot_objects(self, first_q, last_q, skip_q,
                     first_t, last_t, skip_t):
        plot_objects1 = []
        markers1 = []
        color_index = 0
        dt = (last_t-first_t+150)/150
        for i in range(first_q, last_q, skip_q):
            data = self.fn[i, first_t:last_t:dt]
            data_t = self.t[first_t:last_t:dt]
            points = N.transpose(N.array([data_t, data]))
            plot_objects1.append(PolyLine(points,
                                          color=self.colors[color_index]))
            markers1.append(self.q[i], color_index)
            color_index = (color_index + 1) % len(self.colors)

        plot_objects2 = []
        markers2 = []
        color_index = 0
        dq = (last_q-first_q+150)/150
        for i in range(first_t, last_t, skip_t):
            data = self.fn[first_q:last_q:dq, i]
            data_q = self.q[first_q:last_q:dq]
            points = N.transpose(N.array([data_q, data]))
            plot_objects2.append(PolyLine(points,
                                          color=self.colors[color_index]))
            markers2.append(self.t[i], color_index)
            color_index = (color_index + 1) % len(self.colors)

        for q, color_index in markers1:
            plot_objects2.append(VerticalLine(q,
                                              color=self.colors[color_index],
                                              stipple='gray50'))
        for t, color_index in markers2:
            plot_objects1.append(VerticalLine(t,
                                              color=self.colors[color_index],
                                              stipple='gray50'))
        return plot_objects1, plot_objects2

    colors = ['red', 'green', 'blue', 'orange', 'brown',
              'yellow', 'violet', 'pink', 'cyan', 'grey']



if __name__ == '__main__':

    import sys
    root = Tk()
    root.title(sys.argv[1])
    plot = ScatteringFunctionPlot(root, sys.argv[1])
    plot.pack(side=TOP, fill=BOTH, expand=YES)
    root.mainloop()
