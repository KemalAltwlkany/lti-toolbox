import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

from scipy import signal

s22 = math.sqrt(2) / 2  # used a lot

class NyquistPlot:
    def __init__(self, tf=None, xmin=-5, xmax=5, ymin=-5, ymax=5, ticks_frequency=1, label_plane=True,
                 planelabel=r'$G\{j\omega\}$', label_axis=True, xlabel=r'$u(\omega)$', ylabel=r'$jv(\omega)$', ticks=True, grid=True,
                 identical_scales=True):
        self.tf = tf
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.ticks_frequency = ticks_frequency
        self.label_plane = label_plane
        self.planelabel = planelabel
        self.label_axis = label_axis
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.ticks = ticks
        self.grid = grid
        self.identical_scales = identical_scales
        self.zero_values, self.pole_values, self.zero_orders, self.pole_orders = None, None, None, None

        self.fig, self.ax = None, None

    # draws the coordinate system, labels the axis and plane.
    def draw_template(self):
        fig, ax = plt.subplots(figsize=(11, 11))

        # identical scales for both axes
        if self.identical_scales is True:
            ax.set(xlim=(self.xmin - 1.5, self.xmax + 1.5), ylim=(self.ymin - 1.5, self.ymax + 1.5), aspect='equal')

        # set bottom and left spines as x and y axes of coordinate system
        ax.spines['bottom'].set_position('zero')
        ax.spines['left'].set_position('zero')

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Create 'sigma' and 'j omega' labels placed at the end of the axes
        if self.label_axis:
            ax.set_xlabel(self.xlabel, size=11, labelpad=-20, y=0.90, x=1.03)
            ax.set_ylabel(self.ylabel, size=11, labelpad=-18, y=1.00, rotation=0)

        # Create '{s}' plane label, located in upper right corner
        if self.label_plane:
            ax.annotate(self.planelabel, xy=(self.xmax, self.ymax), xycoords='data', fontsize=11)

        # Create custom major ticks to determine position of tick labels
        if self.ticks is True:
            x_ticks = np.arange(self.xmin, self.xmax + 1, self.ticks_frequency)
            y_ticks = np.arange(self.ymin, self.ymax + 1, self.ticks_frequency)
            ax.set_xticks(x_ticks[x_ticks != 0])
            ax.set_yticks(y_ticks[y_ticks != 0])
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        # Create minor ticks placed at each integer to enable drawing of minor grid
        # lines: note that this has no effect in this example with ticks_frequency=1
        ax.set_xticks(np.arange(self.xmin, self.xmax + 1), minor=True)
        ax.set_yticks(np.arange(self.ymin, self.ymax + 1), minor=True)

        # Draw major and minor grid lines
        if self.grid is True:
            ax.grid(which='both', color='grey', linewidth=1, linestyle='-', alpha=0.2)

        # Draw arrows
        arrow_fmt = dict(markersize=4, color='black', clip_on=False)
        ax.plot((1), (0), marker='>', transform=ax.get_yaxis_transform(), **arrow_fmt)
        ax.plot((0), (1), marker='^', transform=ax.get_xaxis_transform(), **arrow_fmt)

        self.fig, self.ax = fig, ax
        self.ax.set_axisbelow(True)  # puts grid "behind" all other figures.
        # https://stackoverflow.com/questions/1726391/matplotlib-draw-grid-lines-behind-other-graph-elements

        return self.fig, self.ax

    def astatism(self, start_theta, end_theta, ls='-'):
        """
        Draws the astatism part of the Nyquist plot. Arc is drawn in clock-wise direction, so start and end theta
            are reversed.
        :param start_pt: starting phase
        :param end_pt: ending phase
        :return:
        """
        patch = patches.Arc(xy=(0, 0), width=2 * self.ymax, height=2 * self.ymax, angle=0, theta1=start_theta, theta2=end_theta,
                            lw=2.5, ls=ls)
        self.ax.add_patch(patch)

    def inf_R_arrow(self, quadrant=4):
        # R -> inf arrow
        xy = (0, 0)
        dx = 0
        dy = 0
        if quadrant == 1:
            xy = (3 * self.ymax / 4, 3 * self.ymax / 7)
            dx = self.ymax * math.sqrt(3) / 2
            dy = self.ymax / 2
        elif quadrant == 2:
            xy=(-3 * self.ymax / 4, 3 * self.ymax / 7)
            dx = -self.ymax * math.sqrt(3) / 2
            dy = self.ymax / 2
        elif quadrant == 3:
            xy=(-3 * self.ymax / 4, -3 * self.ymax / 7)
            dx = -self.ymax * math.sqrt(3) / 2
            dy = -self.ymax / 2
        elif quadrant == 4:
            xy=(3 * self.ymax / 4, -3 * self.ymax / 7)
            dx = self.ymax * math.sqrt(3) / 2
            dy = -self.ymax / 2
        plt.arrow(x=0, y=0, dx=dx, dy=dy,
                  shape='full', ls='--', lw=0.5, length_includes_head=True, head_width=.35, color='k')
        self.ax.annotate(r'$R\rightarrow \infty$', xy=xy, xycoords='data',
                         fontsize=12)

    def arc_direction_arrow(self, quadrant=4, reversed=False):
        x, y, dx, dy = 0, 0, 0, 0
        if quadrant == 1:
            x = self.ymax * s22
            y = self.ymax * s22
            dx = 0.01
            dy = -0.01
        elif quadrant == 2:
            x = -self.ymax * s22
            y = self.ymax * s22
            dx = 0.01
            dy = 0.01
        elif quadrant == 3:
            x = -self.ymax * s22
            y = -self.ymax * s22
            dx = -0.01
            dy = 0.01
        elif quadrant == 4:
            x=self.ymax * s22
            y=-self.ymax * s22
            dx = -0.01
            dy = -0.01
        if reversed is True:
            dx = -1 * dx
            dy = -1 * dy
        plt.arrow(x=x, y=y, dx=dx, dy=dy,
                  shape='full', lw=4.5, length_includes_head=True, head_width=.15, color='k', zorder=25)


def t05_z01():


    # PLOT 1, Astatism
    # ----
    # Astatism
    # o1 = NyquistPlot(xmin=-12, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    # o1.draw_template()
    # o1.astatism(-180, 0)
    # o1.arc_direction_arrow(quadrant=3)
    # o1.arc_direction_arrow(quadrant=4)
    # o1.inf_R_arrow(4)
    #
    # o1.ax.annotate(r'$A^\prime$', (o1.xmax, 0.1), fontsize=14)
    # o1.ax.annotate(r'$B^\prime$', (-o1.xmax, 0.1), fontsize=14)
    # plt.savefig(os.getcwd() + '/T05/figs/z01/pt1.png', bbox_inches='tight')
    # plt.show()
    # ----

    # PLOT 2, u (w) >> v(w)
    # ----
    # o2 = NyquistPlot(xmin=-12, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    # o2.draw_template()
    # o2.ax.annotate(r'$u(\omega) \gg v(\omega)$', (-o2.xmax, 3), fontsize=14, color='red')
    # plt.plot([-o2.xmax], [-1], marker='o', color='red', markersize=12)
    # plt.arrow(x=-10, y=2.5, dx=-1.5, dy=-3,
    #           shape='full', ls='--', lw=0.5, length_includes_head=True, head_width=.35, color='r')
    # plt.savefig(os.getcwd() + '/T05/figs/z01/pt2.png', bbox_inches='tight')
    # plt.show()
    # ----


    # PLOT 3, jw beginning, u(w) >> v(w)
    # ----
    # o3 = NyquistPlot(xmin=-12, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    # o3.draw_template()
    # o3.astatism(-180, 0)
    # o3.arc_direction_arrow(quadrant=3)
    # o3.arc_direction_arrow(quadrant=4)
    # o3.inf_R_arrow(4)
    #
    # o3.ax.annotate(r'$A^\prime$', (o3.xmax, 0.1), fontsize=14)
    # o3.ax.annotate(r'$B^\prime$', (-o3.xmax, 0.1), fontsize=14)
    # plt.savefig(os.getcwd() + '/T05/figs/z01/pt3.png', bbox_inches='tight')
    # plt.plot([-o3.xmax, -o3.xmax + 2], [-0.07, -0.07], lw='2.5', color='k')
    # patch = patches.Arc(xy=(-9, 0), width=2, height=4, angle=0, theta1=-180,
    #                     theta2=-90,
    #                     lw=2.5)
    # o3.ax.add_patch(patch)
    # plt.plot([-o3.xmax+3, -o3.xmax + 4.5], [-2, -2], lw='2.5', color='k')
    # o3.ax.annotate(r'$B^\prime$', (-o3.xmax+2, 0.1), fontsize=14)
    # o3.ax.annotate(r'$u(\omega) \gg v(\omega)$', (-9, -5), fontsize=14, color='k')
    # plt.arrow(x=-8.7, y=-4.5, dx=0.5, dy=2.5,
    #            shape='full', ls='--', lw=0.5, length_includes_head=True, head_width=.35, color='k')
    # plt.savefig(os.getcwd() + '/T05/figs/z01/pt3.png', bbox_inches='tight')
    # plt.show()
    # ----


    # PLOT 4, final diagram for positive omega
    # ----
    # o4 = NyquistPlot(xmin=-12, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    # o4.draw_template()
    # o4.astatism(-180, 0)
    # o4.arc_direction_arrow(quadrant=3)
    # o4.arc_direction_arrow(quadrant=4)
    # o4.inf_R_arrow(4)
    #
    # o4.ax.annotate(r'$A^\prime$', (o4.xmax, 0.1), fontsize=14)
    # o4.ax.annotate(r'$B^\prime$', (-o4.xmax, 0.1), fontsize=14)
    # plt.savefig(os.getcwd() + '/T05/figs/z01/pt3.png', bbox_inches='tight')
    # plt.plot([-o4.xmax, -o4.xmax + 2], [-0.07, -0.07], lw='2.5', color='k')
    # patch = patches.Arc(xy=(-9, 0), width=2, height=4, angle=0, theta1=-180,
    #                     theta2=-90,
    #                     lw=2.5)
    # o4.ax.add_patch(patch)
    # #plt.plot([-o4.xmax + 3, -o4.xmax + 4.5], [-2, -2], lw='2.5', color='k')
    # o4.ax.annotate(r'$B^\prime$', (-o4.xmax + 2, 0.1), fontsize=14)
    #
    # # Bezier curve
    # verts = [
    #     (-9, -2), # start point P0
    #     (-2.9, -2.4), # control point P1
    #     (-2.2, -1.2), # end point P2
    # ]
    # codes = [
    #     Path.MOVETO,
    #     Path.CURVE3,
    #     Path.CURVE3
    # ]
    # path = Path(verts, codes)
    # patch = patches.PathPatch(path, facecolor='none', lw=2.5)
    # o4.ax.add_patch(patch)
    #
    # verts = [
    #     (-2.2, -1.2),
    #     (-1.9, 0.1),
    #     (-1.1, -0.05),
    # ]
    # path = Path(verts, codes)
    # patch = patches.PathPatch(path, facecolor='none', lw=2.5)
    # o4.ax.add_patch(patch)
    #
    # plt.plot([-1.1, 0], [-0.025, -0.025], color='k', lw=2.5)
    #
    #
    # plt.arrow(x=-2.31, y=-1.24, dx=0.01, dy=0.015,
    #           shape='full', lw=4, length_includes_head=True, head_width=.15, color='k', zorder=25)
    #
    #
    # o4.ax.annotate(r'$\phi_{ul} = -180\degree$', (1.5, 1.5), fontsize=10)
    # plt.arrow(x=1.4, y=1.4, dx=-1.3, dy=-1.3,
    #           shape='full', ls='--', lw=0.5, length_includes_head=True, head_width=.35, color='k', zorder=25)
    # plt.savefig(os.getcwd() + '/T05/figs/z01/pt4.png', bbox_inches='tight')
    # plt.show()
    # ----
    pass


# def example1():
#     tf = signal.lti([10], [1, 4, 0])
#     w, H = signal.freqresp(tf)
#     plt.plot(H.real, H.imag, 'b')
#     plt.plot(H.real, -H.imag, 'r')
#     plt.show()


if __name__ == '__main__':
    # t05_z01()
    pass

