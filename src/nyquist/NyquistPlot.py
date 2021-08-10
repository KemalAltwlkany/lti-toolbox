import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


s22 = math.sqrt(2) / 2  # used a lot


class NyquistPlot:
    def __init__(self, tf=None, xmin=-5, xmax=5, ymin=-5, ymax=5, ticks_frequency=1, label_plane=True,
                 planelabel=r'$\{s\}$', label_axis=True, xlabel=r'$\sigma$', ylabel=r'$j\omega$', ticks=True, grid=True,
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
        self.compute_zeros_poles()

        self.fig, self.ax = None, None

    # TODO:
    # If future upgrades allow entering the tf in a zpk format, this step should obviously be skipped.
    def compute_zeros_poles(self):
        # assumes tf is a list of lists, 1st are the numerator coeffs, 2nd the denom coeffs
        self.zero_values = np.roots(self.tf[0])
        self.zero_values = np.around(self.zero_values, decimals=2)
        # round zeros to 2 decimals. "around" used for complex numbers.
        # https://numpy.org/doc/stable/reference/generated/numpy.around.html
        self.zero_values, self.zero_orders = np.unique(self.zero_values, return_counts=True)

        self.pole_values = np.roots(self.tf[1])
        self.pole_values = np.around(self.pole_values, decimals=2)
        self.pole_values, self.pole_orders = np.unique(self.pole_values, return_counts=True)

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
            ax.set_xlabel(self.xlabel, size=11, labelpad=-20, x=1.00)
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

    def draw_roots(self):
        print(self.zero_values, self.zero_orders)
        print(self.pole_values, self.pole_orders)

        for z_val, z_ord in zip(self.zero_values, self.zero_orders):
            self.draw_root(rtype='zero', coords=z_val, multiplicity=z_ord)

        for p_val, p_ord in zip(self.pole_values, self.pole_orders):
            self.draw_root(rtype='pole', coords=p_val, multiplicity=p_ord)

    # Marks a zero or a pole on the Nyquist plot, also supports multiplicity
    def draw_root(self, rtype='zero', coords=(0, 0), multiplicity=0):
        # if zero/pole {s} is real, the y coordinate is 0, otherwise (x, y) = Re{s}, Im{s}
        if np.iscomplex(coords):
            coords = (np.real(coords), np.imag(coords))
        else:
            coords = (coords, 0)  # Im{s} = 0

        # determine marker type (zero - 'o', pole - 'x')
        marker = 'o'
        if rtype != 'zero':
            marker = 'x'

        # plot the zero/pole
        plt.plot(coords[0], coords[1], lw=2, marker=marker, markersize=8, markeredgecolor='k', markerfacecolor='k',
                 fillstyle='full', markeredgewidth=2.25)

        # plot circles if zero/pole is of higher order (has a multiplicity > 1)
        i = 0
        while multiplicity > 1:
            multiplicity = multiplicity - 1
            patch = patches.Arc(xy=(coords[0], coords[1]), width=0.5 + i * 0.4, height=0.5 + i * 0.4, angle=0, theta1=0,
                                theta2=360, lw=2.25)
            self.ax.add_patch(patch)
            i = i + 1

    # TODO:
    # this function does not work perfectly. needs to be fixed.
    # Assumes system does NOT have more than one pair of complex-conjugate poles on the imaginary axis.
    # Problems with contour size may happen as well: e.g., a pole in origin and a zero/pole at s=0.5 or s=1.
    def draw_contour(self):
        # pole in origin
        # 1st part of contour (potentially)
        origin_pole = False
        for pole in self.pole_values:
            if pole == 0:
                origin_pole = True
                patch = patches.Arc(xy=(0, 0), width=2.0, height=2.0, angle=0, theta1=-90, theta2=90, lw=2.25)
                self.ax.add_patch(patch)
                plt.arrow(x=s22, y=s22, dx=-0.01, dy=0.01,
                          shape='full', lw=4.5, length_includes_head=True, head_width=.15, color='k', zorder=25)
                break

        # (only for systems with complex-conjugate poles on imaginary axis)
        imag_value = 0
        for pole in self.pole_values:
            if np.iscomplex(pole) and np.real(pole) == 0:
                imag_value = np.abs(np.imag(pole))
                patch = patches.Arc(xy=(np.real(pole), np.imag(pole)), width=2.0, height=2.0, angle=0, theta1=-90, theta2=90, lw=2.25)
                self.ax.add_patch(patch)
                plt.arrow(x=np.real(pole)+s22, y=np.imag(pole)+s22, dx=-0.01, dy=0.01,
                          shape='full', lw=4.5, length_includes_head=True, head_width=.15, color='k', zorder=25)

                # repeat above for complex-conjugate pole
                patch = patches.Arc(xy=(np.real(pole), -np.imag(pole)), width=2.0, height=2.0, angle=0, theta1=-90,
                                    theta2=90, lw=2.25)
                self.ax.add_patch(patch)
                plt.arrow(x=np.real(pole) + s22, y=-np.imag(pole) + s22, dx=-0.01, dy=0.01,
                          shape='full', lw=4.5, length_includes_head=True, head_width=.15, color='k', zorder=25)
                break

        # 2nd part of contour
        omega_axis_start = 0
        if origin_pole is True:
            omega_axis_start = 1
        print(imag_value)
        plt.plot([0, 0], [omega_axis_start, omega_axis_start+imag_value-1], color='k', lw=2.25)
        plt.plot([0, 0], [omega_axis_start + imag_value + 1, self.ymax], color='k', lw=2.25)
        plt.arrow(x=0, y=self.ymax/2, dx=0.00, dy=0.01,
                  shape='full', lw=4.5, length_includes_head=True, head_width=.15, color='k', zorder=25)

        plt.plot([0, 0], [-omega_axis_start, -omega_axis_start-imag_value+1], color='k', lw=2.25)
        plt.plot([0, 0], [-omega_axis_start - imag_value-1, -self.ymax], color='k', lw=2.25)
        plt.arrow(x=0, y=-self.ymax/2, dx=0.00, dy=0.01,
                  shape='full', lw=4.5, length_includes_head=True, head_width=.15, color='k', zorder=25)


        # 3rd part of the contour
        patch = patches.Arc(xy=(0, 0), width=2*self.ymax, height=2*self.ymax, angle=0, theta1=-90, theta2=90, lw=2.5)
        self.ax.add_patch(patch)
        plt.arrow(x=self.ymax * s22, y=self.ymax * s22, dx=0.01, dy=-0.01,
                  shape='full', lw=4.5, length_includes_head=True, head_width=.15, color='k', zorder=25)
        plt.arrow(x=self.ymax * s22, y=-self.ymax * s22, dx=-0.01, dy=-0.01,
                  shape='full', lw=4.5, length_includes_head=True, head_width=.15, color='k', zorder=25)

        # R -> inf arrow
        plt.arrow(x=omega_axis_start, y=omega_axis_start, dx=self.ymax*math.sqrt(3)/2, dy=-self.ymax/2,
                  shape='full', ls='--', lw=0.5, length_includes_head=True, head_width=.35, color='k')
        self.ax.annotate(r'$R\rightarrow \infty$', xy=(3*self.ymax/4, -3*self.ymax/7), xycoords='data', fontsize=12)


def example1():
    tf = [[1, 3], [1, 3, 2]]
    nyq = NyquistPlot(tf, xmin=-5, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    nyq.draw_template()
    nyq.draw_contour()
    nyq.draw_roots()
    plt.show()


def example2():
    tf = [[1, -3], [1, 1, 0]]
    nyq = NyquistPlot(tf, xmin=-5, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    nyq.draw_template()
    nyq.draw_contour()
    nyq.draw_roots()
    plt.show()


def example3():
    tf = [[1, 4], [1, 0, 4]]
    nyq = NyquistPlot(tf, xmin=-5, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    nyq.draw_template()
    nyq.draw_contour()
    nyq.draw_roots()
    plt.show()


if __name__ == '__main__':
    example1()
    example2()
    example3()
    pass

