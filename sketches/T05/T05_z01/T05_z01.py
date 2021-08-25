import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from scipy import signal

from sketches.templates.nyquistplot import NyquistPlot
from sketches.templates.nyquistcontour import NyquistContour


def pt01():
    # PLOT 1, Astatism
    # ----
    o1 = NyquistPlot(xmin=-12, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    o1.draw_template()
    o1.astatism(-180, 0)
    o1.arc_direction_arrow(quadrant=3)
    o1.arc_direction_arrow(quadrant=4)
    o1.inf_R_arrow(4)

    o1.ax.annotate(r'$A^\prime$', (o1.xmax, 0.1), fontsize=14)
    o1.ax.annotate(r'$B^\prime$', (-o1.xmax, 0.1), fontsize=14)
    plt.savefig(os.getcwd() + '/pt1.png', bbox_inches='tight')
    plt.show()
    # ----


def pt02():
    # PLOT 2, u (w) >> v(w)
    # ----
    o2 = NyquistPlot(xmin=-12, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    o2.draw_template()
    o2.ax.annotate(r'$u(\omega) \gg v(\omega)$', (-o2.xmax, 3), fontsize=14, color='red')
    plt.plot([-o2.xmax], [-1], marker='o', color='red', markersize=12)
    plt.arrow(x=-10, y=2.5, dx=-1.5, dy=-3,
              shape='full', ls='--', lw=0.5, length_includes_head=True, head_width=.35, color='r')
    plt.savefig(os.getcwd() + '/pt2.png', bbox_inches='tight')
    plt.show()
    # ----


def pt03():
    # PLOT 3, jw beginning, u(w) >> v(w)
    # ----
    o3 = NyquistPlot(xmin=-12, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    o3.draw_template()
    o3.astatism(-180, 0)
    o3.arc_direction_arrow(quadrant=3)
    o3.arc_direction_arrow(quadrant=4)
    o3.inf_R_arrow(4)

    o3.ax.annotate(r'$A^\prime$', (o3.xmax, 0.1), fontsize=14)
    o3.ax.annotate(r'$B^\prime$', (-o3.xmax, 0.1), fontsize=14)
    #plt.savefig(os.getcwd() + '/T05/figs/z01/pt3.png', bbox_inches='tight')
    plt.plot([-o3.xmax, -o3.xmax + 2], [-0.07, -0.07], lw='2.5', color='k')
    patch = patches.Arc(xy=(-9, 0), width=2, height=4, angle=0, theta1=-180,
                        theta2=-90,
                        lw=2.5)
    o3.ax.add_patch(patch)
    plt.plot([-o3.xmax+3, -o3.xmax + 4.5], [-2, -2], lw='2.5', color='k')
    o3.ax.annotate(r'$B^\prime$', (-o3.xmax+2, 0.1), fontsize=14)
    o3.ax.annotate(r'$u(\omega) \gg v(\omega)$', (-9, -5), fontsize=14, color='k')
    plt.arrow(x=-8.7, y=-4.5, dx=0.5, dy=2.5,
               shape='full', ls='--', lw=0.5, length_includes_head=True, head_width=.35, color='k')
    plt.savefig(os.getcwd() + '/pt3.png', bbox_inches='tight')
    plt.show()
    # ----


def pt04():
    # PLOT 4, final diagram for positive omega
    # ----
    o4 = NyquistPlot(xmin=-12, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    o4.draw_template()
    o4.astatism(-180, 0)
    o4.arc_direction_arrow(quadrant=3)
    o4.arc_direction_arrow(quadrant=4)
    o4.inf_R_arrow(4)

    o4.ax.annotate(r'$A^\prime$', (o4.xmax, 0.1), fontsize=14)
    o4.ax.annotate(r'$B^\prime$', (-o4.xmax, 0.1), fontsize=14)
    plt.plot([-o4.xmax, -o4.xmax + 2], [-0.07, -0.07], lw='2.5', color='k')
    patch = patches.Arc(xy=(-9, 0), width=2, height=4, angle=0, theta1=-180,
                        theta2=-90,
                        lw=2.5)
    o4.ax.add_patch(patch)
    #plt.plot([-o4.xmax + 3, -o4.xmax + 4.5], [-2, -2], lw='2.5', color='k')
    o4.ax.annotate(r'$B^\prime$', (-o4.xmax + 2, 0.1), fontsize=14)

    # Bezier curve
    verts = [
        (-9, -2), # start point P0
        (-2.9, -2.4), # control point P1
        (-2.2, -1.2), # end point P2
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE3,
        Path.CURVE3
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=2.5)
    o4.ax.add_patch(patch)

    verts = [
        (-2.2, -1.2),
        (-1.9, 0.1),
        (-1.1, -0.05),
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=2.5)
    o4.ax.add_patch(patch)

    plt.plot([-1.1, 0], [-0.025, -0.025], color='k', lw=2.5)


    plt.arrow(x=-2.31, y=-1.24, dx=0.01, dy=0.015,
              shape='full', lw=4, length_includes_head=True, head_width=.15, color='k', zorder=25)


    o4.ax.annotate(r'$\phi_{ul} = -180\degree$', (1.5, 1.5), fontsize=10)
    plt.arrow(x=1.4, y=1.4, dx=-1.3, dy=-1.3,
              shape='full', ls='--', lw=0.5, length_includes_head=True, head_width=.35, color='k', zorder=25)
    plt.savefig(os.getcwd() + '/pt4.png', bbox_inches='tight')
    plt.show()
    # ----


def pt05():
    # PLOT 5, final diagram for BOTH positive and negatives frequencies
    # ----
    # astatism part
    o5 = NyquistPlot(xmin=-12, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    o5.draw_template()

    # positive freq
    o5.astatism(-180, 0)
    o5.arc_direction_arrow(quadrant=3)
    o5.arc_direction_arrow(quadrant=4)

    # negative freq
    o5.astatism(0, 180, ls='--')
    o5.arc_direction_arrow(quadrant=2, reversed=True)
    o5.arc_direction_arrow(quadrant=1, reversed=True)
    o5.inf_R_arrow(4)


    # part jw -> 0+
    # positive freq
    plt.plot([-o5.xmax, -o5.xmax + 2], [-0.07, -0.07], lw='2.5', color='k')
    patch = patches.Arc(xy=(-9, 0), width=2, height=4, angle=0, theta1=-180,
                        theta2=-90,
                        lw=2.5)
    o5.ax.add_patch(patch)

    # negative freq
    plt.plot([-o5.xmax, -o5.xmax + 2], [0.07, 0.07], lw='2.5', color='k', ls='--')
    patch = patches.Arc(xy=(-9, 0), width=2, height=4, angle=0, theta1=90,
                        theta2=180,
                        lw=2.5, ls='--')
    o5.ax.add_patch(patch)

    # Part jw -> +oo
    # Bezier curve, positive freq
    verts = [
        (-9, -2), # start point P0
        (-2.9, -2.4), # control point P1
        (-2.2, -1.2), # end point P2
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE3,
        Path.CURVE3
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=2.5)
    o5.ax.add_patch(patch)

    verts = [
        (-2.2, -1.2),
        (-1.9, 0.1),
        (-1.1, -0.05),
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=2.5)
    o5.ax.add_patch(patch)

    plt.plot([-1.1, 0], [-0.025, -0.025], color='k', lw=2.5)


    plt.arrow(x=-2.31, y=-1.24, dx=0.01, dy=0.015,
              shape='full', lw=4, length_includes_head=True, head_width=.15, color='k', zorder=25)

    # Bezier, negative freq
    verts = [
        (-9, 2), # start point P0
        (-2.9, 2.4), # control point P1
        (-2.2, 1.2), # end point P2
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE3,
        Path.CURVE3
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=2.5, ls='--')
    o5.ax.add_patch(patch)

    verts = [
        (-2.2, 1.2),
        (-1.9, 0.1),
        (-1.1, 0.05),
    ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=2.5, ls='--')
    o5.ax.add_patch(patch)

    plt.plot([-1.1, 0], [0.025, 0.025], color='k', lw=2.5, ls='--')


    plt.arrow(x=-2.31, y=1.24, dx=0.01, dy=-0.015,
              shape='full', lw=4, length_includes_head=True, head_width=.15, color='k', zorder=25)

    o5.ax.annotate(r'$\omega_{-}$', (o5.xmax, 5), fontsize=14)
    o5.ax.annotate(r'$\omega_{+}$', (o5.xmax, -5), fontsize=14)
    plt.savefig(os.getcwd() + '/pt5.png', bbox_inches='tight')
    plt.show()
    # ----


def contour():
    tf = [[3, 3], [1, 3, 0, 0]]
    nyq = NyquistContour(tf, xmin=-5, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    nyq.draw_template()
    nyq.draw_contour()
    nyq.draw_roots()
    nyq.annotate_zero_pole('-1', (-1.5, -1), fontsize=12)
    nyq.annotate_zero_pole('-3', (-3.5, -1), fontsize=12)
    plt.savefig(os.getcwd() + '/contour.png', bbox_inches='tight')
    plt.show()


def pole_zero():
    tf = [[3, 3], [1, 3, 0, 0]]
    nyq = NyquistContour(tf, xmin=-5, xmax=5, ymin=-4, ymax=4, label_axis=True, grid=True, ticks=False, figsize=(5, 5))
    nyq.draw_template()
    nyq.draw_roots()
    nyq.annotate_zero_pole('-1', (-1.5, -1), fontsize=12)
    nyq.annotate_zero_pole('-3', (-3.5, -1), fontsize=12)
    plt.savefig(os.getcwd() + '/PZdiagram.png', bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    pole_zero()
    # contour()
    # pt01()
    # pt02()
    # pt03()
    # pt04()
    # pt05()
    pass

