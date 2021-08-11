from src.nyquist.NyquistPlot import NyquistPlot
import matplotlib.pyplot as plt


def problem1():
    tf = [[10], [1, 4, 0]]
    nyq = NyquistPlot(tf, xmin=-5, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    nyq.draw_template()
    nyq.draw_contour()
    nyq.draw_roots()
    nyq.annotate_zero_pole('-4', (-4.5, -1), fontsize=12)
    plt.show()


def problem2():
    tf = [[3, 3], [1, 3, 0, 0]]
    nyq = NyquistPlot(tf, xmin=-5, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    nyq.draw_template()
    nyq.draw_contour()
    nyq.draw_roots()
    nyq.annotate_zero_pole('-1', (-1.5, -1), fontsize=12)
    nyq.annotate_zero_pole('-3', (-3.5, -1), fontsize=12)
    plt.show()


def problem3():
    tf = [[10], [1, 6, 13, 0]]
    nyq = NyquistPlot(tf, xmin=-5, xmax=10, ymin=-10, ymax=10, label_axis=True, grid=True, ticks=False)
    nyq.draw_template()
    nyq.draw_contour()
    nyq.draw_roots()
    nyq.annotate_zero_pole('+2j', (-2.5, 2.5), fontsize=12)
    nyq.annotate_zero_pole('-2j', (-2.5, -2.5), fontsize=12)
    nyq.annotate_zero_pole('-3', (-3.5, -0.5), fontsize=12)
    plt.plot([-3, -3], [-2, 2], lw=1, ls='--', color='k')
    plt.plot([0, -3], [2, 2], lw=1, ls='--', color='k')
    plt.plot([0, -3], [-2, -2], lw=1, ls='--', color='k')
    plt.show()



def problem4():
    tf = [[1], [1/9, 0, 1]]
    nyq = NyquistPlot(tf, xmin=-5, xmax=12, ymin=-12, ymax=12, label_axis=True, grid=True, ticks=False)
    nyq.draw_template()
    nyq.draw_contour()
    nyq.draw_roots()
    nyq.annotate_zero_pole('+3j', (-1.5, +2.8), fontsize=12)
    nyq.annotate_zero_pole('-3j', (-1.3, -3.2), fontsize=12)
    plt.show()


if __name__ == '__main__':
    problem1()
    problem2()
    problem3()
    problem4()

