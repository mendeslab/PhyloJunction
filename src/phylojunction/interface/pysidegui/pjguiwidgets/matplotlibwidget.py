import matplotlib.pyplot as plt # type: ignore

from PySide6.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure # type: ignore

class MatplotlibWidget(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.fig = Figure(figsize=(11,4.5), constrained_layout=True)
        # self.fig = Figure(figsize=(15,6))
        
        # populate self.axes
        self.initialize_axes()
        
        self.canvas = FigureCanvasQTAgg(self.fig)  # widget
        self.canvas.setParent(self)

        # The containing page determines canvas size; a shared canvas must not impose
        # a page-specific minimum that can exceed its available space.
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        fig_layout = QVBoxLayout()
        fig_layout.addWidget(self.canvas)

        self.setLayout(fig_layout)
        

    def initialize_axes(
            self,
            disabled_yticks: bool = True,
            disabled_xticks: bool = True) -> None:
        # Some pages initialize twice while preparing a plot. Reset the figure here
        # so each initialization owns one axes and leaves no old axes or colorbars.
        self.fig.clear()
        # Let Matplotlib allocate room for labels and colorbars as the canvas resizes.
        # Subplots participate in constrained layout; fixed add_axes rectangles do not.
        ax = self.fig.add_subplot(111)
        ax.patch.set_alpha(0.0)

        if disabled_xticks:
            ax.xaxis.set_ticks([])

        if disabled_yticks:
            ax.yaxis.set_ticks([])

        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        self.axes = ax