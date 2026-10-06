import matplotlib.pyplot as plt # type: ignore

from PySide6.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure # type: ignore

class MatplotlibWidget(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.fig = Figure(figsize=(11,4.5))
        # self.fig = Figure(figsize=(15,6))
        
        # populate self.axes
        self.initialize_axes(self.fig)
        
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
        # horiz coord of lower-left corner
        # vertical coord of lower-left corner
        # subplott width
        # subplot height        
        ax = self.fig.add_axes([0.075, 0.25, 0.6, 0.7])
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