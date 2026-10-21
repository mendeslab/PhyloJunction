import os

from pathlib import Path
from PySide6.QtWidgets import QPushButton
from PySide6.QtCore import Qt, QRect
from PySide6.QtGui import QPainter, QPixmap, QIcon, QEnterEvent
# from phylojunction.interface.pysidegui.images.icons import resources

my_dir_path = Path(__file__)


class PJPushButton(QPushButton):
    def __init__(self,
                 text="",
                 height=40,
                 minimum_width=50,
                 text_padding=55,
                 text_color="white",
                 icon_path="",
                 icon_color="white",
                 btn_color="#20639B",
                 btn_hover="#1a4971",
                 btn_pressed="#6ed7d3",
                 is_active=False):
        super().__init__()

        self.setText(text)
        self.setMaximumHeight(height)
        self.setMinimumHeight(height)
        self.setCursor(Qt.PointingHandCursor)

        # custom parameters #
        self.minimum_width = minimum_width
        self.text_padding = text_padding
        self.text_color = text_color
        self.icon_path = icon_path
        self.icon_color = icon_color
        self.btn_color = btn_color
        self.btn_hover = btn_hover
        self.btn_pressed = btn_pressed
        self.is_active = is_active

        self.set_style(
            text_padding=self.text_padding,
            text_color=self.text_color,
            btn_color=self.btn_color,
            btn_hover=self.btn_hover,
            btn_pressed=self.btn_pressed,
            is_active=self.is_active
        )


    def set_active(self, is_menu_active):
        self.set_style(
            text_padding=self.text_padding,
            text_color=self.text_color,
            btn_color=self.btn_color,
            btn_hover=self.btn_hover,
            btn_pressed=self.btn_pressed,
            is_active=is_menu_active
        )

    # text_color = "white",
    # btn_color = "#20639B",
    # btn_hover="#1a4971"
    def set_style(
        self,
        text_padding=55,
        text_color="white",
        btn_color="#20639B",
        btn_hover="#1a4971",
        btn_pressed="#6ed7d3",
        is_active=False
    ):

        default_style_str = f"""
        QPushButton {{
            color: {text_color};
            background-color: {btn_color};
            padding-left: {text_padding}px;
            text-align: left;
            border: none;
        }}
        QPushButton:hover {{
            background-color: {btn_hover};
        }}
        QPushButton:pressed {{
            background-color: {btn_pressed};
        }}
        """

        active_style_str = f"""
        QPushButton {{
            background-color: {btn_hover};
            border-right: 5px solid white
        }}
        """

        if not is_active:
            self.setStyleSheet(default_style_str)

        else:
            self.setStyleSheet(default_style_str + active_style_str)


    def paintEvent(self, event):
        # recover button style defined in set_style() #
        QPushButton.paintEvent(self, event)

        # painter #
        qp = QPainter()
        qp.begin(self)
        qp.setRenderHint(QPainter.Antialiasing)
        qp.setPen(Qt.NoPen)  # no border around icon

        # container for icon #
        rect = QRect(0, 0,
                     self.minimum_width, self.height())

        # draw icon #
        self.draw_icon(
            qp,
            self.icon_path,
            rect,
            self.icon_color
        )

        # must end qpainter #
        qp.end()


    def draw_icon(self, qp, image_name, rect, icon_color):
        # get path #
        icon_dir_path = os.path.join(my_dir_path.parent.parent, "images/icons")
        rel_icon_path = os.path.normpath(
            os.path.join(icon_dir_path, image_name))

        # draw icon #
        icon = QPixmap(rel_icon_path)
        painter = QPainter(icon)
        painter.setCompositionMode(QPainter.CompositionMode_SourceIn)
        painter.fillRect(icon.rect(), icon_color)
        qp.drawPixmap(
            ((rect.width() - icon.width()) / 2),  # hor coord
            ((rect.height() - icon.height()) / 2),  # ver coord
            icon  # draw icon
        )
        painter.end()


class PJWhiteBackgroundQPushButton(QPushButton):

    icon_dir_path: str
    icon_normal: QIcon
    icon_over: QIcon
    icon_pressed: QIcon

    def __init__(self, *args, **kwargs):
        # super(ClearPGMPushButton, self).__init__(*a, **kw)
        super().__init__()
        # self.icon_dir_path = os.path.join(my_dir_path.parent.parent, "images/icons")
        self.icon_normal = kwargs.get("icon_normal", "")
        self.icon_over = kwargs.get("icon_over", "")
        self.icon_pressed = kwargs.get("icon_pressed", "")

    def enterEvent(self, event):
        self.setIcon(self.icon_over)
        self.setStyleSheet("color: #ec4a8a;")

    def leaveEvent(self, event):
        self.setIcon(self.icon_normal)
        self.setStyleSheet("color: black;")

    def mousePressEvent(self, event):
        self.setIcon(self.icon_pressed)
        self.setStyleSheet("color: #6ed7d3;")
        
        # need to handle the original event
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event):
        if self.underMouse():
            self.setIcon(self.icon_over)
            self.setStyleSheet("color: #ec4a8a;")
        
        else:
            self.setIcon(self.icon_normal)
            self.setStyleSheet("color: black;")

        # need to handle the original event
        super().mouseReleaseEvent(event)


class PJClearDAGQPushButton(PJWhiteBackgroundQPushButton):

    clear_dag_normal: QIcon
    clear_pgm_over: QIcon
    clear_pgm_pressed: QIcon

    def __init__(self, *a, **kw):
        self.icon_dir_path = \
            os.path.join(my_dir_path.parent.parent, "images/icons/")
        clear_dag_normal = QIcon(QPixmap(":/icon_clear.svg"))
        clear_dag_over = QIcon(QPixmap(":/icon_clear_over.svg"))
        clear_dag_pressed = QIcon(QPixmap(":/icon_clear_pressed.svg"))
        
        super().__init__(
            icon_normal=clear_dag_normal,
            icon_over=clear_dag_over,
            icon_pressed=clear_dag_pressed
        )

    def enterEvent(self, event):
        super().enterEvent(event)

    def leaveEvent(self, event):
        super().leaveEvent(event)

    def mousePressEvent(self, event):
        # need to handle the original event
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event):
        # need to handle the original event
        super().mouseReleaseEvent(event)
