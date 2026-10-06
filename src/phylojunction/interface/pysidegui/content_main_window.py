import os
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QFrame  # type: ignore
from PySide6.QtWidgets import QHBoxLayout  # type: ignore
from PySide6.QtWidgets import QVBoxLayout  # type: ignore
from PySide6.QtWidgets import QStackedWidget  # type: ignore
from PySide6.QtWidgets import QLabel  # type: ignore
from PySide6.QtWidgets import QSpacerItem  # type: ignore
from PySide6.QtWidgets import QSizePolicy  # type: ignore
from PySide6.QtWidgets import QMessageBox  # type: ignore

# pj imports #
from phylojunction.interface.pysidegui.pjguipages.gui_pages \
    import Ui_PJGUIPages  # type: ignore
from phylojunction.interface.pysidegui.pjguiwidgets.pj_buttons \
    import PJPushButton  # type: ignore

__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


my_dir_path = Path(__file__)


class ContentGUIMainWindow(object):
    def setup_ui(self, parent):
        if not parent.objectName():
            parent.setObjectName("MainWindow")

        # window size #
        parent.resize(1200, 720)
        # parent.setMinimumSize(1020, 800)
        parent.setMinimumSize(1020, 785)

        ####################
        # Creating widgets #
        ####################

        # central widget #
        self.central_frame = QFrame()
        self.central_frame.setStyleSheet("background-color: white")

        # main layout #
        self.main_layout = QHBoxLayout(self.central_frame)
        # removes white margins around main_layout content
        self.main_layout.setContentsMargins(0, 0, 0, 0)
        # removes spacing between widgets inside main_layout
        self.main_layout.setSpacing(0)

        # left menu #
        self.left_menu = QFrame()
        self.left_menu.setObjectName("left_menu")
        self.left_menu.setStyleSheet("background-color: #20639B")
        self.left_menu.setMaximumWidth(50)
        self.left_menu.setMinimumWidth(50)

        # left menu layout #
        self.left_menu_layout = QVBoxLayout(self.left_menu)
        self.left_menu_layout.setContentsMargins(0, 0, 0, 0)
        self.left_menu_layout.setSpacing(0)

        # content #
        self.content = QFrame()
        self.content.setStyleSheet("background-color: white")

        # content layout #
        self.content_layout = QVBoxLayout(self.content)  # V: vertical layout
        self.content_layout.setContentsMargins(0, 0, 0, 0)
        self.content_layout.setSpacing(0)

        # top frame menu #
        self.left_menu_top_frame = QFrame()
        self.left_menu_top_frame.setObjectName("left_menu_top_frame")
        self.left_menu_top_frame.setMinimumHeight(50)
        self.left_menu_top_frame.setStyleSheet(
            "#left_menu_top_frame { background-color: #20639B; }")

        # top frame layout #
        self.left_menu_top_frame_layout = QVBoxLayout(self.left_menu_top_frame)
        self.left_menu_top_frame_layout.setContentsMargins(0, 0, 0, 0)
        self.left_menu_top_frame_layout.setSpacing(0)

        # menu spacer #
        self.left_menu_spacer = \
            QSpacerItem(20, 20, QSizePolicy.Minimum, QSizePolicy.Expanding)

        # bottom frame menu #
        self.left_menu_bottom_frame = QFrame()
        self.left_menu_bottom_frame.setObjectName("left_menu_bottom_frame")
        self.left_menu_bottom_frame.setMinimumHeight(50)
        self.left_menu_bottom_frame.setStyleSheet(
            "#left_menu_bottom_frame { background-color: #20639B; }")

        # bottom frame layout #
        self.left_menu_bottom_frame_layout = \
            QVBoxLayout(self.left_menu_bottom_frame)
        self.left_menu_bottom_frame_layout.setContentsMargins(0, 0, 0, 0)
        self.left_menu_bottom_frame_layout.setSpacing(0)

        # version label #
        self.left_menu_version_label = QLabel("v0.0.2")
        self.left_menu_version_label.setAlignment(Qt.AlignCenter)
        self.left_menu_version_label.setMinimumHeight(30)
        self.left_menu_version_label.setMaximumHeight(30)
        self.left_menu_version_label.setStyleSheet("color: white")

        # top bar #
        self.top_bar = QFrame()
        self.top_bar.setMaximumHeight(30)
        self.top_bar.setMinimumHeight(30)
        self.top_bar.setStyleSheet("background-color: #f1f3f5; color: #495057")

        # top bar layout #
        self.top_bar_layout = QHBoxLayout(self.top_bar)
        # add a bit of padding to the left and right
        self.top_bar_layout.setContentsMargins(10, 0, 10, 0)

        # top-left label #
        self.top_label_left = QLabel("MODEL SPECIFICATION")
        self.top_label_left.setStyleSheet(
            "font: 700 10pt 'Segoe UI'")

        # top spacer (like sg.Push in PySimpleGUI, this forces things
        # to max. left and right positioning)
        self.top_spacer = \
            QSpacerItem(20, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        # top-right label #
        self.top_label_right = QLabel("")
        self.top_label_right.setStyleSheet(
            "font: 700 10pt 'Segoe UI'")

        # bottom bar #
        self.bottom_bar = QFrame()
        self.bottom_bar.setMaximumHeight(30)
        self.bottom_bar.setMinimumHeight(30)
        self.bottom_bar.setStyleSheet(
            "background-color: #f1f3f5; color: black")

        # bottom bar layout #
        self.bottom_bar_layout = QHBoxLayout(self.bottom_bar)
        # add a bit of padding to the left and right
        self.bottom_bar_layout.setContentsMargins(10, 0, 10, 0)

        # bottom-left label #
        self.bottom_label_left = QLabel("")
        self.bottom_label_left.setStyleSheet("font: 11pt; color: red")

        # about/licensing #
        self.icon_dir_path = os.path.join(my_dir_path.parent, "images/icons/")
        self.about_licensing = QMessageBox()
        self.about_licensing.setWindowTitle("PhyloJunction")
        self.about_licensing.setText(
            """
            <p style='text-align: center;'>
            <div style='font: 14px;'><b>PhyloJunction</b></div>
            <div style='font: 12px;'>(v0.0.1)</div>
            <br>
            <div style='font: 12px;'><strong>Developed by</strong>
            <div>Fabio K. Mendes and Michael J. Landis</div><br>
            <div style='font: 12px;'><strong>License</strong>
            <div style='font: 12px;'>GNU General Public License v3.0</div>
            """
        )
        self.about_licensing.setIconPixmap(
            QPixmap(self.icon_dir_path + "pj_logo_icon.png").scaled(
                120, 120, Qt.KeepAspectRatio))

        ###########
        # Buttons #
        ###########

        self.menu_button = PJPushButton(
            text="",
            icon_path="icon_menu.svg"
        )

        self.pgm_button = PJPushButton(
            text="Model",
            icon_path="icon_pgm.svg",
            is_active=True
        )

        self.compare_button = PJPushButton(
            text="Compare",
            icon_path="icon_compare.svg"
        )

        self.covg_button = PJPushButton(
            text="Coverage",
            icon_path="icon_covg.svg"
        )

        self.cmd_log_button = PJPushButton(
            text="Command log",
            icon_path="icon_cmd_log.svg"
        )

        # bottom bins #
        self.settings_button = PJPushButton(
            text="Settings",
            icon_path="icon_settings.svg"
        )

        self.warning_button = PJPushButton(
            text="Warnings",
            icon_path="icon_warning.svg"
        )

        #########
        # Pages #
        #########

        # This is a place where Qt Designer can be used to start a
        # QStackedWidget, and then we can get .py code from there,
        # Form > View Code, and paste it here)
        #
        # Once you create a widget with Qt Designer, it might be the case that
        # your Qt Designer won't have 'Form > View Python Code...'
        # If that's the case, just save what you created as a .ui file
        # somewhere, and then on your Terminal, run:
        #
        # $ pyuic6 -x created_file.ui -o created_file.py
        #
        # (the -x flag allows you to run created_file.py immediately to see
        # things)
        #
        # NOTE: the .py out file will have imports involving PyQt6; if PyQt6
        # and PySide6 imports are used in the same project, it seems to cause
        # a 'version `Qt_6.3' not found' error
        #
        # So we fix this by replacing 'from PyQt6 import' with
        # 'from PySide6 import' everywhere

        self.pages = QStackedWidget()
        self.pages.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.pages.setStyleSheet("font-size: 14pt; color: black")

        self.ui_pages = Ui_PJGUIPages()
        self.ui_pages.setupUi(self.pages)  # self.pages is the parent
        self.ui_pages.sample_idx_spin.setDisabled(True)
        self.ui_pages.repl_idx_spin.setDisabled(True)

        # icon_dir_path = os.path.join(my_dir_path.parent, "images/icons/")
        # clear_model_icon = QIcon(icon_dir_path + "icon_covg.svg")
        # self.ui_pages.clear_model.setIconSize(QSize(36, 36))

        ##################
        # Adding widgets #
        ##################

        self.main_layout.addWidget(self.left_menu)
        self.main_layout.addWidget(self.content)

        self.top_bar_layout.addWidget(self.top_label_left)
        # note that spacer is an item
        self.top_bar_layout.addItem(self.top_spacer)
        self.top_bar_layout.addWidget(self.top_label_right)

        self.bottom_bar_layout.addWidget(self.bottom_label_left)

        self.left_menu_top_frame_layout.addWidget(self.menu_button)
        self.left_menu_top_frame_layout.addWidget(self.pgm_button)
        self.left_menu_top_frame_layout.addWidget(self.compare_button)
        self.left_menu_top_frame_layout.addWidget(self.covg_button)
        self.left_menu_top_frame_layout.addWidget(self.cmd_log_button)

        self.left_menu_bottom_frame_layout.addWidget(self.settings_button)
        self.left_menu_bottom_frame_layout.addWidget(self.warning_button)

        self.left_menu_layout.addWidget(self.left_menu_top_frame)
        self.left_menu_layout.addItem(self.left_menu_spacer)
        self.left_menu_layout.addWidget(self.left_menu_bottom_frame)
        self.left_menu_layout.addWidget(self.left_menu_version_label)

        self.content_layout.addWidget(self.top_bar)
        self.content_layout.addWidget(self.pages, stretch=1)
        self.content_layout.addWidget(self.bottom_bar)

        ###################
        # Setting widgets #
        ###################

        # central widget #
        parent.setCentralWidget(self.central_frame)
