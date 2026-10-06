# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'pjgui_pages.ui'
##
## Created by: Qt User Interface Compiler version 6.11.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QCheckBox, QFrame, QGridLayout,
    QHBoxLayout, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QPushButton, QRadioButton, QSizePolicy,
    QSpacerItem, QSpinBox, QStackedWidget, QTabWidget,
    QTextEdit, QVBoxLayout, QWidget)

from phylojunction.interface.pysidegui.images.icons import resources
from phylojunction.interface.pysidegui.pjguiwidgets.matplotlibwidget import MatplotlibWidget
from phylojunction.interface.pysidegui.pjguiwidgets.pj_buttons import (PJClearDAGQPushButton, PJReDrawQPushButton)

class Ui_PJGUIPages(object):
    def setupUi(self, PJGUIPages):
        if not PJGUIPages.objectName():
            PJGUIPages.setObjectName(u"PJGUIPages")
        PJGUIPages.resize(980, 700)
        PJGUIPages.setMinimumSize(QSize(960, 650))
        self.settings_page = QWidget()
        self.settings_page.setObjectName(u"settings_page")
        self.gridLayoutWidget_4 = QWidget(self.settings_page)
        self.gridLayoutWidget_4.setObjectName(u"gridLayoutWidget_4")
        self.gridLayoutWidget_4.setGeometry(QRect(10, 10, 961, 681))
        self.settings_page_grid_layout = QGridLayout(self.gridLayoutWidget_4)
        self.settings_page_grid_layout.setObjectName(u"settings_page_grid_layout")
        self.settings_page_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.settings_vert_spacer = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.settings_page_grid_layout.addItem(self.settings_vert_spacer, 7, 0, 1, 1)

        self.filename_prefix_label = QLabel(self.gridLayoutWidget_4)
        self.filename_prefix_label.setObjectName(u"filename_prefix_label")
        self.filename_prefix_label.setMinimumSize(QSize(110, 24))
        self.filename_prefix_label.setMaximumSize(QSize(110, 24))

        self.settings_page_grid_layout.addWidget(self.filename_prefix_label, 3, 0, 1, 1)

        self.filename_prefix_textbox = QTextEdit(self.gridLayoutWidget_4)
        self.filename_prefix_textbox.setObjectName(u"filename_prefix_textbox")
        self.filename_prefix_textbox.setMinimumSize(QSize(0, 24))
        self.filename_prefix_textbox.setMaximumSize(QSize(16777215, 24))
        self.filename_prefix_textbox.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.filename_prefix_textbox.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.settings_page_grid_layout.addWidget(self.filename_prefix_textbox, 3, 1, 1, 1)

        self.tree_dag_node_label = QLabel(self.gridLayoutWidget_4)
        self.tree_dag_node_label.setObjectName(u"tree_dag_node_label")

        self.settings_page_grid_layout.addWidget(self.tree_dag_node_label, 5, 0, 1, 1)

        self.random_seed_prefix_label = QLabel(self.gridLayoutWidget_4)
        self.random_seed_prefix_label.setObjectName(u"random_seed_prefix_label")
        self.random_seed_prefix_label.setMinimumSize(QSize(110, 24))
        self.random_seed_prefix_label.setMaximumSize(QSize(110, 24))

        self.settings_page_grid_layout.addWidget(self.random_seed_prefix_label, 1, 0, 1, 1)

        self.stoch_map_label = QLabel(self.gridLayoutWidget_4)
        self.stoch_map_label.setObjectName(u"stoch_map_label")
        self.stoch_map_label.setMinimumSize(QSize(0, 24))
        self.stoch_map_label.setMaximumSize(QSize(16777215, 24))
        font = QFont()
        font.setBold(True)
        self.stoch_map_label.setFont(font)

        self.settings_page_grid_layout.addWidget(self.stoch_map_label, 4, 0, 1, 1)

        self.settings_label = QLabel(self.gridLayoutWidget_4)
        self.settings_label.setObjectName(u"settings_label")
        self.settings_label.setMinimumSize(QSize(0, 24))
        self.settings_label.setMaximumSize(QSize(16777215, 24))
        self.settings_label.setFont(font)

        self.settings_page_grid_layout.addWidget(self.settings_label, 2, 0, 1, 1)

        self.random_seed_prefix_textbox = QTextEdit(self.gridLayoutWidget_4)
        self.random_seed_prefix_textbox.setObjectName(u"random_seed_prefix_textbox")
        self.random_seed_prefix_textbox.setMinimumSize(QSize(0, 24))
        self.random_seed_prefix_textbox.setMaximumSize(QSize(16777215, 24))
        self.random_seed_prefix_textbox.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.random_seed_prefix_textbox.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.settings_page_grid_layout.addWidget(self.random_seed_prefix_textbox, 1, 1, 1, 1)

        self.line1 = QFrame(self.gridLayoutWidget_4)
        self.line1.setObjectName(u"line1")
        self.line1.setFrameShape(QFrame.Shape.HLine)
        self.line1.setFrameShadow(QFrame.Shadow.Sunken)

        self.settings_page_grid_layout.addWidget(self.line1, 0, 1, 1, 1)

        self.line2 = QFrame(self.gridLayoutWidget_4)
        self.line2.setObjectName(u"line2")
        self.line2.setFrameShape(QFrame.Shape.HLine)
        self.line2.setFrameShadow(QFrame.Shadow.Sunken)

        self.settings_page_grid_layout.addWidget(self.line2, 2, 1, 1, 1)

        self.line3 = QFrame(self.gridLayoutWidget_4)
        self.line3.setObjectName(u"line3")
        self.line3.setFrameShape(QFrame.Shape.HLine)
        self.line3.setFrameShadow(QFrame.Shadow.Sunken)

        self.settings_page_grid_layout.addWidget(self.line3, 4, 1, 1, 1)

        self.sim_config_label = QLabel(self.gridLayoutWidget_4)
        self.sim_config_label.setObjectName(u"sim_config_label")
        self.sim_config_label.setMinimumSize(QSize(0, 24))
        self.sim_config_label.setFont(font)

        self.settings_page_grid_layout.addWidget(self.sim_config_label, 0, 0, 1, 1)

        self.attr_name_label = QLabel(self.gridLayoutWidget_4)
        self.attr_name_label.setObjectName(u"attr_name_label")
        self.attr_name_label.setMinimumSize(QSize(0, 24))
        self.attr_name_label.setMaximumSize(QSize(16777215, 24))

        self.settings_page_grid_layout.addWidget(self.attr_name_label, 6, 0, 1, 1)

        self.tr_dag_node_textbox = QTextEdit(self.gridLayoutWidget_4)
        self.tr_dag_node_textbox.setObjectName(u"tr_dag_node_textbox")
        self.tr_dag_node_textbox.setMinimumSize(QSize(0, 24))
        self.tr_dag_node_textbox.setMaximumSize(QSize(16777215, 24))
        self.tr_dag_node_textbox.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.tr_dag_node_textbox.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.settings_page_grid_layout.addWidget(self.tr_dag_node_textbox, 5, 1, 1, 1)

        self.attr_name_textbox = QTextEdit(self.gridLayoutWidget_4)
        self.attr_name_textbox.setObjectName(u"attr_name_textbox")
        self.attr_name_textbox.setMinimumSize(QSize(0, 24))
        self.attr_name_textbox.setMaximumSize(QSize(16777215, 24))
        self.attr_name_textbox.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.attr_name_textbox.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.settings_page_grid_layout.addWidget(self.attr_name_textbox, 6, 1, 1, 1)

        PJGUIPages.addWidget(self.settings_page)
        self.pgm_page = QWidget()
        self.pgm_page.setObjectName(u"pgm_page")
        self.verticalLayout = QVBoxLayout(self.pgm_page)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.pgm_page_frame = QFrame(self.pgm_page)
        self.pgm_page_frame.setObjectName(u"pgm_page_frame")
        self.pgm_page_frame.setStyleSheet(u"QFrame#pgm_page_frame {\n"
"background-color: white;\n"
"border: 0;\n"
"}")
        self.pgm_page_frame.setFrameShape(QFrame.StyledPanel)
        self.pgm_page_frame.setFrameShadow(QFrame.Raised)
        self.model_frame_layout = QVBoxLayout(self.pgm_page_frame)
        self.model_frame_layout.setSpacing(6)
        self.model_frame_layout.setObjectName(u"model_frame_layout")
        self.model_frame_layout.setContentsMargins(6, 6, 6, 6)
        self.node_content_tabs = QTabWidget(self.pgm_page_frame)
        self.node_content_tabs.setObjectName(u"node_content_tabs")
        self.node_content_tabs.setStyleSheet(u"color: black;\n"
"font: 14pt ;\n"
"")
        self.node_content_tabs.setTabShape(QTabWidget.Rounded)
        self.node_content_tabs.setIconSize(QSize(16, 16))
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.node_content_tabs.sizePolicy().hasHeightForWidth())
        self.node_content_tabs.setSizePolicy(sizePolicy)
        self.values_tab = QWidget()
        self.values_tab.setObjectName(u"values_tab")
        self.values_tab_layout = QVBoxLayout(self.values_tab)
        self.values_tab_layout.setObjectName(u"values_tab_layout")
        self.values_tab_layout.setContentsMargins(0, 0, 0, 0)
        self.values_content = QTextEdit(self.values_tab)
        self.values_content.setObjectName(u"values_content")
        self.values_content.setAcceptDrops(False)
        self.values_content.setStyleSheet(u"border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"font: 11pt \"Courier\";")
        self.values_content.setReadOnly(True)

        self.values_tab_layout.addWidget(self.values_content)

        self.node_content_tabs.addTab(self.values_tab, "")
        self.summary_tab = QWidget()
        self.summary_tab.setObjectName(u"summary_tab")
        self.summary_tab_layout = QVBoxLayout(self.summary_tab)
        self.summary_tab_layout.setObjectName(u"summary_tab_layout")
        self.summary_tab_layout.setContentsMargins(0, 0, 0, 0)
        self.summary_content = QTextEdit(self.summary_tab)
        self.summary_content.setObjectName(u"summary_content")
        self.summary_content.setStyleSheet(u"border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"font: 11pt \"Courier\";")

        self.summary_tab_layout.addWidget(self.summary_content)

        self.node_content_tabs.addTab(self.summary_tab, "")

        self.model_frame_layout.addWidget(self.node_content_tabs)

        self.plot_controls = QWidget(self.pgm_page_frame)
        self.plot_controls.setObjectName(u"plot_controls")
        self.one_sample_radio = QRadioButton(self.plot_controls)
        self.one_sample_radio.setObjectName(u"one_sample_radio")
        self.one_sample_radio.setEnabled(False)
        self.one_sample_radio.setStyleSheet(u"color: black;\n"
"")
        self.one_sample_radio.setCheckable(False)
        self.one_sample_radio.setChecked(False)
        self.all_samples_radio = QRadioButton(self.plot_controls)
        self.all_samples_radio.setObjectName(u"all_samples_radio")
        self.all_samples_radio.setEnabled(False)
        self.all_samples_radio.setStyleSheet(u"color: black;")
        self.all_samples_radio.setCheckable(False)
        self.all_samples_radio.setAutoExclusive(True)
        self.reconstructed_tree_check = QCheckBox(self.plot_controls)
        self.reconstructed_tree_check.setObjectName(u"reconstructed_tree_check")
        self.sample_idx_spin = QSpinBox(self.plot_controls)
        self.sample_idx_spin.setObjectName(u"sample_idx_spin")
        self.sample_idx_spin.setEnabled(False)
        self.sample_idx_spin.setStyleSheet(u"color: black;")
        self.sample_idx_spin.setAccelerated(True)
        self.repl_idx_spin = QSpinBox(self.plot_controls)
        self.repl_idx_spin.setObjectName(u"repl_idx_spin")
        self.repl_idx_spin.setEnabled(False)
        self.repl_idx_spin.setStyleSheet(u"color: black;\n"
"")
        self.repl_idx_spin.setReadOnly(False)
        self.repl_idx_spin.setAccelerated(True)
        self.save_pgm_node_plot = QPushButton(self.plot_controls)
        self.save_pgm_node_plot.setObjectName(u"save_pgm_node_plot")
        self.save_pgm_node_plot.setEnabled(True)
        self.save_pgm_node_plot.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.save_pgm_node_plot.setStyleSheet(u"QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")

        self.model_frame_layout.addWidget(self.plot_controls)

        self.plot_and_nodes_layout = QHBoxLayout()
        self.plot_and_nodes_layout.setObjectName(u"plot_and_nodes_layout")
        self.pgm_page_matplotlib_widget = MatplotlibWidget(self.pgm_page_frame)
        self.pgm_page_matplotlib_widget.setObjectName(u"pgm_page_matplotlib_widget")
        sizePolicy1 = QSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.pgm_page_matplotlib_widget.sizePolicy().hasHeightForWidth())
        self.pgm_page_matplotlib_widget.setSizePolicy(sizePolicy1)
        self.pgm_page_matplotlib_widget.setMinimumSize(QSize(320, 240))

        self.plot_and_nodes_layout.addWidget(self.pgm_page_matplotlib_widget)

        self.node_list_vert_layout = QVBoxLayout()
        self.node_list_vert_layout.setObjectName(u"node_list_vert_layout")
        self.model_label = QLabel(self.pgm_page_frame)
        self.model_label.setObjectName(u"model_label")
        self.model_label.setStyleSheet(u"font: 14pt \"Ubuntu\";\n"
"color: black;\n"
"")
        self.model_label.setAlignment(Qt.AlignCenter)

        self.node_list_vert_layout.addWidget(self.model_label)

        self.node_list = QListWidget(self.pgm_page_frame)
        self.node_list.setObjectName(u"node_list")
        self.node_list.setStyleSheet(u"QListWidget {\n"
"	background-color: #f1f3f5;\n"
"	color: #495057;\n"
"	border-radius: 5px;\n"
"	border: 0;\n"
"}\n"
"QListWidget::item:selected{\n"
"	background-color: #495057;\n"
"	color: white;\n"
"	border: 0;\n"
"}\n"
"")

        self.node_list_vert_layout.addWidget(self.node_list)

        self.redraw_node = PJReDrawQPushButton(self.pgm_page_frame)
        self.redraw_node.setObjectName(u"redraw_node")
        self.redraw_node.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.redraw_node.setMouseTracking(True)
        icon = QIcon()
        icon.addFile(u":/draw.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.redraw_node.setIcon(icon)

        self.node_list_vert_layout.addWidget(self.redraw_node)

        self.clear_model = PJClearDAGQPushButton(self.pgm_page_frame)
        self.clear_model.setObjectName(u"clear_model")
        self.clear_model.setEnabled(True)
        self.clear_model.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.clear_model.setStyleSheet(u"QPushButton:hover {\n"
"    color: #ec4a8a;\n"
"}")
        icon1 = QIcon()
        icon1.addFile(u":/icon_clear.svg", QSize(), QIcon.Mode.Normal, QIcon.State.Off)
        self.clear_model.setIcon(icon1)
        self.clear_model.setIconSize(QSize(20, 20))

        self.node_list_vert_layout.addWidget(self.clear_model)

        self.node_list_vert_layout.setStretch(1, 1)

        self.plot_and_nodes_layout.addLayout(self.node_list_vert_layout)

        self.plot_and_nodes_layout.setStretch(0, 1)

        self.model_frame_layout.addLayout(self.plot_and_nodes_layout)

        self.cmd_prompt_vert_layout = QVBoxLayout()
        self.cmd_prompt_vert_layout.setObjectName(u"cmd_prompt_vert_layout")
        self.cmd_prompt_label = QLabel(self.pgm_page_frame)
        self.cmd_prompt_label.setObjectName(u"cmd_prompt_label")
        self.cmd_prompt_label.setStyleSheet(u"font: 14pt \"Ubuntu\";\n"
"color: black;")
        self.cmd_prompt_label.setAlignment(Qt.AlignBottom|Qt.AlignLeading|Qt.AlignLeft)

        self.cmd_prompt_vert_layout.addWidget(self.cmd_prompt_label)

        self.cmd_prompt = QLineEdit(self.pgm_page_frame)
        self.cmd_prompt.setObjectName(u"cmd_prompt")
        font1 = QFont()
        font1.setFamilies([u"Courier"])
        font1.setPointSize(14)
        font1.setBold(True)
        font1.setItalic(False)
        self.cmd_prompt.setFont(font1)
        self.cmd_prompt.setStyleSheet(u"background-color: black;\n"
"font: 14pt \"Courier\";\n"
"font-weight: bold;\n"
"color: white;\n"
"border-radius: 5px;")
        self.cmd_prompt.setCursorPosition(0)
        self.cmd_prompt.setAlignment(Qt.AlignLeading|Qt.AlignLeft|Qt.AlignVCenter)
        self.cmd_prompt.setDragEnabled(True)

        self.cmd_prompt_vert_layout.addWidget(self.cmd_prompt)


        self.model_frame_layout.addLayout(self.cmd_prompt_vert_layout)

        self.model_frame_layout.setStretch(2, 1)

        self.verticalLayout.addWidget(self.pgm_page_frame)

        PJGUIPages.addWidget(self.pgm_page)
        self.compare_page = QWidget()
        self.compare_page.setObjectName(u"compare_page")
        self.compare_page.setEnabled(True)
        self.gridLayoutWidget_2 = QWidget(self.compare_page)
        self.gridLayoutWidget_2.setObjectName(u"gridLayoutWidget_2")
        self.gridLayoutWidget_2.setGeometry(QRect(10, 10, 961, 681))
        self.compare_page_grid_layout = QGridLayout(self.gridLayoutWidget_2)
        self.compare_page_grid_layout.setObjectName(u"compare_page_grid_layout")
        self.compare_page_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.node_stat_vert_layout = QVBoxLayout()
        self.node_stat_vert_layout.setObjectName(u"node_stat_vert_layout")
        self.compare_node_frame = QFrame(self.gridLayoutWidget_2)
        self.compare_node_frame.setObjectName(u"compare_node_frame")
        self.compare_node_frame.setMinimumSize(QSize(195, 350))
        self.compare_node_frame.setMaximumSize(QSize(195, 350))
        self.compare_node_frame.setStyleSheet(u"border-radius: 5px;\n"
"background-color: #e6f4f4;\n"
"border: 0;")
        self.compare_node_frame.setFrameShape(QFrame.StyledPanel)
        self.compare_node_frame.setFrameShadow(QFrame.Raised)
        self.compare_node_list = QListWidget(self.compare_node_frame)
        self.compare_node_list.setObjectName(u"compare_node_list")
        self.compare_node_list.setGeometry(QRect(10, 30, 172, 280))
        self.compare_node_list.setMinimumSize(QSize(172, 280))
        self.compare_node_list.setMaximumSize(QSize(172, 280))
        self.compare_node_list.setStyleSheet(u"QListWidget {\n"
"	background-color: #f2f9f9;\n"
"	color: #495057;\n"
"	border: 0;\n"
"    border-radius: 0px;\n"
"}\n"
"QListWidget::item:selected{\n"
"	background-color: #495057;\n"
"	color: white;\n"
"	border: 0;\n"
"}")
        self.compare_node_label = QLabel(self.compare_node_frame)
        self.compare_node_label.setObjectName(u"compare_node_label")
        self.compare_node_label.setGeometry(QRect(20, 10, 155, 16))
        self.compare_node_label.setMinimumSize(QSize(155, 16))
        self.compare_node_label.setMaximumSize(QSize(155, 16))
        self.compare_node_label.setAlignment(Qt.AlignCenter)
        self.avg_replicate_check_button = QCheckBox(self.compare_node_frame)
        self.avg_replicate_check_button.setObjectName(u"avg_replicate_check_button")
        self.avg_replicate_check_button.setGeometry(QRect(40, 320, 155, 20))
        self.avg_replicate_check_button.setMinimumSize(QSize(155, 20))
        self.avg_replicate_check_button.setMaximumSize(QSize(155, 20))
        font2 = QFont()
        font2.setPointSize(9)
        font2.setItalic(False)
        self.avg_replicate_check_button.setFont(font2)
        self.avg_replicate_check_button.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))

        self.node_stat_vert_layout.addWidget(self.compare_node_frame)

        self.compare_vert_spacer = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.node_stat_vert_layout.addItem(self.compare_vert_spacer)

        self.compare_stats_label = QLabel(self.gridLayoutWidget_2)
        self.compare_stats_label.setObjectName(u"compare_stats_label")
        self.compare_stats_label.setMinimumSize(QSize(0, 16))
        self.compare_stats_label.setMaximumSize(QSize(16777215, 16))

        self.node_stat_vert_layout.addWidget(self.compare_stats_label, 0, Qt.AlignHCenter)

        self.summary_stats_list = QListWidget(self.gridLayoutWidget_2)
        self.summary_stats_list.setObjectName(u"summary_stats_list")
        self.summary_stats_list.setMinimumSize(QSize(195, 250))
        self.summary_stats_list.setMaximumSize(QSize(195, 250))
        self.summary_stats_list.setStyleSheet(u"QListWidget{\n"
"    background-color: #f7f7f7;\n"
"    border: 2px solid lightgray;\n"
"}\n"
"QListWidget::item:selected{\n"
"	background-color: #495057;\n"
"	color: white;\n"
"	border: 0;\n"
"}")

        self.node_stat_vert_layout.addWidget(self.summary_stats_list, 0, Qt.AlignHCenter)


        self.compare_page_grid_layout.addLayout(self.node_stat_vert_layout, 0, 0, 1, 1)

        self.csv_violinplot_vert_layout = QVBoxLayout()
        self.csv_violinplot_vert_layout.setObjectName(u"csv_violinplot_vert_layout")
        self.compare_csv_button = QPushButton(self.gridLayoutWidget_2)
        self.compare_csv_button.setObjectName(u"compare_csv_button")
        self.compare_csv_button.setMinimumSize(QSize(170, 24))
        self.compare_csv_button.setMaximumSize(QSize(170, 24))
        self.compare_csv_button.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.compare_csv_button.setMouseTracking(True)
        self.compare_csv_button.setStyleSheet(u"QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")

        self.csv_violinplot_vert_layout.addWidget(self.compare_csv_button, 0, Qt.AlignHCenter)

        self.compare_csv_textbox = QTextEdit(self.gridLayoutWidget_2)
        self.compare_csv_textbox.setObjectName(u"compare_csv_textbox")
        self.compare_csv_textbox.setMinimumSize(QSize(0, 188))
        self.compare_csv_textbox.setMaximumSize(QSize(16777215, 188))
        self.compare_csv_textbox.setStyleSheet(u"background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"font: 11pt \"Courier\";")
        self.compare_csv_textbox.setReadOnly(True)

        self.csv_violinplot_vert_layout.addWidget(self.compare_csv_textbox)

        self.compare_page_matplotlib_widget = MatplotlibWidget(self.gridLayoutWidget_2)
        self.compare_page_matplotlib_widget.setObjectName(u"compare_page_matplotlib_widget")
        self.compare_page_matplotlib_widget.setMinimumSize(QSize(740, 400))
        self.compare_page_matplotlib_widget.setMaximumSize(QSize(740, 400))

        self.csv_violinplot_vert_layout.addWidget(self.compare_page_matplotlib_widget)

        self.draw_save_hor_layout = QHBoxLayout()
        self.draw_save_hor_layout.setObjectName(u"draw_save_hor_layout")
        self.draw_violins_button = QPushButton(self.gridLayoutWidget_2)
        self.draw_violins_button.setObjectName(u"draw_violins_button")
        self.draw_violins_button.setMinimumSize(QSize(95, 24))
        self.draw_violins_button.setMaximumSize(QSize(95, 24))
        self.draw_violins_button.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.draw_violins_button.setStyleSheet(u"QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")

        self.draw_save_hor_layout.addWidget(self.draw_violins_button)

        self.save_violins_button = QPushButton(self.gridLayoutWidget_2)
        self.save_violins_button.setObjectName(u"save_violins_button")
        self.save_violins_button.setEnabled(True)
        self.save_violins_button.setMinimumSize(QSize(130, 24))
        self.save_violins_button.setMaximumSize(QSize(130, 24))
        self.save_violins_button.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        self.save_violins_button.setStyleSheet(u"QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")

        self.draw_save_hor_layout.addWidget(self.save_violins_button)


        self.csv_violinplot_vert_layout.addLayout(self.draw_save_hor_layout)


        self.compare_page_grid_layout.addLayout(self.csv_violinplot_vert_layout, 0, 1, 1, 1)

        PJGUIPages.addWidget(self.compare_page)
        self.coverage_page = QWidget()
        self.coverage_page.setObjectName(u"coverage_page")
        self.gridLayoutWidget_3 = QWidget(self.coverage_page)
        self.gridLayoutWidget_3.setObjectName(u"gridLayoutWidget_3")
        self.gridLayoutWidget_3.setGeometry(QRect(10, 10, 963, 684))
        self.coverage_page_grid_layout = QGridLayout(self.gridLayoutWidget_3)
        self.coverage_page_grid_layout.setObjectName(u"coverage_page_grid_layout")
        self.coverage_page_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.csv_covplot_vert_layout = QVBoxLayout()
        self.csv_covplot_vert_layout.setObjectName(u"csv_covplot_vert_layout")
        self.cov_top_hor_layout = QHBoxLayout()
        self.cov_top_hor_layout.setObjectName(u"cov_top_hor_layout")
        self.read_hpd_csv_button = QPushButton(self.gridLayoutWidget_3)
        self.read_hpd_csv_button.setObjectName(u"read_hpd_csv_button")
        self.read_hpd_csv_button.setMinimumSize(QSize(170, 24))
        self.read_hpd_csv_button.setMaximumSize(QSize(170, 24))
        self.read_hpd_csv_button.setStyleSheet(u"QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")

        self.cov_top_hor_layout.addWidget(self.read_hpd_csv_button, 0, Qt.AlignHCenter)

        self.read_logfile_button = QPushButton(self.gridLayoutWidget_3)
        self.read_logfile_button.setObjectName(u"read_logfile_button")
        self.read_logfile_button.setMinimumSize(QSize(150, 24))
        self.read_logfile_button.setMaximumSize(QSize(150, 24))
        self.read_logfile_button.setStyleSheet(u"QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")

        self.cov_top_hor_layout.addWidget(self.read_logfile_button, 0, Qt.AlignHCenter)

        self.node_name_label = QLabel(self.gridLayoutWidget_3)
        self.node_name_label.setObjectName(u"node_name_label")
        self.node_name_label.setMinimumSize(QSize(150, 24))
        self.node_name_label.setMaximumSize(QSize(150, 24))

        self.cov_top_hor_layout.addWidget(self.node_name_label)

        self.log_node_name_lineedit = QLineEdit(self.gridLayoutWidget_3)
        self.log_node_name_lineedit.setObjectName(u"log_node_name_lineedit")
        self.log_node_name_lineedit.setMinimumSize(QSize(150, 24))
        self.log_node_name_lineedit.setMaximumSize(QSize(0, 24))

        self.cov_top_hor_layout.addWidget(self.log_node_name_lineedit)


        self.csv_covplot_vert_layout.addLayout(self.cov_top_hor_layout)

        self.cov_grid_layout = QGridLayout()
        self.cov_grid_layout.setObjectName(u"cov_grid_layout")
        self.coverage_csv_textbox = QTextEdit(self.gridLayoutWidget_3)
        self.coverage_csv_textbox.setObjectName(u"coverage_csv_textbox")
        self.coverage_csv_textbox.setMinimumSize(QSize(570, 170))
        self.coverage_csv_textbox.setMaximumSize(QSize(570, 170))
        self.coverage_csv_textbox.setStyleSheet(u"background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"font: 11pt \"Courier\";")
        self.coverage_csv_textbox.setReadOnly(True)

        self.cov_grid_layout.addWidget(self.coverage_csv_textbox, 1, 0, 1, 1)

        self.coverage_textbox = QListWidget(self.gridLayoutWidget_3)
        self.coverage_textbox.setObjectName(u"coverage_textbox")
        self.coverage_textbox.setMinimumSize(QSize(170, 170))
        self.coverage_textbox.setMaximumSize(QSize(170, 170))
        self.coverage_textbox.setStyleSheet(u"background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;")

        self.cov_grid_layout.addWidget(self.coverage_textbox, 1, 1, 1, 1)

        self.coverage_label = QLabel(self.gridLayoutWidget_3)
        self.coverage_label.setObjectName(u"coverage_label")

        self.cov_grid_layout.addWidget(self.coverage_label, 0, 1, 1, 1, Qt.AlignHCenter)


        self.csv_covplot_vert_layout.addLayout(self.cov_grid_layout)

        self.coverage_page_matplotlib_widget = MatplotlibWidget(self.gridLayoutWidget_3)
        self.coverage_page_matplotlib_widget.setObjectName(u"coverage_page_matplotlib_widget")
        self.coverage_page_matplotlib_widget.setMinimumSize(QSize(740, 400))
        self.coverage_page_matplotlib_widget.setMaximumSize(QSize(740, 400))

        self.csv_covplot_vert_layout.addWidget(self.coverage_page_matplotlib_widget)

        self.cov_draw_save_hor_layout = QHBoxLayout()
        self.cov_draw_save_hor_layout.setObjectName(u"cov_draw_save_hor_layout")
        self.draw_cov_button = QPushButton(self.gridLayoutWidget_3)
        self.draw_cov_button.setObjectName(u"draw_cov_button")
        self.draw_cov_button.setMinimumSize(QSize(95, 24))
        self.draw_cov_button.setMaximumSize(QSize(95, 24))
        self.draw_cov_button.setStyleSheet(u"QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")

        self.cov_draw_save_hor_layout.addWidget(self.draw_cov_button, 0, Qt.AlignHCenter)

        self.save_cov_button = QPushButton(self.gridLayoutWidget_3)
        self.save_cov_button.setObjectName(u"save_cov_button")
        self.save_cov_button.setMinimumSize(QSize(130, 24))
        self.save_cov_button.setMaximumSize(QSize(130, 24))
        self.save_cov_button.setStyleSheet(u"QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")

        self.cov_draw_save_hor_layout.addWidget(self.save_cov_button, 0, Qt.AlignHCenter)


        self.csv_covplot_vert_layout.addLayout(self.cov_draw_save_hor_layout)


        self.coverage_page_grid_layout.addLayout(self.csv_covplot_vert_layout, 0, 1, 1, 1)

        self.cov_node_stat_vert_layout = QVBoxLayout()
        self.cov_node_stat_vert_layout.setObjectName(u"cov_node_stat_vert_layout")
        self.coverage_frame = QFrame(self.gridLayoutWidget_3)
        self.coverage_frame.setObjectName(u"coverage_frame")
        self.coverage_frame.setMinimumSize(QSize(195, 320))
        self.coverage_frame.setMaximumSize(QSize(195, 320))
        self.coverage_frame.setStyleSheet(u"border-radius: 5px;\n"
"background-color: #fcf5e3;\n"
"border: 0;")
        self.coverage_frame.setFrameShape(QFrame.StyledPanel)
        self.coverage_frame.setFrameShadow(QFrame.Raised)
        self.coverage_node_list = QListWidget(self.coverage_frame)
        self.coverage_node_list.setObjectName(u"coverage_node_list")
        self.coverage_node_list.setGeometry(QRect(10, 30, 172, 280))
        self.coverage_node_list.setMinimumSize(QSize(172, 280))
        self.coverage_node_list.setMaximumSize(QSize(172, 280))
        self.coverage_node_list.setStyleSheet(u"QListWidget {\n"
"	background-color: #f7f5f0;\n"
"	color: #495057;\n"
"	border: 0;\n"
"    border-radius: 0px;\n"
"}\n"
"QListWidget::item:selected{\n"
"	background-color: #495057;\n"
"	color: white;\n"
"	border: 0;\n"
"}")
        self.coverage_node_label = QLabel(self.coverage_frame)
        self.coverage_node_label.setObjectName(u"coverage_node_label")
        self.coverage_node_label.setGeometry(QRect(20, 10, 140, 16))
        self.coverage_node_label.setMinimumSize(QSize(140, 16))
        self.coverage_node_label.setMaximumSize(QSize(140, 16))
        self.coverage_node_label.setAlignment(Qt.AlignCenter)

        self.cov_node_stat_vert_layout.addWidget(self.coverage_frame, 0, Qt.AlignHCenter)

        self.cov_vert_spacer = QSpacerItem(20, 40, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.cov_node_stat_vert_layout.addItem(self.cov_vert_spacer)

        self.compare_stats_label_2 = QLabel(self.gridLayoutWidget_3)
        self.compare_stats_label_2.setObjectName(u"compare_stats_label_2")
        self.compare_stats_label_2.setMinimumSize(QSize(0, 16))
        self.compare_stats_label_2.setMaximumSize(QSize(16777215, 16))

        self.cov_node_stat_vert_layout.addWidget(self.compare_stats_label_2, 0, Qt.AlignHCenter)

        self.cov_summary_stats_list = QListWidget(self.gridLayoutWidget_3)
        self.cov_summary_stats_list.setObjectName(u"cov_summary_stats_list")
        self.cov_summary_stats_list.setMinimumSize(QSize(195, 250))
        self.cov_summary_stats_list.setMaximumSize(QSize(195, 250))
        self.cov_summary_stats_list.setStyleSheet(u"QListWidget{\n"
"    background-color: #f7f7f7;\n"
"    border: 2px solid lightgray;\n"
"}\n"
"QListWidget::item:selected{\n"
"	background-color: #495057;\n"
"	color: white;\n"
"	border: 0;\n"
"}")

        self.cov_node_stat_vert_layout.addWidget(self.cov_summary_stats_list, 0, Qt.AlignHCenter)


        self.coverage_page_grid_layout.addLayout(self.cov_node_stat_vert_layout, 0, 0, 1, 1)

        PJGUIPages.addWidget(self.coverage_page)
        self.cmd_log_page = QWidget()
        self.cmd_log_page.setObjectName(u"cmd_log_page")
        self.verticalLayoutWidget_3 = QWidget(self.cmd_log_page)
        self.verticalLayoutWidget_3.setObjectName(u"verticalLayoutWidget_3")
        self.verticalLayoutWidget_3.setGeometry(QRect(10, 10, 962, 681))
        self.cmd_log_page_frame = QVBoxLayout(self.verticalLayoutWidget_3)
        self.cmd_log_page_frame.setObjectName(u"cmd_log_page_frame")
        self.cmd_log_page_frame.setContentsMargins(0, 0, 0, 0)
        self.cmd_log_textbox = QTextEdit(self.verticalLayoutWidget_3)
        self.cmd_log_textbox.setObjectName(u"cmd_log_textbox")
        self.cmd_log_textbox.setEnabled(True)
        self.cmd_log_textbox.setMinimumSize(QSize(960, 650))
        self.cmd_log_textbox.setMaximumSize(QSize(16777215, 650))
        self.cmd_log_textbox.setAcceptDrops(False)
        self.cmd_log_textbox.setStyleSheet(u"background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;")
        self.cmd_log_textbox.setReadOnly(True)

        self.cmd_log_page_frame.addWidget(self.cmd_log_textbox)

        PJGUIPages.addWidget(self.cmd_log_page)
        self.warnings_page = QWidget()
        self.warnings_page.setObjectName(u"warnings_page")
        self.verticalLayoutWidget_2 = QWidget(self.warnings_page)
        self.verticalLayoutWidget_2.setObjectName(u"verticalLayoutWidget_2")
        self.verticalLayoutWidget_2.setGeometry(QRect(10, 10, 962, 681))
        self.warnings_page_frame = QVBoxLayout(self.verticalLayoutWidget_2)
        self.warnings_page_frame.setObjectName(u"warnings_page_frame")
        self.warnings_page_frame.setContentsMargins(0, 0, 0, 0)
        self.warnings_textbox = QTextEdit(self.verticalLayoutWidget_2)
        self.warnings_textbox.setObjectName(u"warnings_textbox")
        self.warnings_textbox.setMinimumSize(QSize(960, 650))
        self.warnings_textbox.setMaximumSize(QSize(16777215, 650))
        self.warnings_textbox.setStyleSheet(u"background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"color: red;")
        self.warnings_textbox.setReadOnly(True)

        self.warnings_page_frame.addWidget(self.warnings_textbox)

        PJGUIPages.addWidget(self.warnings_page)

        self.retranslateUi(PJGUIPages)

        PJGUIPages.setCurrentIndex(1)
        self.node_content_tabs.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(PJGUIPages)
    # setupUi

    def retranslateUi(self, PJGUIPages):
        PJGUIPages.setWindowTitle(QCoreApplication.translate("PJGUIPages", u"StackedWidget", None))
        self.filename_prefix_label.setText(QCoreApplication.translate("PJGUIPages", u"File name prefix:", None))
        self.tree_dag_node_label.setText(QCoreApplication.translate("PJGUIPages", u"Tree DAG node(s)", None))
        self.random_seed_prefix_label.setText(QCoreApplication.translate("PJGUIPages", u"Random seed:", None))
        self.stoch_map_label.setText(QCoreApplication.translate("PJGUIPages", u"Stochastic mapping", None))
        self.settings_label.setText(QCoreApplication.translate("PJGUIPages", u"Output configuration", None))
        self.sim_config_label.setText(QCoreApplication.translate("PJGUIPages", u"Simulation configuration", None))
        self.attr_name_label.setText(QCoreApplication.translate("PJGUIPages", u"Attribute name", None))
        self.node_content_tabs.setTabText(self.node_content_tabs.indexOf(self.values_tab), QCoreApplication.translate("PJGUIPages", u"Value(s)", None))
        self.node_content_tabs.setTabText(self.node_content_tabs.indexOf(self.summary_tab), QCoreApplication.translate("PJGUIPages", u"Summary stats.", None))
        self.one_sample_radio.setText(QCoreApplication.translate("PJGUIPages", u"One sample", None))
        self.all_samples_radio.setText(QCoreApplication.translate("PJGUIPages", u"All samples", None))
        self.reconstructed_tree_check.setText(QCoreApplication.translate("PJGUIPages", u"Reconstructed", None))
        self.sample_idx_spin.setPrefix(QCoreApplication.translate("PJGUIPages", u"Sample #", None))
        self.repl_idx_spin.setPrefix(QCoreApplication.translate("PJGUIPages", u"Replicate #", None))
        self.save_pgm_node_plot.setText(QCoreApplication.translate("PJGUIPages", u"Save plot as", None))
        self.model_label.setText(QCoreApplication.translate("PJGUIPages", u"Model nodes", None))
        self.redraw_node.setText(QCoreApplication.translate("PJGUIPages", u"Redraw", None))
        self.clear_model.setText(QCoreApplication.translate("PJGUIPages", u" Clear model", None))
        self.cmd_prompt_label.setText(QCoreApplication.translate("PJGUIPages", u"Command prompt", None))
#if QT_CONFIG(tooltip)
        self.cmd_prompt.setToolTip("")
#endif // QT_CONFIG(tooltip)
        self.cmd_prompt.setPlaceholderText("")
        self.compare_node_label.setText(QCoreApplication.translate("PJGUIPages", u"Node to compare", None))
        self.avg_replicate_check_button.setText(QCoreApplication.translate("PJGUIPages", u"Replicates", None))
        self.compare_stats_label.setText(QCoreApplication.translate("PJGUIPages", u"Summary statistics", None))
        self.compare_csv_button.setText(QCoreApplication.translate("PJGUIPages", u"Compare to .csv (...)", None))
        self.draw_violins_button.setText(QCoreApplication.translate("PJGUIPages", u"Draw", None))
        self.save_violins_button.setText(QCoreApplication.translate("PJGUIPages", u"Save plot as", None))
        self.read_hpd_csv_button.setText(QCoreApplication.translate("PJGUIPages", u"Read HPDs .csv (...)", None))
        self.read_logfile_button.setText(QCoreApplication.translate("PJGUIPages", u"Directory to .log's (...)", None))
        self.node_name_label.setText(QCoreApplication.translate("PJGUIPages", u"Parameter name in .log", None))
        self.coverage_label.setText(QCoreApplication.translate("PJGUIPages", u"Coverage", None))
        self.draw_cov_button.setText(QCoreApplication.translate("PJGUIPages", u"Draw", None))
        self.save_cov_button.setText(QCoreApplication.translate("PJGUIPages", u"Save plot as", None))
        self.coverage_node_label.setText(QCoreApplication.translate("PJGUIPages", u"Non-det. nodes", None))
        self.compare_stats_label_2.setText(QCoreApplication.translate("PJGUIPages", u"Summary statistics", None))
    # retranslateUi
