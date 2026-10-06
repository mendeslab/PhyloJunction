import sys
from contextlib import ExitStack
import os
import re
import typing as ty
import numpy as np
import pandas as pd
import arviz as az  # type: ignore
from natsort import natsorted  # type: ignore
from tabulate import tabulate  # type: ignore
from PySide6.QtWidgets import \
    QApplication, QMainWindow, QPushButton, QFileDialog, QMessageBox, QStyle
from PySide6.QtGui import QAction
from PySide6.QtCore import QPropertyAnimation, QEasingCurve, QTimer, QSignalBlocker, Qt

# pj imports #
from phylojunction.interface.pysidegui.content_main_window \
    import ContentGUIMainWindow
from phylojunction.pgm.pgm import DirectedAcyclicGraph
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.plotting.pj_organize as pjorg
import phylojunction.plotting.pj_draw as pjdraw
import phylojunction.readwrite.pj_read as pjread
import phylojunction.readwrite.pj_write as pjwrite
import phylojunction.utility.helper_functions as pjh
import phylojunction.utility.exception_classes as ec
import phylojunction.data.tree as pjdt


__author__ = "Fabio K. Mendes"
__email__ = "fmendes@lsu.edu"


class GUIModeling():

    dag_obj: DirectedAcyclicGraph
    cmd_log_list: ty.List[str]

    def __init__(self):
        self.dag_obj = DirectedAcyclicGraph()
        self.cmd_log_list = []
        # None marks an incomplete or unverified source; history also contains reset notes.
        self.active_script_lines: ty.Optional[ty.List[str]] = []

    def update_dag_random_seed(self, a_random_seed: int):
        self.dag_obj.random_seed = a_random_seed

    def parse_cmd_update_pgm(
            self,
            cmd_line_list,
            gui_main_window_obj,
            clear_cmd_log_list: bool = False):
        
        # in case we read two scripts
        # consecutively, we want to clear
        # the comand line history list
        if clear_cmd_log_list:
            self.cmd_log_list = []

        print("\nReading following command lines:")
        execution_failed = False

        for line in cmd_line_list:
            valid_cmd_line = None
            # removing whitespaces from left and right
            line = line.strip()
            print("  " + line)

            # Parse only after setting an explicitly requested initial seed.
            try:
                # will set random seed once if DAG does not already have one
                if not isinstance(self.dag_obj.random_seed, int):
                    a_random_seed: str = \
                        gui_main_window_obj.ui.ui_pages. \
                            random_seed_prefix_textbox.toPlainText().strip()
                
                    # random seed is not None and not empty string
                    if a_random_seed:
                        random_seed: int = int(a_random_seed)
                        self.dag_obj.random_seed = random_seed

                # cmdline2dag updates the current DAG in place.
                valid_cmd_line = cmdp.cmdline2dag(self.dag_obj, line)

            except Exception as e:
                execution_failed = True
                self.active_script_lines = None
                gui_main_window_obj.ui.bottom_label_left.setText("Warning produced")

                # if not event == "Simulate":
                if e.__context__:
                    # __context__ catches innermost exception
                    # update warnings page later, with e.__context__
                    # print(e.__context__)
                    # gui_main_window_obj.ui.ui_pages.warnings_textbox.setText(str(e.__context__))
                    print(e)
                    gui_main_window_obj.ui.ui_pages.warnings_textbox.setText(str(e))

                else:
                    print(e)  # update warnings page with e
                    gui_main_window_obj.ui.ui_pages.warnings_textbox.setText(str(e))

            if valid_cmd_line:
                # try:
                if len(self.cmd_log_list) >= 1 and self.cmd_log_list[-1] in ("\n#"):
                    self.cmd_log_list.append("")  # adding new line

                # except Exception as e:
                #     print("PhyloJunction's GUI: An error occurred: ", type(e).__name__, " - ", e)
                #     pass

                self.cmd_log_list.append(valid_cmd_line)
                if self.active_script_lines is not None:
                    self.active_script_lines.append(valid_cmd_line)

        gui_main_window_obj.refresh_resample_button(
            "Execution failed. Clear the model and enter commands, or load a script, to enable resampling."
            if execution_failed else None)

    def cmd_log(self):
        return "\n".join(self.cmd_log_list) + "\n\n"

    def clear(self):
        self.dag_obj = DirectedAcyclicGraph()
        self.active_script_lines = []


class GUIMainWindow(QMainWindow):

    input_script_filepath: str
    gui_modeling: GUIModeling
    is_avg_repl_check: bool

    pj_comparison_df: pd.DataFrame
    other_comparison_df: pd.DataFrame
    coverage_df: pd.DataFrame
    hpd_df: pd.DataFrame

    def __init__(self):
        super().__init__()

        self.setWindowTitle("PhyloJunction")

        self.init_top_menus()

        # setup main window #
        self.ui = ContentGUIMainWindow()

        # passing myself as argument to auxiliary setup class
        self.ui.setup_ui(self)

        #######################
        # Members for storage #
        # and operations      #
        #######################

        self.input_script_filepath = ""

        # model-related members #
        self.gui_modeling = GUIModeling()

        self.is_avg_repl_check = False

        self.pj_comparison_df = pd.DataFrame()
        self.other_comparison_df = pd.DataFrame()
        self.coverage_df = pd.DataFrame()
        self.hpd_df = pd.DataFrame()

        #########################
        # Below we specify what #
        # happens when buttons  #
        # are pressed           #
        #########################

        # toggle button #
        self.ui.menu_button.clicked.connect(self.toggle_button_animation)

        # menu button #
        self.ui.pgm_button.clicked.connect(self.show_pgm_page)

        # pgm page button #
        self.ui.cmd_log_button.clicked.connect(self.show_cmd_log_page)

        # compare page button #
        self.ui.compare_button.clicked.connect(self.show_compare_page)

        # coverage page button #
        self.ui.covg_button.clicked.connect(self.show_coverage_page)

        # settings page button #
        self.ui.settings_button.clicked.connect(self.show_settings_page)

        # warnings page button #
        self.ui.warning_button.clicked.connect(self.show_warnings_page)

        ###################################
        # gui_pages buttons               #
        # (created by QtCreator/Designer) #
        ###################################

        ############
        # PGM page #
        ############

        # cmd line enter #
        self.ui.ui_pages.cmd_prompt.returnPressed.connect(
            self.parse_cmd_update_gui)

        # Selection changes include keyboard navigation as well as mouse clicks.
        self.ui.ui_pages.node_list.currentItemChanged.connect(
            lambda current, previous: self.do_selected_node_dag_page())
        self.ui.ui_pages.reconstructed_tree_check.toggled.connect(
            self.refresh_selected_node_display_plot_spin)

        # radio button update #
        self.ui.ui_pages.one_sample_radio.clicked.connect(
            self.one_sample_clicked)
        self.ui.ui_pages.all_samples_radio.clicked.connect(
            self.all_samples_clicked)

        # spin boxes update #
        self.ui.ui_pages.sample_idx_spin.valueChanged.connect(
            self.refresh_selected_node_display_plot_spin)
        self.ui.ui_pages.repl_idx_spin.valueChanged.connect(
            self.refresh_selected_node_display_plot_spin)

        # save plot #
        self.ui.ui_pages.save_pgm_node_plot.clicked.connect(
            lambda plot2fig: self.write_plot_to_file(
                self.ui.ui_pages.pgm_page_matplotlib_widget.fig,
                self.ui.ui_pages.pgm_page_matplotlib_widget.axes
            )
        )

        # Replay the specification to generate a new model, then display it.
        self.ui.ui_pages.resample_model.clicked.connect(self.resample_current_script)
        self.ui.ui_pages.resample_model.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_BrowserReload))
        self.refresh_resample_button()

        # clear model button #
        self.ui.ui_pages.clear_model.clicked.connect(
            lambda clear_model:
                self.clean_disable_everything(user_reset=True))

        ################
        # Compare page #
        ################
        self.ui.ui_pages.compare_csv_button.clicked.connect(self.read_compare_csv)
        self.ui.ui_pages.compare_node_list.itemClicked.connect(self.do_selected_node_compare_page)
        self.ui.ui_pages.avg_replicate_check_button.clicked.connect(self.avg_repl_check)
        self.ui.ui_pages.draw_violins_button.clicked.connect(self.draw_violin)
        self.ui.ui_pages.save_violins_button.clicked.connect(
            lambda plot2fig: self.write_plot_to_file(
                self.ui.ui_pages.compare_page_matplotlib_widget.fig,
                self.ui.ui_pages.compare_page_matplotlib_widget.axes
            )
        )

        #################
        # Coverage page #
        #################
        self.ui.ui_pages.read_hpd_csv_button.clicked.connect(self.read_coverage_hpd_csv)
        self.ui.ui_pages.read_logfile_button.clicked.connect(self.read_log_files_dir_build_hpd)
        self.ui.ui_pages.coverage_node_list.itemClicked.connect(self.do_selected_node_coverage_page)

        # side-effect: sets xticklabels in
        self.ui.ui_pages.draw_cov_button.clicked.connect(self.draw_cov)
        self.ui.ui_pages.save_cov_button.clicked.connect(
            lambda plot2fig: self.write_plot_to_file(
                self.ui.ui_pages.coverage_page_matplotlib_widget.fig,
                self.ui.ui_pages.coverage_page_matplotlib_widget.axes
            )
        )

        # show #
        self.show()

    # set up all top menus
    def init_top_menus(self):
        # setup menu bar #
        menubar = self.menuBar()

        # first item (pj menu)
        pj_menu = menubar.addMenu("PhyloJunction")

        # pj menu items
        # certain keywords cannot be used in
        # menubars on Mac, like "quit", "about"
        # https://stackoverflow.com/questions/73293306/string-containing-about-cannot-be-added-to-the-menu-in-pyside6
        about_action = pj_menu.addAction("About", self.print_about)
        about_action.setMenuRole(QAction.MenuRole.NoRole)
        pj_menu.addAction(about_action)

        quit_action = pj_menu.addAction("Quit", self.quit_app_button_clicked)
        quit_action.setMenuRole(QAction.MenuRole.NoRole)
        pj_menu.addAction(quit_action)

        # second item (file menu)
        file_menu = menubar.addMenu("File")

        # file menu items
        file_menu.addAction("Read script", self.read_execute_script)

        file_menu.addAction("Load model", self.load_model)
        file_menu.addSeparator()
        file_menu.addAction("Save data as", self.write_data_to_dir)
        file_menu.addAction("Save stochastic mappings as", self.write_smaps_to_dir)
        file_menu.addAction("Save model as", self.write_model_to_file)

    # verify which menu buttons are selected #
    def reset_selection(self):
        for btn in self.ui.left_menu.findChildren(QPushButton):
            try:
                btn.set_active(False)

            except Exception as e:
                print("An error occurred: ", type(e).__name__, " - ", e)
                pass

    def toggle_button_animation(self):
        # get menu width #
        menu_width = self.ui.left_menu.width()

        # check width #
        width = 50  # base width
        if menu_width == 50:
            width = 180

        # start animation #
        self.animation = QPropertyAnimation(self.ui.left_menu, b"minimumWidth")
        self.animation.setStartValue(menu_width)
        self.animation.setEndValue(width)
        self.animation.setDuration(500)  # 500 ms
        self.animation.setEasingCurve(QEasingCurve.InOutCirc)
        self.animation.start()

    #########################################
    # Functions for activating page buttons #
    #########################################

    def show_pgm_page(self):
        self.reset_selection()
        # ui_pages is a member of self.ui,
        # which is a ContentGUIMainWindow instance
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.pgm_page)
        self.ui.pgm_button.set_active(True)
        self.ui.top_label_left.setText("MODEL SPECIFICATION")

    def show_compare_page(self):
        self.reset_selection()
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.compare_page)
        self.ui.compare_button.set_active(True)
        self.ui.top_label_left.setText("IMPLEMENTATION COMPARISON")

    def show_coverage_page(self):
        self.reset_selection()
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.coverage_page)
        self.ui.covg_button.set_active(True)
        self.ui.top_label_left.setText("COVERAGE VALIDATION")

    def show_cmd_log_page(self):
        self.reset_selection()
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.cmd_log_page)
        self.ui.cmd_log_button.set_active(True)
        self.ui.top_label_left.setText("COMMAND LOG")

    def show_settings_page(self):
        self.reset_selection()
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.settings_page)
        self.ui.warning_button.set_active(True)
        self.ui.top_label_left.setText("SETTINGS")

    def show_warnings_page(self):
        self.reset_selection()
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.warnings_page)
        self.ui.warning_button.set_active(True)
        self.ui.top_label_left.setText("WARNINGS")
        self.ui.bottom_label_left.setText("")

    ###################################
    # Various event-related functions #
    ###################################

    def parse_cmd_update_gui(self):
        cmd_line = self.ui.ui_pages.cmd_prompt.text()  # grab input

        # in case the user enters many lines
        # (also ignores empty lines)
        cmd_line_list = [cl for cl in re.split("\n", cmd_line) if cl]

        # read commands
        # (side-effect: gui_modeling stores cmd hist)
        self.gui_modeling.parse_cmd_update_pgm(cmd_line_list, self)

        # update GUI cmd history
        self.refresh_cmd_history()

        self.ui.ui_pages.cmd_prompt.clear()  # clear

        # updating node lists
        self.refresh_node_lists()

    def selected_node_display(
            self,
            node_dag,
            do_all_samples,
            sample_idx=None,
            repl_idx=0,
            repl_size=1,
            display_reconstructed=False):

        display_node_dag_value_str = str()
        display_node_dag_stat_str = str()

        is_tree = False
        # try:
        if isinstance(node_dag.value, list) and \
                isinstance(node_dag.value[0], pjdt.AnnotatedTree):
            is_tree = True

        # if deterministic, not subscriptable
        # except Exception as e:
            # print("(pj_gui.py) An error occurred: ", type(e).__name__, " - ", e)
            # pass

        # first we do values
        # we care about a specific sample and maybe a specific replicate
        if sample_idx is not None and not do_all_samples:
            start = sample_idx * repl_size
            end = start + repl_size

            # values
            display_node_dag_value_str = \
                node_dag.get_start2end_str(
                    start,
                    end,
                    repl_idx=repl_idx,
                    is_tree=is_tree,
                    display_reconstructed=display_reconstructed)

            # summary stats
            display_node_dag_stat_str = \
                node_dag.get_node_stats_str(
                    start,
                    end,
                    repl_idx)

        # we get all samples
        else:
            # just calling __str__
            display_node_dag_value_str = \
                self.gui_modeling.dag_obj \
                .get_display_str_by_name(node_dag.node_name)

            # getting all values
            display_node_dag_stat_str = \
                node_dag.get_node_stats_str(
                    0, len(node_dag.value), repl_idx)  # summary stats

        # print("Set values_content QLineEdit widget with text: " + display_node_dag_value_str)
        # print("Set summary_content QLineEdit widget with text: " + display_node_dag_stat_str)
        self.ui.ui_pages.values_content.setText(display_node_dag_value_str)
        self.ui.ui_pages.summary_content.setText(display_node_dag_stat_str)

    def selected_node_plot(
            self,
            fig_obj,
            fig_axes,
            node_dag,
            do_all_samples,
            sample_idx=None,
            repl_idx=0,
            repl_size=1,
            draw_reconstructed=False):
        """
        Plot DAG node on 'node_display_fig_axes' (Axes object) scoped to 'call_gui()',
        then update canvas with new plot
        """

        if isinstance(node_dag.value, list):
            # if a tree
            if isinstance(node_dag.value[0], pjdt.AnnotatedTree):
                self.draw_node_dag(fig_axes,
                                   node_dag,
                                   sample_idx=sample_idx,
                                   repl_idx=repl_idx,
                                   repl_size=repl_size,
                                   draw_reconstructed=draw_reconstructed)

            # when not a tree
            else:
                if do_all_samples:
                    self.draw_node_dag(fig_axes,
                                       node_dag,
                                       repl_size=repl_size)

                else:
                    self.draw_node_dag(fig_axes,
                                        node_dag,
                                        sample_idx=sample_idx,
                                        repl_size=repl_size)

        # if deterministic, not subscriptable
        # except Exception as e:
        #     print("An error occurred: ", type(e).__name__, " - ", e)

        # deterministic node
        else:
            self.draw_node_dag(fig_axes, node_dag)

        fig_obj.canvas.draw()

    def selected_node_read(self, node_name: str):
        node_dag = self.gui_modeling.dag_obj.get_node_dag_by_name(node_name)
        # this is n_sim inside sampling distribution classes
        sample_size = len(node_dag)
        repl_size = node_dag.repl_size

        return node_dag, sample_size, repl_size

    def do_selected_node_dag_page(self, spin_buttons_clicked: bool = False):
        """
        Display selected node's string representation and
        plot it on canvas if possible, for pgm page
        (model specification menu)
        """

        # first get selected node's name #
        selected_node_name: str = ""
        active_item = self.ui.ui_pages.node_list.currentItem()
        if active_item:
            selected_node_name = active_item.text()

        if selected_node_name not in ("", None):
            # reading node information #
            node_dag, sample_size, repl_size = \
                self.selected_node_read(selected_node_name)

            # Range and mode initialization must not trigger intermediate plots.
            with QSignalBlocker(self.ui.ui_pages.sample_idx_spin), \
                    QSignalBlocker(self.ui.ui_pages.repl_idx_spin):
                self.ui.ui_pages.repl_idx_spin.setMaximum(repl_size)
                self.ui.ui_pages.sample_idx_spin.setMaximum(sample_size)
                if not spin_buttons_clicked:
                    self.init_and_refresh_radio_spin(node_dag, sample_size, repl_size)

            fig_obj = self.ui.ui_pages.pgm_page_matplotlib_widget.fig
            fig_axes = self.ui.ui_pages.pgm_page_matplotlib_widget.axes

            do_all_samples = \
                self.ui.ui_pages.all_samples_radio.isChecked()

            do_reconstructed = \
                self.ui.ui_pages.reconstructed_tree_check.isChecked()

            #################################
            # Collect information from spin #
            #################################

            sample_idx = self.ui.ui_pages.sample_idx_spin.value() - 1
            repl_idx = self.ui.ui_pages.repl_idx_spin.value() - 1

            ###############
            # Now do node #
            ###############

            self.selected_node_plot(fig_obj,
                                    fig_axes,
                                    node_dag,
                                    do_all_samples,
                                    sample_idx=sample_idx,
                                    repl_idx=repl_idx,
                                    repl_size=repl_size,
                                    draw_reconstructed=do_reconstructed)

            self.selected_node_display(node_dag,
                                       do_all_samples,
                                       sample_idx=sample_idx,
                                       repl_idx=repl_idx,
                                       repl_size=repl_size,
                                       display_reconstructed=do_reconstructed)

        # no nodes to do
        else:
            pass

    def do_selected_node_compare_page(self):
        """
        Chose node to get summary statistics
        and compare against user implementation
        (implementation comparison menu)
        """

        # first get selected node's name #
        selected_node_name: str = ""
        active_item = self.ui.ui_pages.compare_node_list.currentItem()
        if active_item:
            selected_node_name = active_item.text()

        # if there is at least one
        # sampled node with a name
        if selected_node_name not in ("", None):

            # reading node information #
            selected_node_dag, selected_node_sample_size, selected_node_repl_size = \
                self.selected_node_read(selected_node_name)

            # could be more efficient, but this
            # makes sure that stashes are always up-to-date
            scalar_output_stash, tree_output_stash = \
                pjwrite.prep_data_df(self.gui_modeling.dag_obj)

            # collecting scalars #
            _, scalar_value_df_dict, scalar_repl_summary_df = \
                scalar_output_stash

            # collecting trees #
            tree_value_df_dict, tree_ann_value_df_dict, tree_rec_value_df_dict, \
                tree_rec_ann_value_df_dict, tree_summary_df_dict, tree_repl_summary_df_dict, \
                tree_living_nd_states_str_dict, tree_living_nd_states_str_nexus_dict, \
                tree_internal_nd_states_str_dict = tree_output_stash

            # adding "program" column to
            # PJ's pandas DataFrame's stored inside dicts
            for repl_size in scalar_value_df_dict:
                scalar_value_df_dict[selected_node_repl_size] \
                    .loc[:, "program"] = "PJ"

            for tr_nd_name in tree_summary_df_dict:
                tree_summary_df_dict[tr_nd_name] \
                    .loc[:, "program"] = "PJ"

            for tr_w_repl_nd_name in tree_repl_summary_df_dict:
                tree_repl_summary_df_dict[tr_w_repl_nd_name] \
                    .loc[:, "program"] = "PJ"

            # scalar was selected #
            if isinstance(selected_node_dag.value[0], (int, float, np.float64)):
                self.ui.ui_pages.summary_stats_list.clear()

                if not self.is_avg_repl_check:
                    self.pj_comparison_df = scalar_value_df_dict[selected_node_repl_size]

                elif selected_node_repl_size > 1:
                    self.ui.ui_pages.summary_stats_list.addItems(
                        ["average", "std. dev."]
                    )
                    self.pj_comparison_df = scalar_repl_summary_df

            # tree was selected #
            elif isinstance(selected_node_dag.value[0], pjdt.AnnotatedTree):
                self.pj_comparison_df = tree_summary_df_dict[selected_node_name]
                self.ui.ui_pages.summary_stats_list.clear()
                self.ui.ui_pages.summary_stats_list.addItems(
                    [tree_stat for tree_stat in self.pj_comparison_df.keys()
                     if tree_stat not in ("program", "sample", "replicate")]
                )

        # no sampled nodes
        else:
            print("If there are no sampled nodes in PJ model, we cannot compare against PJ. Later write an Exception for this")
            pass

    def do_selected_node_coverage_page(self):
        """
        Chose node to get values or summary
        statistics for coverage plot
        (coverage menu)
        """

        # first get selected node's name #
        selected_node_name: str = ""
        active_item = self.ui.ui_pages.coverage_node_list.currentItem()
        if active_item:
            selected_node_name = active_item.text()

        # if there is at least one
        # sampled node with a name
        if selected_node_name not in ("", None):

            # reading node information #
            selected_node_dag, selected_node_sample_size, selected_node_repl_size = \
                self.selected_node_read(selected_node_name)

            # could be more efficient, but this
            # makes sure that stashes are always up-to-date
            scalar_output_stash, tree_output_stash = \
                pjwrite.prep_data_df(self.gui_modeling.dag_obj)

            # collecting scalars #
            scalar_constant_value_df, scalar_value_df_dict, scalar_repl_summary_df = \
                scalar_output_stash

            # collecting trees #
            tree_value_df_dict, tree_ann_value_df_dict, tree_rec_value_df_dict, \
                tree_rec_ann_value_df_dict, tree_summary_df_dict, tree_repl_summary_df_dict, \
                tree_living_nd_states_str_dict, tree_living_nd_states_str_nexus_dict, \
                tree_internal_nd_states_str_dict = tree_output_stash

            # str because could be constant set by hand
            if isinstance(selected_node_dag.value[0], (str, int, float, np.float64)):
                if selected_node_repl_size <= 1:
                    self.ui.ui_pages.cov_summary_stats_list.clear()

                    # sampled, non-deterministic
                    if selected_node_dag.is_sampled:
                        self.coverage_df = scalar_value_df_dict[selected_node_repl_size]

                    # constant, non-deterministic
                    else:
                        self.coverage_df = scalar_constant_value_df

                else:
                    self.ui.ui_pages.cov_summary_stats_list.addItems(
                        ["average", "std. dev."]
                    )
                    self.coverage_df = scalar_repl_summary_df

        else:
            print("If there are no non-deterministic nodes in PJ model, no sensical coverage can be calculated. Later write an Exception for this")
            pass

    ##################
    # Functions for  #
    # reading and    #
    # loading events #
    ##################

    def read_execute_script(self):

        # read file path #
        script_fp, filter = \
            QFileDialog.getOpenFileName(
                parent=self,
                caption="Read script",
                dir=".",
                filter="*.pj")

        if script_fp:
            # reset everything #
            self.clean_disable_everything()

            # set seed
            a_random_seed = self.ui.ui_pages.random_seed_prefix_textbox.toPlainText().strip()
            # is not None and is not an empty string
            if a_random_seed:
                random_seed: int = int(a_random_seed)
                self.gui_modeling.update_dag_random_seed(random_seed)

            # read all command line strings
            cmd_line_list = pjread.read_text_file(script_fp)

            # (side-effect: gui_modeling stores cmd hist) #
            self.gui_modeling.parse_cmd_update_pgm(cmd_line_list, self, clear_cmd_log_list=True)

            # update GUI cmd history #
            self.refresh_cmd_history(user_reset=True)

            # add node names to...
            # pgm page node list
            self.ui.ui_pages.node_list.addItems(
                [node_name for node_name in
                 self.gui_modeling.dag_obj.name_node_dict]
            )

            # compare page node list
            self.ui.ui_pages.compare_node_list.addItems(
                [nd.node_name for nd in
                 self.gui_modeling.dag_obj.get_sorted_node_dag_list() if
                 nd.is_sampled]
            )

            # coverage page node list
            self.ui.ui_pages.coverage_node_list.addItems(
                [nd.node_name for nd in
                 self.gui_modeling.dag_obj.get_sorted_node_dag_list() if
                 not nd.is_deterministic]
            )

    # Availability depends on executable source, not on the seed: a seeded click explains
    # why it cannot resample. Preserve a more specific unavailable reason when supplied.
    def refresh_resample_button(self, unavailable_reason=None):
        button = self.ui.ui_pages.resample_model
        lines = self.gui_modeling.active_script_lines
        button.setEnabled(bool(lines))
        if lines:
            button.setToolTip("Rerun the current script to generate new samples.")
        elif lines is not None:
            button.setToolTip("Load a script or enter commands to enable resampling.")
        elif unavailable_reason:
            button.setToolTip(unavailable_reason)

    # Evaluate into a separate DAG before replacing the displayed model. Retain source and
    # view preferences across the existing cleanup routine, which resets GUIModeling.
    def resample_current_script(self):
        u = self.ui.ui_pages
        if u.random_seed_prefix_textbox.toPlainText().strip():
            QMessageBox.warning(self, "Resampling disabled",
                                "Resampling is disabled while a random seed is specified. "
                                "Clear the Random seed field in Settings to generate new samples.")
            return
        if not self.gui_modeling.active_script_lines or not u.resample_model.isEnabled():
            return

        lines = list(self.gui_modeling.active_script_lines)
        history = list(self.gui_modeling.cmd_log_list)
        item = u.node_list.currentItem()
        selected_name = item.text() if item else None
        sample, replicate = u.sample_idx_spin.value(), u.repl_idx_spin.value()
        all_samples = u.all_samples_radio.isChecked()
        reconstructed = u.reconstructed_tree_check.isChecked()
        u.resample_model.setEnabled(False)
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            try:
                candidate = cmdp.script2dag("\n".join(lines), in_pj_file=False, random_seed=None)
            except Exception as error:
                self.ui.bottom_label_left.setText("Resampling failed; previous results retained")
                u.warnings_textbox.setText(str(error))
                return

            # No callbacks may draw an intermediate selection or out-of-range sample.
            with ExitStack() as blocked:
                for widget in (u.node_list, u.sample_idx_spin, u.repl_idx_spin,
                               u.one_sample_radio, u.all_samples_radio, u.reconstructed_tree_check):
                    blocked.enter_context(QSignalBlocker(widget))
                self.clean_disable_everything()
                u.resample_model.setEnabled(False)
                self.gui_modeling.dag_obj = candidate
                self.gui_modeling.active_script_lines = lines
                self.gui_modeling.cmd_log_list = history
                self.refresh_cmd_history()
                self.refresh_node_lists()
                matches = u.node_list.findItems(selected_name, Qt.MatchFlag.MatchExactly) if selected_name else []
                if matches:
                    u.node_list.setCurrentItem(matches[0])
                elif u.node_list.count():
                    u.node_list.setCurrentRow(0)
                if u.node_list.currentItem():
                    node, sample_size, repl_size = self.selected_node_read(u.node_list.currentItem().text())
                    u.sample_idx_spin.setMaximum(sample_size)
                    u.repl_idx_spin.setMaximum(repl_size)
                    self.init_and_refresh_radio_spin(node, sample_size, repl_size)
                    if u.all_samples_radio.isEnabled():
                        if all_samples or not u.one_sample_radio.isEnabled():
                            u.all_samples_radio.setChecked(True)
                            u.sample_idx_spin.setMinimum(0)
                            u.sample_idx_spin.setValue(0)
                            u.sample_idx_spin.setEnabled(False)
                        else:
                            u.one_sample_radio.setChecked(True)
                            u.sample_idx_spin.setMinimum(1)
                            u.sample_idx_spin.setEnabled(True)
                    if u.sample_idx_spin.isEnabled():
                        u.sample_idx_spin.setValue(sample)
                    if u.repl_idx_spin.isEnabled():
                        u.repl_idx_spin.setValue(replicate)
                    u.reconstructed_tree_check.setChecked(reconstructed)
            try:
                self.do_selected_node_dag_page(spin_buttons_clicked=True)
                self.ui.bottom_label_left.setText("Resampled")
            except Exception as error:
                self.ui.bottom_label_left.setText("Resampled, but displaying the new model failed")
                u.warnings_textbox.setText(str(error))
        finally:
            QApplication.restoreOverrideCursor()
            self.refresh_resample_button()

    def read_compare_csv(self):

        # read file path #
        compare_csv_fp, filter = QFileDialog.getOpenFileName(parent=self, caption="Read compare .csv", dir=".", filter="*.csv")

        if compare_csv_fp:
            if not os.path.isfile(compare_csv_fp):
                print("Could not find " + compare_csv_fp)
                # TODO: raise specific Error later

            try:
                self.other_comparison_df = \
                    pjread.read_csv_into_dataframe(compare_csv_fp)

                other_comparison_df_str = \
                    tabulate(
                        self.other_comparison_df,
                        self.other_comparison_df.head(),
                        tablefmt="plain",
                        showindex=False).lstrip()

                self.ui.ui_pages.compare_csv_textbox.setText(other_comparison_df_str)

            except Exception as e:
                print("An error occurred: ", type(e).__name__, " - ", e)
                print(("Could not load .csv file into pandas DataFrame. "
                       "Later write an Exception for this"))

        else:
            pass  # event canceled by user

    def read_coverage_hpd_csv(self):

        # read file path #
        coverage_hpd_csv_fp, filter = \
            QFileDialog.getOpenFileName(
                parent=self, caption="Read HPDs .csv", dir=".",
                filter="*.csv")

        if coverage_hpd_csv_fp:
            if not os.path.isfile(coverage_hpd_csv_fp):
                print("Could not find " + coverage_hpd_csv_fp)
                # TODO: raise specific Error later

            try:
                self.hpd_df = \
                    pjread.read_csv_tsv_into_dataframe(coverage_hpd_csv_fp)

                coverage_df_str = \
                    tabulate(
                        self.hpd_df,
                        self.hpd_df.head(),
                        tablefmt="plain",
                        showindex=False).lstrip()

                self.ui.ui_pages.coverage_csv_textbox.setText(coverage_df_str)

            except Exception as e:
                print("Could not load .csv file into pandas DataFrame. Later write an Exception for this")
                print("An error occurred: ", type(e).__name__)
                pass

        else:
            pass  # event canceled by user

    def read_log_files_dir_build_hpd(self):
        
        # read variable name #
        log_par_name = self.ui.ui_pages.log_node_name_lineedit.text()
        if not log_par_name:
            # try selected node's name #
            active_item = self.ui.ui_pages.coverage_node_list.currentItem()

            if active_item:
                log_par_name = active_item.text()
        
        # read directory path #
        logs_dir = \
            QFileDialog.getExistingDirectory(caption="Read .log directory",
                                             dir=".")
        
        if not logs_dir.endswith("/"):
            logs_dir += "/"

        if os.path.isdir(logs_dir):
            summary_dict = {"posterior_mean": [],
                            "lower_95hpd": [],
                            "higher_95hpd": []}
            log_files = natsorted(
                [f for f in os.listdir(logs_dir) if f.endswith(".log")])
            n_logs = len(log_files)

            if n_logs <= 1:
                raise ec.PJIOFileDoesNotExistError("read_log_files_dir_build_hpd",
                                                   ".log files inside " + logs_dir)
            
            for lf in log_files:
                lf_fp = logs_dir + lf
                log_df = pjread.read_csv_tsv_into_dataframe(lf_fp, is_file_csv=False)

                # check if log_par_name is in the .log file
                if log_par_name not in list(log_df.columns.values):
                    raise ec.MissingColumnName(log_par_name,
                                               "Could not parse .log file")

                mcmc_samples_np_arr = log_df[log_par_name].to_numpy()
                posterior_mean = mcmc_samples_np_arr.mean()
                lower_95hpd, higher_95hpd = az.hdi(mcmc_samples_np_arr, hdi_prob=.95)

                summary_dict["posterior_mean"].append(posterior_mean)
                summary_dict["lower_95hpd"].append(lower_95hpd)
                summary_dict["higher_95hpd"].append(higher_95hpd)
                
                # debugging
                # print(list(log_df.columns.values))
                # print(log_df[log_par_name])

            self.hpd_df = pd.DataFrame(summary_dict)
            # the rest is carried out by self.draw_cov()

        else:
            pass  # event canceled by user

    def load_model(self):
        # read file path #
        model_fp, filter = QFileDialog.getOpenFileName(
            parent=self,
            caption="Load model",
            dir=".",
            filter="*.pickle")

        if model_fp:
            self.clean_disable_everything()
            self.ui.ui_pages.cmd_log_textbox.clear()
            self.gui_modeling.dag_obj, self.gui_modeling.cmd_log_list = \
                pjread.read_serialized_pgm(model_fp)
            self.gui_modeling.active_script_lines = None
            self.refresh_resample_button("This saved model has no verified script. Load its script to resample.")
            self.refresh_node_lists()
            self.ui.ui_pages.cmd_log_textbox.setText(self.gui_modeling.cmd_log())

    ##################
    # Functions for  #
    # writing events #
    ##################

    def write_plot_to_file(self, fig_obj, axes_obj):
        is_plot = True

        if len(axes_obj.get_xticklabels()) == 0 \
                and len(axes_obj.get_yticklabels()) == 0:
            is_plot = False

        if is_plot:
            # get fp
            fig_fp, filter = \
                QFileDialog.getSaveFileName(self, "Save file", "", ".png")

            head, tail = os.path.split(fig_fp)
            prefix = self.ui.ui_pages.filename_prefix_textbox.toPlainText()

            # should user cancel, then no error!
            if head and tail:
                if prefix:
                    pjwrite.write_fig_to_file(
                        head + "/" + prefix + "_" + tail,
                        fig_obj
                    )

                else:
                    pjwrite.write_fig_to_file(
                        head + "/" + tail,
                        fig_obj
                    )

    def write_model_to_file(self, prefix: str = ""):
        # get file path
        pickle_fp, filter = QFileDialog.getSaveFileName(self, "Save file", "", "")

        if not prefix:
            prefix = self.ui.ui_pages.filename_prefix_textbox.toPlainText()

        # pickling and saving PGM
        pjwrite.dump_serialized_pgm(
            pickle_fp,
            self.gui_modeling.dag_obj,
            self.gui_modeling.cmd_log_list,
            prefix=prefix)

    def write_data_to_dir(self):
        # get dir path
        data_out_dir = QFileDialog.getExistingDirectory(
            self, "Save in folder", "./", QFileDialog.ShowDirsOnly
        )

        prefix = self.ui.ui_pages.filename_prefix_textbox.toPlainText()

        pjwrite.dump_pgm_data(
            data_out_dir,
            self.gui_modeling.dag_obj,
            prefix=prefix)

    def write_smaps_to_dir(self):
        # get dir path
        data_out_dir = QFileDialog.getExistingDirectory(
            self, "Save in folder", "./", QFileDialog.ShowDirsOnly
        )

        prefix = self.ui.ui_pages.filename_prefix_textbox.toPlainText()
        mapped_attr_name = \
            self.ui.ui_pages.attr_name_textbox.toPlainText()
        tr_dag_node_name_list = \
            self.ui.ui_pages.tr_dag_node_textbox.toPlainText()

        pjwrite.dump_trees_rb_smap_dfs(
            data_out_dir,
            self.gui_modeling.dag_obj,
            tr_dag_node_name_list,
            mapped_attr_name,
            prefix=prefix)

    #################
    # Functions for #
    # drawing       #
    #################

    def draw_node_dag(
            self,
            axes,
            node_dag,
            sample_idx=None,
            repl_idx=0,
            repl_size=1,
            draw_reconstructed=False) -> None:

        return node_dag.plot_node(
            axes,
            sample_idx=sample_idx,
            repl_idx=repl_idx,
            repl_size=repl_size,
            draw_reconstructed=draw_reconstructed)

    def draw_violin(self):
        thing_to_compare: str = ""
        compare_node_list_empty = self.ui.ui_pages.compare_node_list.count() == 0
        summary_stat_list_empty = self.ui.ui_pages.summary_stats_list.count() == 0

        # before draw, we reset violins
        self.ui.ui_pages.compare_page_matplotlib_widget.fig.clf()
        self.ui.ui_pages.compare_page_matplotlib_widget.initialize_axes()
        self.ui.ui_pages.compare_page_matplotlib_widget.fig.canvas.draw()

        # first get selected node's name #
        selected_node_name: str = ""
        active_item = self.ui.ui_pages.compare_node_list.currentItem()
        if active_item:
            selected_node_name = active_item.text()

        if selected_node_name not in ("", None):
            # reading node information #
            node_dag, sample_size, repl_size = \
                self.selected_node_read(selected_node_name)

            if not self.pj_comparison_df.empty and \
                    not self.other_comparison_df.empty:

                # scalar
                if isinstance(node_dag.value[0], (int, float, np.float64)):
                    if not self.is_avg_repl_check:
                        thing_to_compare = selected_node_name

                    # averaging over replicates
                    else:
                        if not summary_stat_list_empty:
                            if self.ui.ui_pages.summary_stats_list.currentItem():
                                thing_to_compare = self.ui.ui_pages.summary_stats_list.currentItem().text()

                            # summary stats were not picked, doing nothing
                            else:
                                pass

                # custom class object (e.g., AnnotatedTree)
                else:
                    if not summary_stat_list_empty:
                        if self.ui.ui_pages.summary_stats_list.currentItem():
                            thing_to_compare = self.ui.ui_pages.summary_stats_list.currentItem().text()

                        # summary stats were not picked, doing nothing
                        else:
                            pass

                # debugging
                # print("node_dag.is_sampled = " + str(node_dag.is_sampled))
                # print(tabulate(self.pj_comparison_df, self.pj_comparison_df.head(), tablefmt="pretty", showindex=False).lstrip())
                # print(tabulate(self.other_comparison_df, self.other_comparison_df.head(), tablefmt="pretty", showindex=False).lstrip())

                joint_dataframe = pjorg.join_dataframes(
                    self.pj_comparison_df,
                    self.other_comparison_df,
                    thing_to_compare,
                    summaries_avg_over_repl=False)

                # something went wrong,
                # we clear comparison figure
                if joint_dataframe.empty or not node_dag.is_sampled:
                    # node list should only contain
                    # sampled nodes already, but
                    # just being sure...
                    print("joint_dataframe was empty or node was constant")
                    pass

                else:
                    # initializing plot
                    self.ui.ui_pages.compare_page_matplotlib_widget.initialize_axes(disabled_yticks=False)
                    self.ui.ui_pages.compare_page_matplotlib_widget.fig.canvas.draw()
                    fig_obj = self.ui.ui_pages.compare_page_matplotlib_widget.fig
                    fig_axes = self.ui.ui_pages.compare_page_matplotlib_widget.axes

                    pjdraw.plot_violins(
                        fig_obj,
                        fig_axes,
                        joint_dataframe,
                        "program",
                        thing_to_compare,
                        xlab="Program",
                        ylab=thing_to_compare)

    def draw_cov(self):
        selected_node_name: str = ""
        thing_to_validate: str = ""
        coverage_node_list_empty = self.ui.ui_pages.coverage_node_list.count() == 0

        # before draw, we reset coverage plot
        self.ui.ui_pages.coverage_page_matplotlib_widget.fig.clf()
        self.ui.ui_pages.coverage_page_matplotlib_widget.initialize_axes()
        self.ui.ui_pages.coverage_page_matplotlib_widget.fig.canvas.draw()

        # first get selected node's name #
        selected_node_name: str = ""
        active_item = self.ui.ui_pages.coverage_node_list.currentItem()
        if active_item:
            selected_node_name = active_item.text()

        if selected_node_name not in ("", None):
            # reading node information #
            node_dag, sample_size, repl_size = \
                self.selected_node_read(selected_node_name)

            if not self.coverage_df.empty and \
                    not self.hpd_df.empty:

                # scalar
                if isinstance(node_dag.value[0], (str, int, float, np.float64)):
                    thing_to_validate = selected_node_name

                # debugging
                # print(tabulate(self.coverage_df, self.coverage_df.head(), tablefmt="plain", showindex=False).lstrip())
                # print(tabulate(self.hpd_df, self.hpd_df.head(), tablefmt="plain", showindex=False).lstrip())

                full_cov_df = pd.concat([self.coverage_df, self.hpd_df], axis=1)
                full_cov_df = pjorg.add_within_hpd_col(full_cov_df, thing_to_validate)

                # debugging
                # print(tabulate(full_cov_df, full_cov_df.head(), tablefmt="plain", showindex=False).lstrip())

                covg = pjh.get_covg(full_cov_df, selected_node_name)

                # something went wrong, we clear validation figure
                if full_cov_df.empty:
                    pass

                else:
                    # initializing plot
                    self.ui.ui_pages.coverage_page_matplotlib_widget.initialize_axes(
                        disabled_yticks=False, disabled_xticks=False)
                    self.ui.ui_pages.coverage_page_matplotlib_widget.fig.canvas.draw()
                    fig_obj = self.ui.ui_pages.coverage_page_matplotlib_widget.fig
                    fig_axes = self.ui.ui_pages.coverage_page_matplotlib_widget.axes

                    # print bars!
                    pjdraw.plot_intervals(fig_obj,
                                          fig_axes,
                                          full_cov_df,
                                          thing_to_validate,
                                          "posterior_mean",
                                          ylab="Posterior mean")
                    
                    # print coverage to top-right window if not already there
                    covg_str = selected_node_name + "\t" + str(covg)
                    qlist_widget = self.ui.ui_pages.coverage_textbox
                    item_str_list = [qlist_widget.item(i).text() \
                                     for i in range(qlist_widget.count())]

                    if covg_str not in item_str_list:
                        qlist_widget.addItem(selected_node_name + "\t" + str(covg))

    #####################
    # Events related to #
    # initializing,     #
    # refreshing,       #
    # cleaning,         #
    # checking          #
    #####################

    def init_and_refresh_radio_spin(self, node_dag, sample_size, repl_size):

        def _prepare_for_tree(potential_repl: bool = False):
            # radio #
            # we always look one tree at a time
            self.ui.ui_pages.all_samples_radio.setChecked(False)
            self.ui.ui_pages.all_samples_radio.setCheckable(False)
            self.ui.ui_pages.all_samples_radio.setDisabled(True)

            self.ui.ui_pages.one_sample_radio.setEnabled(True)
            self.ui.ui_pages.one_sample_radio.setCheckable(True)
            self.ui.ui_pages.one_sample_radio.setChecked(True)

            # spin #
            self.ui.ui_pages.sample_idx_spin.setEnabled(True)
            self.ui.ui_pages.sample_idx_spin.setMinimum(1)
            self.ui.ui_pages.sample_idx_spin.setValue(1)

            self.ui.ui_pages.repl_idx_spin.setMinimum(1)
            self.ui.ui_pages.repl_idx_spin.setValue(1)
            if potential_repl:
                self.ui.ui_pages.repl_idx_spin.setEnabled(True)

        def _prepare_for_scalar(potential_repl: bool = False):
            # radio #
            self.ui.ui_pages.all_samples_radio.setEnabled(True)
            self.ui.ui_pages.all_samples_radio.setCheckable(True)
            self.ui.ui_pages.all_samples_radio.setChecked(True)

            # spin #
            self.ui.ui_pages.sample_idx_spin.setMinimum(1)
            self.ui.ui_pages.sample_idx_spin.setValue(1)
            self.ui.ui_pages.sample_idx_spin.setDisabled(True)

            # if there are replicates, then it makes sense to
            # even look at each sample individually
            if potential_repl:
                self.ui.ui_pages.one_sample_radio.setEnabled(True)
                self.ui.ui_pages.one_sample_radio.setCheckable(True)
                self.ui.ui_pages.one_sample_radio.setChecked(True)        

                # spin #
                self.ui.ui_pages.sample_idx_spin.setEnabled(True)

        def _nothing_to_spin_through():
            self.ui.ui_pages.one_sample_radio.setChecked(False)
            self.ui.ui_pages.one_sample_radio.setCheckable(False)
            self.ui.ui_pages.one_sample_radio.setDisabled(True)
            self.ui.ui_pages.all_samples_radio.setChecked(False)
            self.ui.ui_pages.all_samples_radio.setCheckable(False)
            self.ui.ui_pages.all_samples_radio.setDisabled(True)
            self.ui.ui_pages.sample_idx_spin.setMinimum(0)
            self.ui.ui_pages.sample_idx_spin.setValue(0)
            self.ui.ui_pages.sample_idx_spin.setDisabled(True)
            self.ui.ui_pages.repl_idx_spin.setMinimum(0)
            self.ui.ui_pages.repl_idx_spin.setValue(0)
            self.ui.ui_pages.repl_idx_spin.setDisabled(True)

        # necessary to avoid infinite recursion
        # otherwise GUI calls spin button actions
        # as their values are adjusted below
        sample_signals_blocked = self.ui.ui_pages.sample_idx_spin.blockSignals(True)
        repl_signals_blocked = self.ui.ui_pages.repl_idx_spin.blockSignals(True)

        # can we even circulate through something
        # (basically: non-deterministic nodes)
        if isinstance(node_dag.value, list):
            if isinstance(node_dag.value[0], pjdt.AnnotatedTree):
                if node_dag.repl_size >= 2:
                    _prepare_for_tree(potential_repl=True)

                else:
                    _prepare_for_tree(potential_repl=False)

            # scalar r.v.
            else:
                if not self.ui.ui_pages.all_samples_radio.isEnabled() and \
                        not self.ui.ui_pages.one_sample_radio.isEnabled():
                    if node_dag.repl_size > 1:
                        _prepare_for_scalar(potential_repl=True)
                    
                    else:
                        _prepare_for_scalar(potential_repl=False)

                # we always bundle replicates together, so
                # we never spin through them
                self.ui.ui_pages.repl_idx_spin.setMinimum(0)
                self.ui.ui_pages.repl_idx_spin.setValue(0)
                self.ui.ui_pages.repl_idx_spin.setDisabled(True) 

                # irrespective of all samples or one sample
                # it should always be possible to check one
                # sample if there are replicates
                if node_dag.repl_size > 1:
                    # radio
                    self.ui.ui_pages.one_sample_radio.setEnabled(True)
                    self.ui.ui_pages.one_sample_radio.setCheckable(True)
                    
                    # spin
                    self.ui.ui_pages.sample_idx_spin.setMinimum(1)
                    self.ui.ui_pages.sample_idx_spin.setValue(1)
                    self.ui.ui_pages.sample_idx_spin.setEnabled(True)
                
                # if there are no replicates, it does not make
                # sense to look at individual samples
                else:
                    # radio
                    self.ui.ui_pages.one_sample_radio.setDisabled(True)
                    self.ui.ui_pages.one_sample_radio.setCheckable(False)       

                    # spin #
                    self.ui.ui_pages.sample_idx_spin.setMinimum(0)
                    self.ui.ui_pages.sample_idx_spin.setValue(0)
                    self.ui.ui_pages.sample_idx_spin.setDisabled(True)

                # all samples checked, so we won't be spinning
                # through individuals samples
                if self.ui.ui_pages.all_samples_radio.isChecked():
                    self.ui.ui_pages.sample_idx_spin.setMinimum(0)
                    self.ui.ui_pages.sample_idx_spin.setValue(0)
                    self.ui.ui_pages.sample_idx_spin.setDisabled(True)

        # can't spin through anything! nothing to do!
        else:
            _nothing_to_spin_through()

        # need to unblock signals from here on
        self.ui.ui_pages.sample_idx_spin.blockSignals(sample_signals_blocked)
        self.ui.ui_pages.repl_idx_spin.blockSignals(repl_signals_blocked)

    def refresh_node_lists(self):
        # pgm page node list #
        pgm_node_list = \
            [node_name for node_name in
             self.gui_modeling.dag_obj.name_node_dict]
        self.ui.ui_pages.node_list.clear()  # clear first
        self.ui.ui_pages.node_list.addItems(pgm_node_list)

        # compare page node list
        compare_nodes_list = \
            [nd.node_name for nd in
             self.gui_modeling.dag_obj.get_sorted_node_dag_list() if
             nd.is_sampled]
        self.ui.ui_pages.compare_node_list.clear()
        self.ui.ui_pages.compare_node_list.addItems(compare_nodes_list)

        # coverage page node list
        coverage_nodes_list = \
            [nd.node_name for nd in
             self.gui_modeling.dag_obj.get_sorted_node_dag_list() if
             not nd.is_deterministic]
        self.ui.ui_pages.coverage_node_list.clear()
        self.ui.ui_pages.coverage_node_list.addItems(coverage_nodes_list)

    def refresh_selected_node_display_plot_radio(self):
        # if nodes have been created and selected #
        if self.ui.ui_pages.node_list.currentItem() is not None:
            self.do_selected_node_dag_page()

    def refresh_selected_node_display_plot_spin(self):
        # if nodes have been created and selected #
        if self.ui.ui_pages.node_list.currentItem() is not None:
            self.do_selected_node_dag_page(spin_buttons_clicked=True)

    def refresh_cmd_history(self, user_reset=False):
        if user_reset:
            self.gui_modeling.cmd_log_list.insert(0, "## Model is being cleared at this point, as a result of (i) reading script, (ii) loading model, or (iii) user reset\n")

        cmd_hist_str = self.gui_modeling.cmd_log().lstrip()
        self.ui.ui_pages.cmd_log_textbox.setText(cmd_hist_str)

    def clean_disable_everything(self, user_reset=False):
        if user_reset:
            cmd_hist_str: str = ""

            try:
                # update GUI cmd history
                self.gui_modeling.cmd_log_list.append(
                    ("\n## Model is being cleared at this point, as a result "
                     "of (i) reading script, (ii) loading model, or (iii) "
                     "user reset"))
                cmd_hist_str = self.gui_modeling.cmd_log()

            except Exception as e:
                print("An error occurred: ", type(e).__name__)
                cmd_hist_str = \
                    ("## Model is being cleared at this point, as a result "
                     "of (i) reading script, (ii) loading model, or (iii) "
                     "user reset")

            self.ui.ui_pages.cmd_log_textbox.setText(cmd_hist_str.lstrip())

        else:
            self.ui.ui_pages.cmd_log_textbox.clear()

        # pgm-related objects are reset
        self.gui_modeling.clear()
        self.refresh_resample_button()

        # remove all node names from list
        self.ui.ui_pages.node_list.clear()
        self.ui.ui_pages.compare_node_list.clear()
        self.ui.ui_pages.summary_stats_list.clear()
        self.ui.ui_pages.coverage_node_list.clear()
        self.ui.ui_pages.cov_summary_stats_list.clear()

        # reset text on display and summary
        # panels
        self.ui.ui_pages.values_content.clear()
        self.ui.ui_pages.summary_content.clear()

        # resetting and disabling buttons
        # radio #
        self.ui.ui_pages.one_sample_radio.setDisabled(True)
        self.ui.ui_pages.one_sample_radio.setCheckable(False)
        self.ui.ui_pages.one_sample_radio.setChecked(False)

        self.ui.ui_pages.all_samples_radio.setDisabled(True)
        self.ui.ui_pages.all_samples_radio.setCheckable(False)
        self.ui.ui_pages.all_samples_radio.setChecked(False)

        # spin #
        self.ui.ui_pages.sample_idx_spin.setMinimum(0)
        self.ui.ui_pages.sample_idx_spin.setValue(0)
        self.ui.ui_pages.sample_idx_spin.setDisabled(True)
        self.ui.ui_pages.repl_idx_spin.setMinimum(0)
        self.ui.ui_pages.repl_idx_spin.setValue(0)
        self.ui.ui_pages.repl_idx_spin.setDisabled(True)

        # reset plot(s)
        self.ui.ui_pages.pgm_page_matplotlib_widget.fig.clf()
        self.ui.ui_pages.pgm_page_matplotlib_widget.initialize_axes()
        self.ui.ui_pages.pgm_page_matplotlib_widget.fig.canvas.draw()
        self.ui.ui_pages.compare_page_matplotlib_widget.fig.clf()
        self.ui.ui_pages.compare_page_matplotlib_widget.initialize_axes()
        self.ui.ui_pages.compare_page_matplotlib_widget.fig.canvas.draw()
        self.ui.ui_pages.coverage_page_matplotlib_widget.fig.clf()
        self.ui.ui_pages.coverage_page_matplotlib_widget.initialize_axes()
        self.ui.ui_pages.coverage_page_matplotlib_widget.fig.canvas.draw()

    def one_sample_clicked(self):
        self.ui.ui_pages.one_sample_radio.setEnabled(True)
        self.ui.ui_pages.one_sample_radio.setChecked(True)
        self.ui.ui_pages.all_samples_radio.setEnabled(True)
        self.ui.ui_pages.all_samples_radio.setChecked(False)
        self.refresh_selected_node_display_plot_radio()

    def all_samples_clicked(self):
        self.ui.ui_pages.one_sample_radio.setEnabled(True)
        self.ui.ui_pages.one_sample_radio.setChecked(False)
        self.ui.ui_pages.all_samples_radio.setEnabled(True)
        self.ui.ui_pages.all_samples_radio.setChecked(True)

        # if we look at all samples, it does not make sense
        # to be able to look at each sample individually
        with QSignalBlocker(self.ui.ui_pages.sample_idx_spin):
            self.ui.ui_pages.sample_idx_spin.setMinimum(0)
            self.ui.ui_pages.sample_idx_spin.setValue(0)
            self.ui.ui_pages.sample_idx_spin.setDisabled(True)

        self.refresh_selected_node_display_plot_radio()

    # side-effect: sets is_avg_repl_check
    # member as True or False
    def avg_repl_check(self):
        chk_button = self.sender()
        self.is_avg_repl_check = chk_button.isChecked()

        # clearing summary stats
        # (should only be there
        # if averaging over repls)
        if not self.is_avg_repl_check:
            self.ui.ui_pages.summary_stats_list.clear()

        # clicking check box should de-select node
        # being compared to force user to select
        self.ui.ui_pages.compare_node_list.clearSelection()

    def quit_app_button_clicked(self):
        sys.exit()

    def print_about(self):
        self.ui.about_licensing.exec()


def call_gui():
    gui_app = QApplication(sys.argv)
    window = GUIMainWindow()

    timer = QTimer()
    timer.timeout.connect(lambda: None)
    timer.start(100)

    sys.exit(gui_app.exec())


# call GUI #
if __name__ == "__main__":
    call_gui()
