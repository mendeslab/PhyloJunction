import sys
from contextlib import ExitStack
import os
import re
from phylojunction.pgm.values import elements, path_text, shape, summaries, plot_values, component_values, component_frame
from phylojunction.pgm.execution import ExecutionResults
from phylojunction.readwrite.execution_export import write_results, write_maps
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
        self.results = DirectedAcyclicGraph().execute()
        self.cmd_log_list = []
        self.active_script_lines: ty.Optional[ty.List[str]] = []
        self._seed = None

    @property
    def dag_obj(self):
        return self.results.model

    def update_dag_random_seed(self, a_random_seed):
        self._seed = a_random_seed

    # Execute new commands transactionally for all current runs. A changed run count
    # takes effect on the next full rerun, so additions never invent missing prior values.
    def parse_cmd_update_pgm(self, cmd_line_list, gui_main_window_obj, clear_cmd_log_list=False):
        lines = [line.strip() for line in cmd_line_list if line.strip() and not line.strip().startswith('#')]
        try:
            if not self.dag_obj.definitions:
                seed = gui_main_window_obj._read_seed()
                candidate = cmdp.compile_model('\n'.join(cmd_line_list), '<interactive>').execute(
                    gui_main_window_obj.ui.ui_pages.execution_count.value(), seed)
            else:
                candidate = self.results.append('\n'.join(cmd_line_list))
        except Exception as error:
            gui_main_window_obj.ui.bottom_label_left.setText('Execution failed; previous results retained')
            gui_main_window_obj.ui.ui_pages.warnings_textbox.setText(str(error))
            return False
        self.results = candidate
        if clear_cmd_log_list:
            self.cmd_log_list = []
        self.cmd_log_list.extend(lines)
        if self.active_script_lines is not None:
            self.active_script_lines.extend(lines)
        gui_main_window_obj.refresh_resample_button()
        return True

    def cmd_log(self):
        return "\n".join(self.cmd_log_list) + "\n\n"

    def clear(self):
        self.results = DirectedAcyclicGraph().execute()
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

        self.ui.ui_pages.random_seed_prefix_textbox.returnPressed.connect(self._seed_entered)
        self.ui.ui_pages.tr_dag_node_textbox.setToolTip(
            "Comma-separated tree variable names, for example: trees, mapped_trees")

        # cmd line enter #
        self.ui.ui_pages.cmd_prompt.returnPressed.connect(
            self.parse_cmd_update_gui)

        # Selection changes include keyboard navigation as well as mouse clicks.
        self.ui.ui_pages.node_list.currentItemChanged.connect(
            lambda current, previous: self._display_selected_node())
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
        if not self.gui_modeling.parse_cmd_update_pgm(cmd_line_list, self):
            return

        # update GUI cmd history
        self.refresh_cmd_history()

        self.ui.ui_pages.cmd_prompt.clear()  # clear

        # updating node lists
        self.refresh_node_lists()

    # Select one explicit nested path, keeping the same path when comparing replicates.
    def selected_values(self, node_name, sample_idx, repl_idx, all_executions=False):
        values = self.gui_modeling.results.values(node_name)
        execution = max(0, sample_idx or 0)
        choices = list(elements(values[execution]))
        path, selected = choices[min(max(0, repl_idx), len(choices) - 1)]
        self.ui.ui_pages.repl_idx_spin.setToolTip(node_name + path_text(path))
        if not all_executions:
            return [selected], path
        if any(shape(value) != shape(values[0]) for value in values[1:]):
            return values, None
        return [dict(elements(value))[path] for value in values], path

    # Display observations separately from summaries across replicates. Tree summaries
    # remain available for individual trees; vector elements are never averaged together.
    def selected_node_display(self, node_dag, do_all_samples, sample_idx=None,
                              repl_idx=0, repl_size=1, display_reconstructed=False):
        values, path = self.selected_values(node_dag.node_name, sample_idx, repl_idx, do_all_samples)
        rendered = []
        for value in values:
            if isinstance(value, pjdt.AnnotatedTree) and display_reconstructed:
                value = value.extract_reconstructed_tree().as_string(schema='newick').strip()
            rendered.append(str(value))
        self.ui.ui_pages.values_content.setText('\n'.join(rendered))
        if len(values) == 1 and isinstance(values[0], pjdt.AnnotatedTree):
            text = tabulate(values[0].get_stats_dict().items(), tablefmt='plain')
        else:
            rows = summaries(values) if path is not None else []
            text = tabulate(rows, headers='keys', tablefmt='plain') if rows else ''
        self.ui.ui_pages.summary_content.setText(text)

    # Plot observations from one selected path, or a single selected tree.
    def selected_node_plot(self, fig_obj, fig_axes, node_dag, do_all_samples,
                           sample_idx=None, repl_idx=0, repl_size=1, draw_reconstructed=False):
        values, _ = self.selected_values(node_dag.node_name, sample_idx, repl_idx, do_all_samples)
        plot_values(fig_axes, values, draw_reconstructed)
        fig_obj.canvas.draw()

    # The first selector chooses a replicate; the second enumerates explicit leaf paths.
    def selected_node_read(self, node_name):
        results = self.gui_modeling.results
        index = max(0, min(self.ui.ui_pages.sample_idx_spin.value() - 1, results.n_executions - 1))
        count = len(list(elements(results.value(node_name, index))))
        return results.model.name_node_dict[node_name], results.n_executions, count

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

    def component_values(self, name):
        return component_values(self.gui_modeling.results, name)

    def component_frame(self, name):
        return component_frame(self.gui_modeling.results, name)

    # Compare the chosen component across replicates, without within-vector averages.
    def do_selected_node_compare_page(self):
        u = self.ui.ui_pages
        item = u.compare_node_list.currentItem()
        u.summary_stats_list.clear()
        if item is None:
            return
        name = item.text()
        self.pj_comparison_df = self.component_frame(name)
        self.pj_comparison_df['program'] = 'PJ'
        if isinstance(self.component_values(name)[0], pjdt.AnnotatedTree):
            u.summary_stats_list.addItems([column for column in self.pj_comparison_df
                                          if column not in ('sample', 'replicate', 'program')])

    # Coverage always compares one numeric component's value per model replicate.
    def do_selected_node_coverage_page(self):
        u = self.ui.ui_pages
        item = u.coverage_node_list.currentItem()
        u.cov_summary_stats_list.clear()
        if item is not None:
            self.coverage_df = self.component_frame(item.text())

    ##################
    # Functions for  #
    # reading and    #
    # loading events #
    ##################

    # Read configuration without changing RNG state. NumPy's existing seed API accepts uint32.
    def _read_seed(self):
        text = self.ui.ui_pages.random_seed_prefix_textbox.text().strip()
        if not text:
            return None
        if text.isascii() and text.isdecimal():
            value = int(text)
            if value <= 4294967295:
                return value
        raise ValueError("Enter a whole number from 0 to 4294967295, or leave the seed blank.")

    # Return checks the field but does not execute or apply a seed; Rerun reads it directly.
    def _seed_entered(self):
        try:
            seed = self._read_seed()
        except ValueError as error:
            QMessageBox.warning(self, "Invalid seed", str(error))
            return
        if seed is None:
            message = "No fixed seed; the next run will continue the random sequence."
        elif self.gui_modeling.active_script_lines:
            message = f"Seed {seed} is ready; click Rerun to apply it."
        else:
            message = f"Seed {seed} is ready; read a script or enter commands."
        self.ui.bottom_label_left.setText(message)

    # Describe every completed run without guessing whether external inputs or results changed.
    def _run_completed(self, seed):
        message = ("Completed without resetting the random seed." if seed is None else
                   f"Completed with seed {seed}. Same seed + inputs reproduce the same samples.")
        self.ui.bottom_label_left.setText(message)

    # Qt selection callbacks cannot propagate exceptions to their caller. Report display failures
    # here, after any execution-status message, without rerunning the model.
    def _display_selected_node(self):
        try:
            self.do_selected_node_dag_page()
        except Exception as error:
            self.ui.bottom_label_left.setText("Displaying the selected node failed; see Warnings")
            self.ui.ui_pages.warnings_textbox.setText(str(error))

    # Read and validate before cleanup; a cancelled chooser or bad input preserves current results.
    def read_execute_script(self):
        script_fp, _ = QFileDialog.getOpenFileName(
            parent=self, caption="Read script", dir=".", filter="*.pj")
        if not script_fp:
            return
        u = self.ui.ui_pages
        item = u.node_list.currentItem()
        selected_name = item.text() if item else None
        try:
            seed = self._read_seed()
            commands = pjread.read_text_file(script_fp, keep_empty=True)
            candidate = cmdp.compile_model('\n'.join(commands), script_fp).execute(u.execution_count.value(), seed)
        except Exception as error:
            self.ui.bottom_label_left.setText("Script failed; previous results retained")
            QMessageBox.warning(self, "Cannot read script", str(error))
            return

        self.clean_disable_everything()
        self.gui_modeling.results = candidate
        self.gui_modeling.active_script_lines = [line.strip() for line in commands
                                                if line.strip() and not line.strip().startswith('#')]
        self.gui_modeling.cmd_log_list = list(self.gui_modeling.active_script_lines)
        self.refresh_resample_button()
        self.refresh_cmd_history(user_reset=True)
        with QSignalBlocker(u.node_list):
            self.refresh_node_lists()
        self.show_pgm_page()

        nodes = list(self.gui_modeling.dag_obj.name_node_dict.values())
        chosen = next((n for n in nodes if n.node_name == selected_name), None)
        if chosen is None:
            chosen = next((n for n in nodes if n.is_sampled and any(isinstance(v, pjdt.AnnotatedTree)
                           for _, v in elements(self.gui_modeling.results.value(n.node_name)))), None)
        if chosen is None:
            chosen = next((n for n in nodes if n.is_sampled), nodes[0] if nodes else None)
        if self.gui_modeling.active_script_lines is None:
            self.ui.bottom_label_left.setText("Script reported errors; see Warnings")
        elif not nodes:
            self.ui.bottom_label_left.setText("Script completed; no model nodes.")
        else:
            self._run_completed(seed)
        if chosen is not None:
            matches = u.node_list.findItems(chosen.node_name, Qt.MatchFlag.MatchExactly)
            u.node_list.setCurrentItem(matches[0])

    # Availability depends on executable source, independently of seed configuration.
    # Preserve a more specific unavailable reason when supplied.
    def refresh_resample_button(self, unavailable_reason=None):
        button = self.ui.ui_pages.resample_model
        lines = self.gui_modeling.active_script_lines
        button.setEnabled(bool(lines))
        if lines:
            button.setToolTip("Rerun the current script and interactive additions using the seed in Settings.")
        elif lines is not None:
            button.setToolTip("Load a script or enter commands to enable rerunning.")
        elif unavailable_reason:
            button.setToolTip(unavailable_reason)

    # Evaluate into a separate DAG before replacing the displayed model. Retain source and
    # view preferences across the existing cleanup routine, which resets GUIModeling.
    def resample_current_script(self):
        u = self.ui.ui_pages
        if not self.gui_modeling.active_script_lines or not u.resample_model.isEnabled():
            return
        try:
            seed = self._read_seed()
        except ValueError as error:
            QMessageBox.warning(self, "Invalid seed", str(error))
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
                candidate = cmdp.compile_model("\n".join(lines), "<rerun>").execute(u.execution_count.value(), seed)
            except Exception as error:
                self.ui.bottom_label_left.setText("Rerun failed; previous results retained")
                u.warnings_textbox.setText(str(error))
                return

            # No callbacks may draw an intermediate selection or out-of-range sample.
            with ExitStack() as blocked:
                for widget in (u.node_list, u.sample_idx_spin, u.repl_idx_spin,
                               u.one_sample_radio, u.all_samples_radio, u.reconstructed_tree_check):
                    blocked.enter_context(QSignalBlocker(widget))
                self.clean_disable_everything()
                u.resample_model.setEnabled(False)
                self.gui_modeling.results = candidate
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
                self._run_completed(seed)
            except Exception as error:
                self.ui.bottom_label_left.setText("Rerun completed, but displaying the new model failed")
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
                    pjread.read_csv_tsv_into_dataframe(compare_csv_fp)

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
            self.gui_modeling.results, self.gui_modeling.cmd_log_list = \
                pjread.read_serialized_pgm(model_fp)
            if not isinstance(self.gui_modeling.results, ExecutionResults):
                raise ValueError("This saved model uses the old format; reload its script.")
            self.gui_modeling.active_script_lines = None
            self.refresh_resample_button("This saved model has no verified script. Load its script to rerun.")
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
            prefix = self.ui.ui_pages.filename_prefix_textbox.text()

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
            prefix = self.ui.ui_pages.filename_prefix_textbox.text()

        # pickling and saving PGM
        pjwrite.dump_serialized_pgm(
            pickle_fp,
            self.gui_modeling.results,
            self.gui_modeling.cmd_log_list,
            prefix=prefix)

    def write_data_to_dir(self):
        # get dir path
        data_out_dir = QFileDialog.getExistingDirectory(
            self, "Save in folder", "./", QFileDialog.ShowDirsOnly
        )

        prefix = self.ui.ui_pages.filename_prefix_textbox.text()

        if data_out_dir:
            write_results(data_out_dir, self.gui_modeling.results, prefix,
                          self.ui.ui_pages.output_format.currentData())

    # Cancel without side effects; the exporter validates all selections before writing.
    def write_smaps_to_dir(self):
        data_out_dir = QFileDialog.getExistingDirectory(
            self, "Save in folder", "./", QFileDialog.ShowDirsOnly
        )
        if not data_out_dir:
            return
        prefix = self.ui.ui_pages.filename_prefix_textbox.text()
        mapped_attr_name = self.ui.ui_pages.attr_name_textbox.text().strip()
        names = [name.strip() for name in self.ui.ui_pages.tr_dag_node_textbox.text().split(',')]
        try:
            write_maps(data_out_dir, self.gui_modeling.results, names, mapped_attr_name, prefix=prefix)
        except Exception as error:
            QMessageBox.warning(self, "Cannot export stochastic maps", str(error))

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

        values, _ = self.selected_values(node_dag.node_name, sample_idx, repl_idx, sample_idx is None)
        return plot_values(axes, values, draw_reconstructed)

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
            node_dag = self.gui_modeling.dag_obj.name_node_dict[selected_node_name.split('[')[0]]

            if not self.pj_comparison_df.empty and \
                    not self.other_comparison_df.empty:

                # scalar
                if isinstance(self.component_values(selected_node_name)[0], (int, float, np.float64)):
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
            node_dag = self.gui_modeling.dag_obj.name_node_dict[selected_node_name.split('[')[0]]

            if not self.coverage_df.empty and \
                    not self.hpd_df.empty:

                # scalar
                if isinstance(self.component_values(selected_node_name)[0], (str, int, float, np.float64)):
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

    # Block selection signals while setting ranges, so layout changes cannot resample
    # or produce intermediate plots. Trees remain a one-replicate display.
    def init_and_refresh_radio_spin(self, node_dag, sample_size, repl_size):
        u = self.ui.ui_pages
        with QSignalBlocker(u.sample_idx_spin), QSignalBlocker(u.repl_idx_spin):
            values = self.gui_modeling.results.values(node_dag.node_name)
            is_tree = any(isinstance(v, pjdt.AnnotatedTree) for _, v in elements(values[0]))
            same_shape = all(shape(value) == shape(values[0]) for value in values[1:])
            u.one_sample_radio.setEnabled(True)
            u.one_sample_radio.setCheckable(True)
            u.all_samples_radio.setCheckable(True)
            u.all_samples_radio.setEnabled(not is_tree and same_shape)
            if is_tree or not same_shape or not u.all_samples_radio.isChecked():
                u.one_sample_radio.setChecked(True)
            u.sample_idx_spin.setRange(1, sample_size)
            u.sample_idx_spin.setEnabled(u.one_sample_radio.isChecked())
            u.repl_idx_spin.setRange(1, repl_size)
            u.repl_idx_spin.setEnabled(repl_size > 1)

    def refresh_node_lists(self):
        # pgm page node list #
        pgm_node_list = \
            [node_name for node_name in
             self.gui_modeling.dag_obj.name_node_dict]
        self.ui.ui_pages.node_list.clear()  # clear first
        self.ui.ui_pages.node_list.addItems(pgm_node_list)

        # List fixed-shape components explicitly; ragged values have no automatic
        # correspondence for comparison/coverage, but remain selectable on the model page.
        comparison, coverage = [], []
        for name, node in self.gui_modeling.dag_obj.name_node_dict.items():
            values = self.gui_modeling.results.values(name)
            if any(shape(value) != shape(values[0]) for value in values[1:]):
                continue
            for path, value in elements(values[0]):
                label = name + path_text(path)
                if node.is_sampled and isinstance(value, (int, float, pjdt.AnnotatedTree)):
                    comparison.append(label)
                if isinstance(value, (int, float)):
                    coverage.append(label)
        self.ui.ui_pages.compare_node_list.clear()
        self.ui.ui_pages.compare_node_list.addItems(comparison)
        self.ui.ui_pages.coverage_node_list.clear()
        self.ui.ui_pages.coverage_node_list.addItems(coverage)

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
