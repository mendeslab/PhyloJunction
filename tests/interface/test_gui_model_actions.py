import unittest
from unittest.mock import patch

from PySide6.QtCore import Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

from phylojunction.interface.pysidegui.pj_gui import GUIMainWindow


class TestModelActions(unittest.TestCase):
    # Exercise application signal wiring absent from layout/parser tests. Keep these checks
    # while the GUI owns selection and replay; move them if that responsibility moves elsewhere.
    def setUp(self):
        self.app = QApplication.instance() or QApplication([])
        self.window = GUIMainWindow()
        self.addCleanup(self.window.close)
        self.controls = self.window.ui.ui_pages

    # Use real commands so the controls exercise the same node types as a loaded script.
    def enter(self, text):
        self.controls.cmd_prompt.setText(text)
        self.window.parse_cmd_update_gui()

    def test_display_controls_redraw_without_sampling(self):
        self.enter('x ~ normal(n=3, nr=2, mean=0, sd=1)')
        self.enter('y ~ normal(n=3, nr=2, mean=0, sd=1)')
        w, u = self.window, self.controls
        dag = w.gui_modeling.dag_obj
        values = list(dag.name_node_dict['x'].value)
        with patch.object(w, 'selected_node_plot', wraps=w.selected_node_plot) as draw:
            u.node_list.setCurrentRow(0)
            self.assertEqual(draw.call_count, 1)
            QTest.keyClick(u.node_list, Qt.Key_Down)
            self.assertEqual(u.node_list.currentItem().text(), 'y')
            self.assertEqual(draw.call_count, 2)
            QTest.mouseClick(u.node_list.viewport(), Qt.LeftButton,
                             pos=u.node_list.visualItemRect(u.node_list.item(0)).center())
            self.assertEqual(draw.call_count, 3)
            u.all_samples_radio.click()
            self.assertEqual(draw.call_count, 4)
            u.one_sample_radio.click()
            self.assertEqual(draw.call_count, 5)
            u.sample_idx_spin.setValue(2)
            self.assertEqual(draw.call_count, 6)
            u.reconstructed_tree_check.setChecked(True)
            self.assertEqual(draw.call_count, 7)
            self.assertEqual(u.sample_idx_spin.value(), 2)
        self.assertIs(w.gui_modeling.dag_obj, dag)
        self.assertEqual(dag.name_node_dict['x'].value, values)

        self.enter('birth := quasse_constant(rate=0)')
        self.enter('death := quasse_constant(rate=0)')
        self.enter('tr ~ quasse(n=3,nr=2,birth_rate=birth,death_rate=death,'
                   'stop="age",stop_value=.1,diffusion=.1)')
        u.node_list.setCurrentRow(u.node_list.count() - 1)
        with patch.object(w, 'selected_node_plot', wraps=w.selected_node_plot) as draw:
            u.sample_idx_spin.setValue(3)
            u.repl_idx_spin.setValue(2)
            u.reconstructed_tree_check.setChecked(False)
            self.assertEqual(draw.call_count, 3)
            self.assertEqual((u.sample_idx_spin.value(), u.repl_idx_spin.value()), (3, 2))
            self.assertFalse(draw.call_args.kwargs['draw_reconstructed'])

    # Replay must replace results only after success and must never replay history annotations.
    # Parser tests do not cover these GUI state transitions or the no-op seed guard.
    def test_resample_replays_commands_and_preserves_view(self):
        import random
        import numpy as np
        from phylojunction.interface.pysidegui.pj_gui import cmdp

        w, u = self.window, self.controls
        commands = ['x ~ normal(n=3, nr=2, mean=0, sd=1)', 'z <- 2']
        # Exercise the same script-load path as File > Read script, followed by prompt input.
        with patch('phylojunction.interface.pysidegui.pj_gui.QFileDialog.getOpenFileName',
                   return_value=('model.pj', '')), \
                patch('phylojunction.interface.pysidegui.pj_gui.pjread.read_text_file',
                      return_value=[commands[0]]):
            w.read_execute_script()
        self.enter(commands[1])
        u.node_list.setCurrentRow(0)
        u.one_sample_radio.click()
        u.sample_idx_spin.setValue(2)
        u.reconstructed_tree_check.setChecked(True)
        dag = w.gui_modeling.dag_obj
        history = w.gui_modeling.cmd_log_list[:]
        old_text = u.values_content.toPlainText()
        self.assertEqual(w.gui_modeling.active_script_lines, commands)
        with patch('phylojunction.interface.pysidegui.pj_gui.QMessageBox.warning') as warning, \
                patch.object(cmdp, 'script2dag', wraps=cmdp.script2dag) as execute:
            python_state, numpy_state = random.getstate(), np.random.get_state()
            u.random_seed_prefix_textbox.setText('12')
            u.resample_model.click()
            warning.assert_called_once()
            execute.assert_not_called()
            self.assertIs(w.gui_modeling.dag_obj, dag)
            self.assertEqual(u.values_content.toPlainText(), old_text)
            self.assertEqual(random.getstate(), python_state)
            np.testing.assert_equal(np.random.get_state(), numpy_state)
            u.random_seed_prefix_textbox.setText(' \n ')
            with patch.object(w, 'selected_node_plot', wraps=w.selected_node_plot) as draw:
                u.resample_model.click()
                self.assertEqual(draw.call_count, 1)
            execute.assert_called_once_with('\n'.join(commands), in_pj_file=False, random_seed=None)
        self.assertIsNot(w.gui_modeling.dag_obj, dag)
        self.assertEqual(u.node_list.currentItem().text(), 'x')
        self.assertEqual(u.sample_idx_spin.value(), 2)
        self.assertTrue(u.one_sample_radio.isChecked())
        self.assertTrue(u.reconstructed_tree_check.isChecked())
        self.assertEqual(w.gui_modeling.cmd_log_list, history)
        self.assertEqual(w.gui_modeling.active_script_lines, commands)
        u.all_samples_radio.click()
        u.resample_model.click()
        self.assertTrue(u.all_samples_radio.isChecked())
        self.assertFalse(u.sample_idx_spin.isEnabled())
        self.assertEqual(w.gui_modeling.active_script_lines, commands)
        self.assertEqual(w.gui_modeling.cmd_log_list, history)

        dag = w.gui_modeling.dag_obj
        old_text = u.values_content.toPlainText()
        old_artists = list(u.pgm_page_matplotlib_widget.axes.get_children())
        with patch.object(cmdp, 'script2dag', side_effect=ValueError('simulation failed')):
            u.resample_model.click()
        self.assertIs(w.gui_modeling.dag_obj, dag)
        self.assertEqual(u.values_content.toPlainText(), old_text)
        self.assertEqual(list(u.pgm_page_matplotlib_widget.axes.get_children()), old_artists)
        self.assertEqual(w.gui_modeling.active_script_lines, commands)
        self.assertTrue(u.resample_model.isEnabled())
        self.assertIsNone(QApplication.overrideCursor())
        self.assertIn('simulation failed', u.warnings_textbox.toPlainText())

        w.clean_disable_everything(user_reset=True)
        self.assertFalse(u.resample_model.isEnabled())
        # A failure following a valid command must neither duplicate that command nor allow
        # replay of a possibly partially changed DAG. Later successes cannot repair the source.
        w.gui_modeling.parse_cmd_update_pgm(['a <- 1', 'invalid'], w, clear_cmd_log_list=True)
        self.assertEqual(w.gui_modeling.cmd_log_list, ['a <- 1'])
        self.assertIsNone(w.gui_modeling.active_script_lines)
        self.enter('b <- 2')
        self.assertFalse(u.resample_model.isEnabled())
