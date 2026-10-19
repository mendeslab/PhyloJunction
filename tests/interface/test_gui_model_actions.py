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
        self.controls.execution_count.setValue(3)

    # Use real commands so the controls exercise the same node types as a loaded script.
    def enter(self, text):
        self.controls.cmd_prompt.setText(text)
        self.window.parse_cmd_update_gui()

    def test_display_controls_redraw_without_sampling(self):
        self.enter('x ~ iid(2,normal(mean=0, sd=1))')
        self.enter('y ~ iid(2,normal(mean=0, sd=1))')
        w, u = self.window, self.controls
        dag = w.gui_modeling.dag_obj
        values = list(w.gui_modeling.results.values('x'))
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
        self.assertEqual(w.gui_modeling.results.values('x'), values)

        self.enter('birth := quasse_constant(rate=0)')
        self.enter('death := quasse_constant(rate=0)')
        self.enter('tr ~ iid(2,quasse(birth_rate=birth,death_rate=death,'
                   'stop="age",stop_value=.1,diffusion=.1))')
        u.node_list.setCurrentRow(u.node_list.count() - 1)
        with patch.object(w, 'selected_node_plot', wraps=w.selected_node_plot) as draw:
            u.sample_idx_spin.setValue(3)
            u.repl_idx_spin.setValue(2)
            u.reconstructed_tree_check.setChecked(False)
            self.assertEqual(draw.call_count, 3)
            self.assertEqual((u.sample_idx_spin.value(), u.repl_idx_spin.value()), (3, 2))
            self.assertFalse(draw.call_args.kwargs['draw_reconstructed'])

    # Replay must replace results only after success and must never replay history annotations.
    # Parser tests do not cover these GUI state transitions or seeded replay and invalid-input handling.
    def test_resample_replays_commands_and_preserves_view(self):
        import random
        import numpy as np
        from phylojunction.interface.pysidegui.pj_gui import cmdp

        w, u = self.window, self.controls
        commands = ['x ~ iid(2,normal(mean=0, sd=1))', 'z <- 2']
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
                patch.object(cmdp, 'compile_model', wraps=cmdp.compile_model) as execute:
            python_state, numpy_state = random.getstate(), np.random.get_state()
            u.random_seed_prefix_textbox.setText('invalid')
            u.resample_model.click()
            warning.assert_called_once()
            execute.assert_not_called()
            self.assertIs(w.gui_modeling.dag_obj, dag)
            self.assertEqual(u.values_content.toPlainText(), old_text)
            self.assertEqual(random.getstate(), python_state)
            np.testing.assert_equal(np.random.get_state(), numpy_state)
            u.random_seed_prefix_textbox.setText('12')
            u.resample_model.click()
            seeded_values = list(w.gui_modeling.results.values('x'))
            u.resample_model.click()
            self.assertEqual(w.gui_modeling.results.values('x'), seeded_values)
            self.assertIn('Completed with seed 12', w.ui.bottom_label_left.text())
            self.assertEqual(warning.call_count, 1)  # Only the invalid input warns.
            execute.reset_mock()
            u.random_seed_prefix_textbox.setText('   ')
            with patch.object(w, 'selected_node_plot', wraps=w.selected_node_plot) as draw:
                with patch('numpy.random.seed') as numpy_seed, patch('random.seed') as python_seed:
                    u.resample_model.click()
                    numpy_seed.assert_not_called()
                    python_seed.assert_not_called()
                self.assertEqual(draw.call_count, 1)
            execute.assert_called_once_with('\n'.join(commands))
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
        with patch.object(cmdp, 'compile_model', side_effect=ValueError('simulation failed')):
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
        # An addition is atomic across its commands and executions. Failure preserves
        # both values and replayable source; later valid commands can still be appended.
        w.gui_modeling.parse_cmd_update_pgm(['a <- 1', 'invalid'], w, clear_cmd_log_list=True)
        self.assertNotIn('a <- 1', w.gui_modeling.cmd_log_list)
        self.assertEqual(w.gui_modeling.active_script_lines, [])
        self.enter('b <- 2')
        self.assertTrue(u.resample_model.isEnabled())

    # Settings commit only configuration; actual execution validates before destroying results.
    # These interactions are not covered by parser or layout tests.
    def test_single_line_settings_and_script_feedback(self):
        import random
        import numpy as np
        from PySide6.QtWidgets import QLineEdit
        w, u = self.window, self.controls
        w.show_settings_page()
        python_state, numpy_state = random.getstate(), np.random.get_state()
        for field in (u.random_seed_prefix_textbox, u.filename_prefix_textbox,
                      u.tr_dag_node_textbox, u.attr_name_textbox):
            self.assertIsInstance(field, QLineEdit)
            field.setFocus()
            QTest.keyClicks(field, '123')
            QTest.keyClick(field, Qt.Key_Return)
            self.assertEqual(field.text(), '123')
        self.assertEqual(random.getstate(), python_state)
        np.testing.assert_equal(np.random.get_state(), numpy_state)
        self.assertFalse(w.gui_modeling.dag_obj.name_node_dict)
        for text, expected in [(' ', None), ('0', 0), ('0001', 1), ('4294967295', 4294967295)]:
            u.random_seed_prefix_textbox.setText(text)
            self.assertEqual(w._read_seed(), expected)
        for text in ['-1', '4294967296', '1.2', '1e3', '+1', 'abc']:
            u.random_seed_prefix_textbox.setText(text)
            with self.assertRaises(ValueError):
                w._read_seed()
        u.random_seed_prefix_textbox.setText('123')
        commands = ['c <- 2', 'x ~ iid(2,normal(mean=0, sd=1))']
        with patch('phylojunction.interface.pysidegui.pj_gui.QFileDialog.getOpenFileName',
                   return_value=('model.pj', '')), \
                patch('phylojunction.interface.pysidegui.pj_gui.pjread.read_text_file',
                      return_value=commands) as read:
            with patch.object(w, 'selected_node_plot', wraps=w.selected_node_plot) as draw:
                w.read_execute_script()
                self.assertEqual(draw.call_count, 1)
            self.assertIs(w.ui.pages.currentWidget(), u.pgm_page)
            self.assertEqual(u.node_list.currentItem().text(), 'x')
            values = list(w.gui_modeling.results.values('x'))
            u.node_list.setCurrentRow(0)
            w.read_execute_script()
            self.assertEqual(u.node_list.currentItem().text(), 'c')
            self.assertEqual(w.gui_modeling.results.values('x'), values)
            dag = w.gui_modeling.dag_obj
            with patch('phylojunction.interface.pysidegui.pj_gui.QMessageBox.warning') as warning:
                u.random_seed_prefix_textbox.setText('bad')
                with patch.object(w, 'clean_disable_everything') as cleanup:
                    w.read_execute_script()
                    cleanup.assert_not_called()
                self.assertIs(w.gui_modeling.dag_obj, dag)
                u.random_seed_prefix_textbox.clear()
                read.side_effect = OSError('unreadable file')
                w.read_execute_script()
                self.assertIs(w.gui_modeling.dag_obj, dag)
                self.assertEqual(warning.call_count, 2)
            read.side_effect = None
            read.return_value = []
            w.read_execute_script()
            self.assertIn('no model nodes', w.ui.bottom_label_left.text())
            read.return_value = ['a <- 1', 'invalid']
            with patch('phylojunction.interface.pysidegui.pj_gui.QMessageBox.warning'):
                w.read_execute_script()
            self.assertIn('Script failed', w.ui.bottom_label_left.text())
            self.assertFalse(u.resample_model.isEnabled())
        w.clean_disable_everything()
        u.random_seed_prefix_textbox.setText('bad')
        with patch('phylojunction.interface.pysidegui.pj_gui.QMessageBox.warning'):
            self.enter('a <- 1')
        self.assertEqual(w.gui_modeling.active_script_lines, [])
        self.assertEqual(u.cmd_prompt.text(), 'a <- 1')
