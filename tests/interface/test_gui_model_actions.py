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
