import unittest

from PySide6.QtWidgets import QApplication
from phylojunction.interface.pysidegui.pj_gui import GUIMainWindow


class TestControlFlow(unittest.TestCase):
    # The real control container must reserve wrapped height and retain state when resized.
    # This covers application integration absent from tree tests, not superqt's algorithm.
    def test_wrap_and_unwrap(self):
        app = QApplication.instance() or QApplication([])
        window = GUIMainWindow()
        controls = window.ui.ui_pages
        container = controls.plot_controls
        container.setParent(None)
        widgets = [controls.one_sample_radio, controls.all_samples_radio,
                   controls.reconstructed_tree_check, controls.sample_idx_spin,
                   controls.repl_idx_spin, controls.save_pgm_node_plot]
        controls.reconstructed_tree_check.setChecked(True)
        controls.sample_idx_spin.setValue(12)
        layout = container.layout()
        width = sum(w.sizeHint().width() for w in widgets) + 6 * (len(widgets) - 1) + 10
        try:
            container.show()
            for available, wraps in [(width, False), (width // 2, True), (width, False)]:
                container.resize(available, layout.heightForWidth(available))
                app.processEvents()
                self.assertEqual(len({w.y() for w in widgets}) > 1, wraps)
                for i, widget in enumerate(widgets):
                    self.assertTrue(container.rect().contains(widget.geometry()))
                    for other in widgets[i + 1:]:
                        self.assertFalse(widget.geometry().intersects(other.geometry()))
                self.assertTrue(controls.reconstructed_tree_check.isChecked())
                self.assertEqual(controls.sample_idx_spin.value(), 12)
        finally:
            container.close()
            window.close()
