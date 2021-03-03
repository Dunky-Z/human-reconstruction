# import sys
#
# from ui import Ui_HumanEstimate
# from maya_widget import MayaviQWidget, myAction
# from PyQt5 import QtWidgets, QtCore
# from PyQt5.QtWidgets import QApplication,QMainWindow
#
# class Demo(QMainWindow,Ui_HumanEstimate):
#     def __init__(self):
#         super(Demo,self).__init__()
#
#         self.setUpMaya(self)
#
#
#         self.set_menu()
#         self.setDoubleSpinBox()
#         self.setHorizontalSlider()
#         self.setWindowTitle("Accurate human body reshape with anthropometry-constrained optimization")
#
#     def setUpMaya(self, parent=None):
#         super(Demo, self).setupUi(self)
#         self.maya_container = QtWidgets.QWidget(self.mayawindows)
#         self.maya_container.setWindowTitle("Embedding Mayavi in a PyQt5 Application")
#         layout = QtWidgets.QGridLayout(self.maya_container)
#         self.viewer3D = MayaviQWidget(self.maya_container)
#         layout.addWidget(self.viewer3D, 1, 1)
#         print("setUp running")
#
#     def set_menu(self):
#         super(Demo, self).setupUi(self)
#         menubar = self.menuBar()
#         fileMenu = menubar.addMenu('&File')
#
#         exit = QtWidgets.QAction("Exit", self)
#         exit.setShortcut("Ctrl+Q")
#         exit.setStatusTip('Exit application')
#         exit.triggered.connect(QApplication.quit)
#         fileMenu.addAction(exit)
#
#         save = QtWidgets.QAction("Save", self)
#         save.setShortcut("Ctrl+S")
#         save.setStatusTip('save obj file')
#         save.triggered.connect(self.viewer3D.save)
#         fileMenu.addAction(save)
#
#         # self.flag_ = 0
#         # self.label_ = "female"
#         # self.mode = {0: "global_mapping", 1: "local_with_mask", 2: "local_with_rfemat"}
#         # for i in range(0, len(self.mode)):
#         #     mode = myAction(i, self.mode[i], self)
#         #     mode.myact.connect(self.select_mode)
#         #     # self.connect(mode, QtCore.SIGNAL('myact(int)'), self.select_mode)
#         #     fileMenu.addAction(mode)
#         # self.setToolTip('This is a window, or <b>something</b>')
#
#     def select_mode(self, id):
#         self.flag_ = id
#         self.setWindowTitle(self.mode[id])
#         self.viewer3D.select_mode(label=self.label_, flag=self.flag_)
#
#     def setDoubleSpinBox(self):
#         self.doubleSpinBox01.setSingleStep(0.05)
#         self.doubleSpinBox01.valueChanged.connect(self.horizontalSlider01.setValue)
#     def setHorizontalSlider(self):
#         self.horizontalSlider01.setSingleStep(0.05)
#         self.horizontalSlider01.valueChanged.connect(self.doubleSpinBox01.setValue)
#
# class CamShow(QMainWindow,Ui_HumanEstimate):
#     def __init__(self,parent = None):
#         super(CamShow,self).__init__(parent)
#         self.setupUi(self)
#
#
# if __name__ == '__main__':
#  app = QApplication(sys.argv)
#  maya = Demo()
#  # maya.setUp()
#  ui=CamShow()
#  maya.show()
#  sys.exit(app.exec_())

import sys

from ui import Ui_HumanEstimate
from maya_widget import MayaviQWidget
from PyQt5 import QtWidgets
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QApplication,QMainWindow
import utils

class Demo(QMainWindow,Ui_HumanEstimate):
    def __init__(self):
        super(Demo,self).__init__()
        self.setUpMaya(self)
        #self.set_menu()

        self.set_radio()
        self.setDoubleSpinBox()
        self.setHorizontalSlider()
        #self.set_slider()

        self.flag_ = 0
        self.label_ = "female"

    def setUpMaya(self,parent = None):
        super(Demo, self).setupUi(self)
        self.maya_container = QtWidgets.QWidget(self.mayawindows)
        self.maya_container.setWindowTitle("Embedding Mayavi in a PyQt5 Application")
        layout  = QtWidgets.QGridLayout(self.maya_container)
        self.viewer3D = MayaviQWidget(self.maya_container)
        layout.addWidget(self.viewer3D, 1, 1)
        self.maya_container.show()
        print("setUp running")

    def set_menu(self):
        super(Demo, self).setupUi(self)
        fileMenu = self.menu

        exit = self.actionExit
        exit.setShortcut("Ctrl+Q")
        exit.setStatusTip('Exit application')
        exit.triggered.connect(QApplication.quit)
        fileMenu.addAction(exit)

        # save = QtWidgets.QAction("Save", self)
        # save.setShortcut("Ctrl+S")
        # save.setStatusTip('save obj file')
        # save.triggered.connect(self.viewer3D.save)
        # fileMenu.addAction(save)
    def find_attributes(self, name_start):
        return [value for name, value in (self.__dict__.items())
                          if name.startswith(name_start)]

    def set_slider(self):
        self.slider = self.find_attributes('horizontalSlider01_')
        self.spin = self.find_attributes('doubleSpinBox01_')
        for i in range(0, utils.M_NUM):
            slider = self.slider[i]
            slider.setStatusTip('%d. %s' % (i, utils.M_STR[i]))
            slider.setRange(-30, 30)
            slider.valueChanged.connect(self.viewer3D.sliderForwardedValueChangeHandler)

            spinBox = self.spin[i]
            spinBox.setRange(-30, 30)
            spinBox.valueChanged.connect(slider.setValue)
            slider.valueChanged.connect(spinBox.setValue)



    def setDoubleSpinBox(self):
        self.doubleSpinBox01_1.setSingleStep(0.05)
        self.doubleSpinBox01_1.valueChanged.connect(self.horizontalSlider01_1.setValue)
    def setHorizontalSlider(self):
        self.horizontalSlider01_1.setSingleStep(0.05)
        self.horizontalSlider01_1.valueChanged.connect(self.doubleSpinBox01_1.setValue)

    def set_radio(self):
        self.radio1 = self.radioButton_female
        self.radio2 = self.radioButton_male
        self.radio1.setChecked(True)
        self.radio1.toggled.connect(self.radio_act)
        self.radio2.toggled.connect(self.radio_act)

    def radio_act(self):
        if self.radio1.isChecked():
            self.label_ = 'female'
        else:
            self.label_ = 'male'
        self.viewer3D.select_mode(label=self.label_, flag=self.flag_)

    def set_button(self):
        self.reset_button = self.pushButton_reset
        self.reset_button.setStatusTip('reset input to mean value')
        self.reset_button.clicked.connect(self.reset)

        # self.pre_button = self.pushButton_predict
        # self.pre_button.setToolTip('model your own shape')
        # self.pre_button.clicked.connect(self.show_dialog)


    def reset(self):
      pass

class CamShow(QMainWindow,Ui_HumanEstimate):
    def __init__(self,parent = None):
        super(CamShow,self).__init__(parent)
        self.setupUi(self)


if __name__ == '__main__':
 app = QApplication(sys.argv)
 maya = Demo()
 # maya.setUp()
 ui=CamShow()
 maya.show()
 sys.exit(app.exec_())