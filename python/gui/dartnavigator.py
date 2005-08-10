# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dartnavigator.ui'
#
# Created: Mon Jun 27 13:47:50 2005
#      by: The PyQt User Interface Compiler (pyuic) 3.13
#
# WARNING! All changes made in this file will be lost!


from qt import *


class DartNavigatorBase(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QDialog.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("DartNavigatorBase")


        DartNavigatorBaseLayout = QGridLayout(self,1,1,11,6,"DartNavigatorBaseLayout")

        self.nextSigmaButton = QPushButton(self,"nextSigmaButton")

        DartNavigatorBaseLayout.addWidget(self.nextSigmaButton,1,0)

        self.nextPhiButton = QPushButton(self,"nextPhiButton")

        DartNavigatorBaseLayout.addWidget(self.nextPhiButton,0,0)

        self.nextAlphaButton = QPushButton(self,"nextAlphaButton")

        DartNavigatorBaseLayout.addWidget(self.nextAlphaButton,1,1)

        self.prevSigmaButton = QPushButton(self,"prevSigmaButton")

        DartNavigatorBaseLayout.addWidget(self.prevSigmaButton,1,2)

        self.prevPhiButton = QPushButton(self,"prevPhiButton")

        DartNavigatorBaseLayout.addWidget(self.prevPhiButton,0,2)

        self.continuousCheckBox = QCheckBox(self,"continuousCheckBox")

        DartNavigatorBaseLayout.addWidget(self.continuousCheckBox,0,1)

        self.dartLabel = QLabel(self,"dartLabel")
        self.dartLabel.setAlignment(QLabel.WordBreak | QLabel.AlignTop | QLabel.AlignLeft)

        DartNavigatorBaseLayout.addMultiCellWidget(self.dartLabel,2,2,0,2)

        self.languageChange()

        self.resize(QSize(585,245).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.setTabOrder(self.nextPhiButton,self.continuousCheckBox)
        self.setTabOrder(self.continuousCheckBox,self.prevPhiButton)
        self.setTabOrder(self.prevPhiButton,self.nextSigmaButton)
        self.setTabOrder(self.nextSigmaButton,self.nextAlphaButton)
        self.setTabOrder(self.nextAlphaButton,self.prevSigmaButton)


    def languageChange(self):
        self.setCaption(self.__tr("DartNavigator"))
        self.nextSigmaButton.setText(self.__tr("nextSigma"))
        self.nextPhiButton.setText(self.__tr("nextPhi"))
        self.nextAlphaButton.setText(self.__tr("nextAlpha"))
        self.prevSigmaButton.setText(self.__tr("prevSigma"))
        self.prevPhiButton.setText(self.__tr("prevPhi"))
        self.continuousCheckBox.setText(self.__tr("continuous"))
        self.dartLabel.setText(QString.null)


    def __tr(self,s,c = None):
        return qApp.translate("DartNavigatorBase",s,c)
