# Form implementation generated from reading ui file 'dartnavigator.ui'
#
# Created: Fri Mar 11 02:47:11 2005
#      by: The PyQt User Interface Compiler (pyuic)
#
# WARNING! All changes made in this file will be lost!


from qt import *


class DartNavigatorBase(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QDialog.__init__(self,parent,name,modal,fl)

        if name == None:
            self.setName("DartNavigatorBase")

        self.resize(418,133)
        self.setCaption(self.trUtf8("DartNavigator"))

        DartNavigatorBaseLayout = QGridLayout(self,1,1,11,6,"DartNavigatorBaseLayout")

        self.nextSigmaButton = QPushButton(self,"nextSigmaButton")
        self.nextSigmaButton.setText(self.trUtf8("nextSigma"))

        DartNavigatorBaseLayout.addWidget(self.nextSigmaButton,1,0)

        self.nextPhiButton = QPushButton(self,"nextPhiButton")
        self.nextPhiButton.setText(self.trUtf8("nextPhi"))

        DartNavigatorBaseLayout.addWidget(self.nextPhiButton,0,0)

        self.nextAlphaButton = QPushButton(self,"nextAlphaButton")
        self.nextAlphaButton.setText(self.trUtf8("nextAlpha"))

        DartNavigatorBaseLayout.addWidget(self.nextAlphaButton,1,1)

        self.prevSigmaButton = QPushButton(self,"prevSigmaButton")
        self.prevSigmaButton.setText(self.trUtf8("prevSigma"))

        DartNavigatorBaseLayout.addWidget(self.prevSigmaButton,1,2)

        self.prevPhiButton = QPushButton(self,"prevPhiButton")
        self.prevPhiButton.setText(self.trUtf8("prevPhi"))

        DartNavigatorBaseLayout.addWidget(self.prevPhiButton,0,2)

        self.continuousCheckBox = QCheckBox(self,"continuousCheckBox")
        self.continuousCheckBox.setText(self.trUtf8("continuous"))

        DartNavigatorBaseLayout.addWidget(self.continuousCheckBox,0,1)

        self.dartLabel = QLabel(self,"dartLabel")
        self.dartLabel.setText(self.trUtf8(""))
        self.dartLabel.setAlignment(QLabel.WordBreak | QLabel.AlignTop | QLabel.AlignLeft)

        DartNavigatorBaseLayout.addMultiCellWidget(self.dartLabel,2,2,0,2)

        self.setTabOrder(self.nextPhiButton,self.continuousCheckBox)
        self.setTabOrder(self.continuousCheckBox,self.prevPhiButton)
        self.setTabOrder(self.prevPhiButton,self.nextSigmaButton)
        self.setTabOrder(self.nextSigmaButton,self.nextAlphaButton)
        self.setTabOrder(self.nextAlphaButton,self.prevSigmaButton)
