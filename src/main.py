import vtkmodules.all as vtk
import sys
from PyQt5 import QtWidgets, QtCore
import organizer
from boolean import boolean_interface
import os
from mu3d.mu3dpy.mu3d import Graph

from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

class Ui_MainWindow(object):

    dpOcclusion = 0.1
    dpNOP = 10

    filenames = []

    dirname = os.path.dirname(__file__)

    def setupUi(self, MainWindow,org,ren):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1920, 1080)
        self.centralWidget = QtWidgets.QWidget(MainWindow)
        self.gridlayout = QtWidgets.QGridLayout(self.centralWidget)
        self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
        self.gridlayout.setSpacing(3)
        ##For flattening.
        self.pickedIds = [[[]],[[]],[[]]]
        self.countOfPickedRegions = 0
        self.inflateStruc = []

        #      init Widgets
        resetCameraButton = QtWidgets.QPushButton("Reset Camera")

        progressBar = QtWidgets.QProgressBar()

        def resetCamera():
            org.resetCamera()
            org.swapViewports(False)
            self.vtkWidget.update()

        resetCameraButton.clicked.connect(resetCamera)

        depthPeelingCheck = QtWidgets.QCheckBox("Depth Peeling")
        depthPeelingCheck.setChecked(True)

        def toggleDepthPeeling():
            if depthPeelingCheck.isChecked():
                ren.SetUseDepthPeeling(True)
                ren.SetOcclusionRatio(self.dpOcclusion)
                ren.SetMaximumNumberOfPeels(self.dpNOP)
            else:
                ren.SetUseDepthPeeling(False)
            self.vtkWidget.update()

        depthPeelingCheck.clicked.connect(toggleDepthPeeling)
        toggleDepthPeeling()

        dpOcclusion = QtWidgets.QLineEdit("0.1")
        dpNOP = QtWidgets.QLineEdit("10")

        def changeDpOcSettings():
            oc = dpOcclusion.text()
            self.dpOcclusion = float(oc)
            ren.SetOcclusionRatio(self.dpOcclusion)
            org.occlusion = self.dpOcclusion

        def changeDpNopSettings():
            nop = dpNOP.text()
            self.dpNOP = int(nop)
            ren.SetMaximumNumberOfPeels(self.dpNOP)
            org.numberOfPeels = self.dpNOP

        dpOcclusion.returnPressed.connect(changeDpOcSettings)
        dpNOP.returnPressed.connect(changeDpNopSettings)

        multiplyingButton = QtWidgets.QPushButton("Multiply")

        def onMultiply():
            org.onMultiply(depthPeelingCheck.isChecked(),filterEnabled.isChecked(),brightenCheck.isChecked())
            org.swapViewports(True)

            self.vtkWidget.update()

        multiplyingButton.clicked.connect(onMultiply)

        multiplyingCheck = QtWidgets.QCheckBox("Multiply")

        def toggleMultiply():
            if multiplyingCheck.isChecked():
                UpdateColorFilter.cam = ren.GetActiveCamera()
                UpdateColorFilter.sr = org
                UpdateColorFilter.dp = depthPeelingCheck.isChecked()
                UpdateColorFilter.filter = filterEnabled.isChecked()
                UpdateColorFilter.brighten = brightenCheck.isChecked()
                UpdateColorFilter.run = True
                ren.GetRenderWindow().GetInteractor().AddObserver("EndInteractionEvent", UpdateColorFilter)
                # self.vtkWidget.update()
            else :
                org.swapViewports(False)
                UpdateColorFilter.run = False

        multiplyingCheck.clicked.connect(toggleMultiply)

        colorDialog = QtWidgets.QColorDialog()

        filterDialog = QtWidgets.QColorDialog()
        filterDialog.hide()
        filterButton = QtWidgets.QPushButton("Filter Color")
        filterEnabled = QtWidgets.QCheckBox("Filter")

        def onFilterColorButton():
            filterDialog.hide()
            filterDialog.show()

        def onFilterColorChange():
            red = colorDialog.currentColor().red()# / 255
            green = colorDialog.currentColor().green()# / 255
            blue = colorDialog.currentColor().blue()# / 255
            #org.changeColor([red, green, blue], idx)
            org.colorFilterImage([red,green,blue])
            self.vtkWidget.update()

        filterButton.clicked.connect(onFilterColorButton)
        filterDialog.currentColorChanged.connect(onFilterColorChange)

        long = QtWidgets.QLineEdit("long")
        lat = QtWidgets.QLineEdit("lat")

        def onBrightMuliplication():
            org.brightenMultiplication()

        brightMultiplicationButton = QtWidgets.QPushButton("Bright Multiplication")
        brightMultiplicationButton.clicked.connect(onBrightMuliplication)

        brightenCheck = QtWidgets.QCheckBox("Brighten")
        brightenCheck.setChecked(True)

        addFileButton = QtWidgets.QPushButton("Add Files")

        def getFile():
            dlg = QtWidgets.QFileDialog()
            dlg.setFileMode(QtWidgets.QFileDialog.ExistingFiles)

            filename = os.path.join(self.dirname, "../Meshes")
            dlg.setDirectory(filename)
            #dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
            #dlg.setFilter("*.stl")
            filenames = []

            if dlg.exec_():
                for i in dlg.selectedFiles():
                    self.filenames.append(i)
                    addMesh(i)
                    idx = self.filenames.index(i)
                    org.addActor(i)
                    if idx % 3 == 0:
                        org.changeColor([0.0, 1.0, 1.0], idx)
                    elif idx % 3 == 1:
                        org.changeColor([1.0, 0.0, 1.0], idx)
                    elif idx % 3 == 2:
                        org.changeColor([1.0, 1.0, 0.0], idx)

            self.vtkWidget.GetRenderWindow().Render()

        addFileButton.clicked.connect(getFile)

        imageName = QtWidgets.QLineEdit("Image Name")

        def onSaveVpToImage():
            wti = vtk.vtkWindowToImageFilter()
            win = self.vtkWidget.GetRenderWindow()
            wti.SetInput(win)
            wti.SetInputBufferTypeToRGB()
            wti.ReadFrontBufferOff()

            writer = vtk.vtkPNGWriter()
            filename = os.path.join(self.dirname, "../tempResources/{}.png")
            writer.SetFileName(filename.format(imageName.text()))
            writer.SetInputConnection(wti.GetOutputPort())
            writer.Write()

        unfoldIterationsTextfield = QtWidgets.QLineEdit()
        unfoldIterationsTextfield.setText("100")

        saveVPtoImageButton = QtWidgets.QPushButton("Save as Image")
        saveVPtoImageButton.clicked.connect(onSaveVpToImage)

        def onCreatePaperMesh():
            try:
                iterations = int(unfoldIterationsTextfield.text())
            except:
                iterations = 100
            org.createPaperMeshPass(iterations,self.inflateStruc)
            self.vtkWidget.update()

        createPaperMeshButton = QtWidgets.QPushButton("Create Papermesh")
        createPaperMeshButton.clicked.connect(onCreatePaperMesh)


        resolutionWidth = QtWidgets.QLineEdit()
        resolutionWidth.setText("500")

        def onProjectPerTriangle():
            try:
                width = int(resolutionWidth.text())
            except:
                width = 500
            org.projectPass(self.inflateStruc, resolution = [width, width])
            #org.projectPassTemp()
            self.vtkWidget.update()
#            self.centralWidget.update()

        projectPerTriangle = QtWidgets.QPushButton("Project")
        projectPerTriangle.clicked.connect(onProjectPerTriangle)

        importTextfield = QtWidgets.QLineEdit()
        importTextfield.setText("model")

        def onImport():
            org.importUnfoldedMeshPass(self.inflateStruc,importTextfield.text())
            self.vtkWidget.GetRenderWindow().Render()
        importButton = QtWidgets.QPushButton("Import .obj Papermesh")
        importButton.clicked.connect(onImport)

        def onFlatten():
            #for i in range(len(self.pickedIds[0])):
            self.resultRenderes = org.onFlatten(self.pickedIds,self.inflateStruc)
            self.countOfPickedRegions = 0
            self.vtkWidget.GetRenderWindow().Render()
            for i in self.pickedIds:
                i.clear()
                i.append([])

        flattenButton = QtWidgets.QPushButton("Flatten")
        flattenButton.clicked.connect(onFlatten)

        def onCreateDedicatedPaperMesh():
            org.createDedicatedPaperMeshesPass(self.inflateStruc)
        createDedicatedPaperMesh_Button = QtWidgets.QPushButton("Dedicated Papermesh")
        createDedicatedPaperMesh_Button.clicked.connect(onCreateDedicatedPaperMesh)

        def onFinish():
            org.finish()
        finishButton = QtWidgets.QPushButton("Finish")
        finishButton.clicked.connect(onFinish)

        def onRegionSelection():
            for i in range(3):
                self.pickedIds[i].append([])
            self.countOfPickedRegions += 1
        selectRegionButton = QtWidgets.QPushButton("Select Regions")
        selectRegionButton.clicked.connect(onRegionSelection)

        booleanButton = QtWidgets.QPushButton("Boolean")

        def onBoolean():
            org.boolean()
        booleanButton.clicked.connect(onBoolean)

        hierarchySlider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        hierarchySlider.setRange(1, 4)

        def onHierarchySlider(value):
            org.draw_level(value)
            self.vtkWidget.GetRenderWindow().Render()

        hierarchySlider.valueChanged.connect(onHierarchySlider)

        hierarchical_difference_button = QtWidgets.QPushButton("Hierarchial Difference")

        def on_hierarchical_difference():
            org.hierarchical_difference()

        hierarchical_difference_button.clicked.connect(on_hierarchical_difference)

        testButton = QtWidgets.QPushButton("TestUnfold")
        def onUnfoldTest():
            org.onUnfoldTest()

        testButton.clicked.connect(onUnfoldTest)

#       Layout  --------------------------------------
        groupRight = QtWidgets.QGroupBox()
        layoutRight = QtWidgets.QVBoxLayout()
        groupRight.setLayout(layoutRight)

        groupLeft = QtWidgets.QGroupBox()
        layoutLeft = QtWidgets.QVBoxLayout()
        groupLeft.setLayout(layoutLeft)

        #DethPeeling Settings
        dpGroup = QtWidgets.QGroupBox()
        layoutDp = QtWidgets.QHBoxLayout()
        dpGroup.setLayout(layoutDp)

        layoutDp.addWidget(depthPeelingCheck)
        layoutDp.addWidget(dpOcclusion)
        layoutDp.addWidget(dpNOP)

        #previewGroup
        previewGroup = QtWidgets.QGroupBox()
        previewGroupLayout = QtWidgets.QVBoxLayout()
        previewGroup.setLayout(previewGroupLayout)

        previewGroupLayout.addWidget(multiplyingButton)
        previewGroupLayout.addWidget(brightenCheck)
        previewGroupLayout.addWidget(multiplyingCheck)
        previewGroupLayout.addWidget(filterEnabled)
        previewGroupLayout.addWidget(filterButton)

        #CameraGroup
        camGroup = QtWidgets.QGroupBox()
        layoutCg = QtWidgets.QVBoxLayout()
        camGroup.setLayout(layoutCg)

        layoutCg.addWidget(resetCameraButton)
        layoutCg.addWidget(dpGroup)
        layoutCg.addWidget(imageName)
        layoutCg.addWidget(saveVPtoImageButton)

        #import
        importGroup = QtWidgets.QGroupBox()
        importLayout = QtWidgets.QVBoxLayout()
        importGroup.setLayout(importLayout)

        importLayout.addWidget(importButton)
        importLayout.addWidget(importTextfield)

        #Insert Groups
        layoutRight.addWidget(progressBar)

        layoutLeft.addWidget(camGroup)
        layoutLeft.addWidget(previewGroup)

        paperCreationBox = QtWidgets.QGroupBox()
        paperCreationLayout = QtWidgets.QVBoxLayout()
        paperCreationBox.setLayout(paperCreationLayout)
        paperCreationLayout.addWidget(createPaperMeshButton)
        paperCreationLayout.addWidget(unfoldIterationsTextfield)
        paperCreationLayout.addWidget(importGroup)
        paperCreationLayout.addWidget(projectPerTriangle)
        paperCreationLayout.addWidget(resolutionWidth)

        editBox = QtWidgets.QGroupBox()
        editLayout = QtWidgets.QVBoxLayout()
        editBox.setLayout(editLayout)
        editLayout.addWidget(flattenButton)
        editLayout.addWidget(selectRegionButton)
        editLayout.addWidget(brightMultiplicationButton)
        editLayout.addWidget(finishButton)

        layoutLeft.addWidget(paperCreationBox)
        layoutLeft.addWidget(editBox)

        #debug
        debugBox = QtWidgets.QGroupBox()
        debugBox_Layout = QtWidgets.QVBoxLayout()
        debugBox.setLayout(debugBox_Layout)

        debugBox_Layout.addWidget(booleanButton)
        debugBox_Layout.addWidget(hierarchySlider)
        debugBox_Layout.addWidget(hierarchical_difference_button)
        debugBox_Layout.addWidget(testButton)

        layoutLeft.addWidget(debugBox)

        layoutRight.addWidget(addFileButton)

#       Append
        self.gridlayout.setColumnStretch(0, 1)
        self.gridlayout.setColumnStretch(1, 4)
        self.gridlayout.setColumnStretch(2, 1)

        self.gridlayout.addWidget(groupLeft,0,0)
        self.gridlayout.addWidget(groupRight,0,2)
        self.gridlayout.addWidget(self.vtkWidget,0,1)

        MainWindow.setCentralWidget(self.centralWidget)

        def addMesh(name,boolIsChild=False,layout=None):

            meshGroup = QtWidgets.QGroupBox()
            meshGroupLayout = QtWidgets.QVBoxLayout()
            meshGroup.setLayout(meshGroupLayout)

            #TODO exchange with regex
            label = QtWidgets.QLabel(name[-8:])

            colorBt = QtWidgets.QPushButton("Color")

            colorDialog = QtWidgets.QColorDialog()

            idx = self.filenames.index(name)

            opacity = QtWidgets.QLineEdit("0.75")

            self.inflateStruc.append(True)

            addMeshButton = QtWidgets.QPushButton("Add Mesh")

            def onOpacitySlider():
                if self.isfloat(opacity.text()):
                    org.changeOpacity(float(opacity.text()),idx)
                    self.vtkWidget.update()

            opacity.textChanged.connect(onOpacitySlider)

            def onColorSelect():
                red = colorDialog.currentColor().red()/255
                green = colorDialog.currentColor().green()/255
                blue = colorDialog.currentColor().blue()/255
                org.changeColor([red,green,blue],idx)
                self.vtkWidget.update()

            def onColorButton():
                colorDialog.hide()
                colorDialog.show()

            inflate = QtWidgets.QPushButton("Inflate")
            inflate.setCheckable(True)
            inflate.setChecked(True)
            clipping = QtWidgets.QPushButton("Clipping")
            clipping.setCheckable(True)

            def onInflate():
                clipping.setChecked(False)
                self.inflateStruc[idx] = True
            inflate.clicked.connect(onInflate)

            def onClipping():
                inflate.setChecked(False)
                self.inflateStruc[idx] = False
            clipping.clicked.connect(onClipping)

            colorDialog.currentColorChanged.connect(onColorSelect)
            #opacitySlider.sliderMoved.connect(onOpacitySlider)
            colorBt.clicked.connect(onColorButton)

            meshGroupLayout.addWidget(label)
            meshGroupLayout.addWidget(colorBt)
            meshGroupLayout.addWidget(opacity)
            meshGroupLayout.addWidget(colorDialog)
            meshGroupLayout.addWidget(inflate)
            meshGroupLayout.addWidget(clipping)

            def onAddMesh():

                dlg = QtWidgets.QFileDialog()
                dlg.setFileMode(QtWidgets.QFileDialog.ExistingFiles)

                filename = os.path.join(self.dirname, "../Meshes")
                dlg.setDirectory(filename)
                # dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
                # dlg.setFilter("*.stl")
                filenames = []

                if dlg.exec_():
                    for i in dlg.selectedFiles():
                        self.filenames.append(i)
                        addMesh(i,boolIsChild=True,layout=meshGroupLayout)
                        idx = self.filenames.index(i)
                        org.addActor(i,True)

                self.vtkWidget.GetRenderWindow().Render()

            addMeshButton.clicked.connect(onAddMesh)
            meshGroupLayout.addWidget(addMeshButton)

            if(boolIsChild==False):
                layoutRight.addWidget(meshGroup)
            else:
                layout.addWidget(meshGroup)

    def isfloat(self, value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    ##Class implementing the vtkInteractorStyleTrackballCamera
    # Adds a custom callback
    class MouseInteractorHighlightCell(vtk.vtkInteractorStyleTrackballCamera):

        def __init__(self, org, parent):
            self.AddObserver("RightButtonPressEvent", self.rightButtonPressEvent)
            self.LastPickedCell = None
            self.org = org
            self.parent = parent

        def rightButtonPressEvent(self, obj, event):
            clickPos = self.GetInteractor().GetEventPosition()

            picker = vtk.vtkCellPicker()
            picker.SetPickTextureData(False)

            pokedRen = self.GetInteractor().FindPokedRenderer(clickPos[0], clickPos[1])

            renId = (self.parent.vtkWidget.GetRenderWindow().GetRenderers().IsItemPresent(pokedRen)) - 2

            self.GetInteractor().GetRenderWindow()

            picker.Pick(clickPos[0], clickPos[1], 0, pokedRen)

            if (picker.GetCellId() != -1):

                self.NewPickedCell = picker.GetCellId()

                if(self.NewPickedCell != self.LastPickedCell):
                    self.parent.pickedIds[renId][self.parent.countOfPickedRegions].append(picker.GetCellId())
                    self.LastPickedCell = self.NewPickedCell
                    self.org.onGenerateColorMesh(renId,self.parent.pickedIds[renId])
                    #print(self.parent.pickedIds)

            self.OnRightButtonDown()
            return

class SimpleView(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.ren = vtk.vtkRenderer()
        org = organizer.Organizer(self.ren)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self,org,self.ren)
        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        style = self.ui.MouseInteractorHighlightCell(org, self.ui)
        style.SetDefaultRenderer(self.ren)
        #self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        self.iren.SetInteractorStyle(style)
        org.window = self.ui.vtkWidget.GetRenderWindow()
        org.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()

        self.ui.vtkWidget.GetRenderWindow().SetMultiSamples(0)
        self.ui.vtkWidget.GetRenderWindow().SetAlphaBitPlanes(True)
        #org.multiplyingActors()
        org.setUp()

def UpdateColorFilter(caller, ev):
    if UpdateColorFilter.run:
         UpdateColorFilter.sr.changeCameraForFilter(UpdateColorFilter.cam.GetPosition(), UpdateColorFilter.cam.GetFocalPoint(), UpdateColorFilter.cam.GetClippingRange(), UpdateColorFilter.cam.GetViewUp(), UpdateColorFilter.cam.GetDistance())
         UpdateColorFilter.sr.onMultiply(UpdateColorFilter.dp,UpdateColorFilter.filter,UpdateColorFilter.brighten)

if __name__ == '__main__':

    # Back up the reference to the exceptionhook
    sys._excepthook = sys.excepthook

    def my_exception_hook(exctype, value, traceback):
        # Print the error and traceback
        print(exctype, value, traceback)

        msgBox = QtWidgets.QMessageBox()
        msgBox.setText("Attribute not found, most likely an object or variable wasn't set.")
        msgBox.exec()

        # Call the normal Exception hook after
        sys._excepthook(exctype, value, traceback)
        #sys.exit(1)


    # Set the exception hook to our wrapping function
    sys.excepthook = my_exception_hook

    try:
        app = QtWidgets.QApplication(sys.argv)
        window = SimpleView()
        window.show()
        window.iren.Initialize()
        sys.exit(app.exec_())
    except:
        print("Exiting")


#    main()