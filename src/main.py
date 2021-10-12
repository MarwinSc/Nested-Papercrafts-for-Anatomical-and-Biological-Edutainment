import vtkmodules.all as vtk
import sys
from PyQt5 import QtWidgets, QtCore
import organizer
import os
from projectionStructure import ProjectionStructure
import random
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from viewpointselection.app import ViewPointComputation

class Ui_MainWindow(object):

    dpOcclusion = 0.1
    dpNOP = 10

    filenames = []

    dirname = os.path.dirname(__file__)

    """
    Initializing the UI elements.
    """
    def setupUi(self, MainWindow,org,ren):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1920, 1080)
        self.centralWidget = QtWidgets.QWidget(MainWindow)
        self.gridlayout = QtWidgets.QGridLayout(self.centralWidget)
        self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
        self.gridlayout.setSpacing(3)
        self.numberOfLoadedStructures = 0
        self.meshGroups = []

        ##deprecated (For flattening)
        self.pickedIds = [[[]],[[]],[[]]]
        self.countOfPickedRegions = 0

        #org.clearOutputDirectory()

        #      init Widgets
        resetCameraButton = QtWidgets.QPushButton("Reset Camera")

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
            org.colorFilterImage([red,green,blue])
            self.vtkWidget.update()

        filterButton.clicked.connect(onFilterColorButton)
        filterDialog.currentColorChanged.connect(onFilterColorChange)

        def onBrightMuliplication():
            org.brightenMultiplication()

        brightMultiplicationButton = QtWidgets.QPushButton("Bright Multiplication & Finish")
        brightMultiplicationButton.clicked.connect(onBrightMuliplication)

        brightenCheck = QtWidgets.QCheckBox("Brighten")
        brightenCheck.setChecked(True)

        addFileButton = QtWidgets.QPushButton("Add Mesh")

        def getFile():
            dlg = QtWidgets.QFileDialog()
            dlg.setFileMode(QtWidgets.QFileDialog.ExistingFiles)

            filename = os.path.join(self.dirname, "../Meshes")
            dlg.setDirectory(filename)
            dlg.setNameFilter("STL (*.stl)")

            if dlg.exec_():
                addMesh(dlg.selectedFiles())

            self.vtkWidget.GetRenderWindow().Render()

        addFileButton.clicked.connect(getFile)


        directImportPapermeshButton = QtWidgets.QPushButton("Direct Import Papermesh")

        def directImportPapermesh():
            dlg = QtWidgets.QFileDialog()
            dlg.setFileMode(QtWidgets.QFileDialog.ExistingFile)
            filename = os.path.join(self.dirname, "../Meshes")
            dlg.setDirectory(filename)
            dlg.setNameFilter("STL (*.stl)")

            if dlg.exec_():
                importPapermesh(dlg.selectedFiles()[0])

            self.vtkWidget.GetRenderWindow().Render()

        directImportPapermeshButton.clicked.connect(directImportPapermesh)


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
        unfoldIterationsTextfield.setText("50000")

        saveVPtoImageButton = QtWidgets.QPushButton("Save as Image")
        saveVPtoImageButton.clicked.connect(onSaveVpToImage)

        def onUnfoldPaperMesh():
            try:
                iterations = int(unfoldIterationsTextfield.text())
            except:
                iterations = 50000
            org.unfoldPaperMeshPass(iterations)
            self.vtkWidget.update()

        unfoldPaperMeshButton = QtWidgets.QPushButton("Unfold Papermesh")
        unfoldPaperMeshButton.clicked.connect(onUnfoldPaperMesh)

        resolutionWidth = QtWidgets.QLineEdit()
        resolutionWidth.setText("500")

        def onProjectPerTriangle():
            try:
                width = int(resolutionWidth.text())
            except:
                width = 500
            org.project(resolution = [width, width])
            #org.projectPassTemp()
            if renderPaperMeshesButton.isChecked():
                org.hierarchical_mesh_anchor.renderPaperMeshes(ren)
            self.vtkWidget.update()
#            self.centralWidget.update()

        projectPerTriangle = QtWidgets.QPushButton("Project")
        projectPerTriangle.clicked.connect(onProjectPerTriangle)

        importTextfield = QtWidgets.QLineEdit()
        importTextfield.setText("model")

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

        testButton = QtWidgets.QPushButton("Temp")
        def onUnfoldTest():
            org.writeObjMeshes()

        testButton.clicked.connect(onUnfoldTest)

        treeToStringButton = QtWidgets.QPushButton("hierarchy toString")
        def onTreeToString():
            org.hierarchical_mesh_anchor.toString()
        treeToStringButton.clicked.connect(onTreeToString)

        renderStructuresButton = QtWidgets.QPushButton("Render Structures")
        renderStructuresButton.setCheckable(True)
        renderStructuresButton.setChecked(True)
        renderPaperMeshesButton = QtWidgets.QPushButton("Render PaperMeshes")
        renderPaperMeshesButton.setCheckable(True)
        renderPaperMeshesButton.setChecked(True)

        def onRenderStructures():
            if not renderStructuresButton.isChecked():
                ren.RemoveAllViewProps()
                if renderPaperMeshesButton.isChecked():
                    org.hierarchical_mesh_anchor.renderPaperMeshes(ren)
            else:
                org.hierarchical_mesh_anchor.renderStructures(ren)
            self.vtkWidget.GetRenderWindow().Render()

        def onRenderPaperMeshes():
            if not renderPaperMeshesButton.isChecked():
                ren.RemoveAllViewProps()
                if renderStructuresButton.isChecked():
                    org.hierarchical_mesh_anchor.renderStructures(ren)
            else:
                org.hierarchical_mesh_anchor.renderPaperMeshes(ren)
            self.vtkWidget.GetRenderWindow().Render()

        renderStructuresButton.clicked.connect(onRenderStructures)
        renderPaperMeshesButton.clicked.connect(onRenderPaperMeshes)


        entropyButton = QtWidgets.QPushButton("Cut")
        def onEntropy():

            hl = org.hierarchical_mesh_anchor.toActorList()
            window = ViewPointComputation(hl)
            viewpoints, bounds = window.getMaxViewPoints_AndBounds()
            org.renderCutPlanes(viewpoints,bounds)
            org.cutHM(viewpoints, bounds)

            self.vtkWidget.update()

        entropyButton.clicked.connect(onEntropy)


        writeUnfoldedButton = QtWidgets.QPushButton("Write Unfolded Meshes")
        def onTestCut():
            org.writeUnfolded()
        writeUnfoldedButton.clicked.connect(onTestCut)

        importUnfoldedButton = QtWidgets.QPushButton("Import Unfolded Meshes")
        def onImportUnfolded():
            org.importUnfolded()
        importUnfoldedButton.clicked.connect(onImportUnfolded)

        """
        Adding UI elements to their respective container.
        """
        layoutRight = QtWidgets.QFormLayout()
        groupBox = QtWidgets.QGroupBox()
        groupBox.setLayout(layoutRight)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(groupBox)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        scroll.setMinimumWidth(350)

        vBoxLayout = QtWidgets.QVBoxLayout()
        vBoxLayout.addWidget(scroll)

        groupRight = QtWidgets.QWidget()
        groupRight.setLayout(vBoxLayout)

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

        checkboxGroup = QtWidgets.QGroupBox()
        checkboxLayout = QtWidgets.QHBoxLayout()
        checkboxGroup.setLayout(checkboxLayout)
        checkboxLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        checkboxLayout.addWidget(brightenCheck)
        checkboxLayout.addWidget(multiplyingCheck)
        checkboxLayout.addWidget(filterEnabled)

        previewGroupLayout.addWidget(multiplyingButton)
        previewGroupLayout.addWidget(checkboxGroup)
        previewGroupLayout.addWidget(filterButton)
        previewGroupLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)

        #ViewGroup
        camGroup = QtWidgets.QGroupBox()
        layoutCg = QtWidgets.QVBoxLayout()
        camGroup.setLayout(layoutCg)

        layoutCg.addWidget(resetCameraButton)
        layoutCg.addWidget(dpGroup)
        layoutCg.addWidget(imageName)
        layoutCg.addWidget(saveVPtoImageButton)

        renderBox = QtWidgets.QGroupBox()
        renderBoxLayout = QtWidgets.QHBoxLayout()
        renderBox.setLayout(renderBoxLayout)
        #renderBoxLayout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
        renderBoxLayout.addWidget(renderStructuresButton)
        renderBoxLayout.addWidget(renderPaperMeshesButton)

        layoutCg.addWidget(renderBox)

        #import


        paperCreationBox = QtWidgets.QGroupBox()
        paperCreationLayout = QtWidgets.QVBoxLayout()
        paperCreationBox.setLayout(paperCreationLayout)
        paperCreationLayout.addWidget(entropyButton)
        paperCreationLayout.addWidget(unfoldPaperMeshButton)
        paperCreationLayout.addWidget(unfoldIterationsTextfield)
        paperCreationLayout.addWidget(projectPerTriangle)
        paperCreationLayout.addWidget(resolutionWidth)
        paperCreationLayout.addWidget(brightMultiplicationButton)

        #Insert Groups
        layoutLeft.addWidget(camGroup)
        layoutLeft.addWidget(previewGroup)
        layoutLeft.addWidget(paperCreationBox)
        #layoutLeft.addWidget(editBox)

        #debug
        debugBox = QtWidgets.QGroupBox()
        debugBox_Layout = QtWidgets.QVBoxLayout()
        debugBox.setLayout(debugBox_Layout)

        #debugBox_Layout.addWidget(booleanButton)
        debugBox_Layout.addWidget(hierarchySlider)
        debugBox_Layout.addWidget(hierarchical_difference_button)
        debugBox_Layout.addWidget(testButton)
        debugBox_Layout.addWidget(treeToStringButton)
        debugBox_Layout.addWidget(writeUnfoldedButton)
        debugBox_Layout.addWidget(importUnfoldedButton)
        #debugBox_Layout.addWidget(importAnchorPapermeshButton)

        layoutLeft.addWidget(debugBox)

        layoutRight.addWidget(addFileButton)
        layoutRight.addWidget(directImportPapermeshButton)


#       Append
        self.gridlayout.setColumnStretch(0, 1)
        self.gridlayout.setColumnStretch(1, 4)
        self.gridlayout.setColumnStretch(2, 1)

        self.gridlayout.addWidget(groupLeft,0,0)
        self.gridlayout.addWidget(groupRight,0,2)
        self.gridlayout.addWidget(self.vtkWidget,0,1)

        MainWindow.setCentralWidget(self.centralWidget)

        def importPapermesh(name):
            mesh = ProjectionStructure(name, self.numberOfLoadedStructures)
            mesh.hierarchicalMesh = org.directImportPapermesh(mesh)
            mesh.hierarchicalMesh.label = QtWidgets.QLabel("Not Unfolded")

            meshGroupBox = QtWidgets.QGroupBox()
            meshGroupLayout = QtWidgets.QVBoxLayout()
            meshGroupBox.setLayout(meshGroupLayout)
            borderColor = '#%02X%02X%02X' % (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
            meshGroupBox.setObjectName("meshGroup")
            meshGroupBox.setStyleSheet("QWidget#meshGroup" + " {border: 1px solid " + borderColor + ";}")

            meshBox, projectBox = setupMeshUiElements(mesh, mesh.filename)
            meshGroupLayout.addWidget(meshBox)
            meshGroupLayout.addWidget(projectBox)
            meshGroupLayout.addWidget(mesh.hierarchicalMesh.label)

            mesh.hierarchicalMesh.setUpAdditionalUiElements()
            meshGroupLayout.addWidget(mesh.hierarchicalMesh.ui_elements_box)

            layoutRight.addWidget(meshGroupBox)

            ren.RemoveAllViewProps()
            if renderStructuresButton.isChecked():
                ren.RemoveAllViewProps()
                org.hierarchical_mesh_anchor.renderStructures(ren)
            if renderPaperMeshesButton.isChecked():
                org.hierarchical_mesh_anchor.renderPaperMeshes(ren)

            self.meshGroups.append(meshGroupLayout)

        def addMesh(names):

            meshes = []
            for i,name in enumerate(names):
                self.numberOfLoadedStructures += 1
                mesh = ProjectionStructure(name, self.numberOfLoadedStructures)
                meshes.append(mesh)

            hierarchicalMesh = org.addMesh(meshes)
            hierarchicalMesh.label = QtWidgets.QLabel("Not Unfolded")

            meshGroupBox = QtWidgets.QGroupBox()
            meshGroupLayout = QtWidgets.QVBoxLayout()
            meshGroupBox.setLayout(meshGroupLayout)

            for mesh in meshes:
                mesh.hierarchicalMesh = hierarchicalMesh
                meshBox, projectBox = setupMeshUiElements(mesh,mesh.filename)

                meshGroupLayout.addWidget(meshBox)
                meshGroupLayout.addWidget(projectBox)

            meshGroupLayout.addWidget(hierarchicalMesh.label)

            hierarchicalMesh.setUpAdditionalUiElements()
            meshGroupLayout.addWidget(hierarchicalMesh.ui_elements_box)

            borderColor = '#%02X%02X%02X' % (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
            meshGroupBox.setObjectName("meshGroup")#{}".format(level))
            meshGroupBox.setStyleSheet("QWidget#meshGroup" + " {border: 1px solid " + borderColor + ";}")
            layoutRight.addWidget(meshGroupBox)

            ren.RemoveAllViewProps()
            if renderStructuresButton.isChecked():
                ren.RemoveAllViewProps()
                org.hierarchical_mesh_anchor.renderStructures(ren)
            if renderPaperMeshesButton.isChecked():
                org.hierarchical_mesh_anchor.renderPaperMeshes(ren)

        def setupMeshUiElements(mesh,name):

            label = QtWidgets.QLabel(name.split("/")[-1])
            colorBt = QtWidgets.QPushButton("Color")
            colorDialog = QtWidgets.QColorDialog()
            opacity = QtWidgets.QLineEdit("1.0")

            def onOpacityChange():
                if self.isfloat(opacity.text()):
                    mesh.setOpacity(float(opacity.text()))
                    self.vtkWidget.update()

            opacity.textChanged.connect(onOpacityChange)

            def onColorSelect():
                red = colorDialog.currentColor().red()
                green = colorDialog.currentColor().green()
                blue = colorDialog.currentColor().blue()
                mesh.setColor([red/ 255, green/ 255, blue/ 255])
                #borderColor = '#%02X%02X%02X' % (red,green,blue)
                #meshGroupBox.setStyleSheet("QWidget#meshGroup{}".format(level) + " {border: 1px solid " + borderColor + ";}")
                self.vtkWidget.update()

            def onColorButton():
                colorDialog.hide()
                colorDialog.show()

            inflate = QtWidgets.QPushButton("Inflate")
            inflate.setCheckable(True)
            inflate.setChecked(True)
            clipping = QtWidgets.QPushButton("Clipping")
            clipping.setCheckable(True)
            cube = QtWidgets.QPushButton("Cube")
            cube.setCheckable(True)

            def onInflate():
                clipping.setChecked(False)
                cube.setChecked(False)
                mesh.projectionMethod = mesh.ProjectionMethod.Inflate

            inflate.clicked.connect(onInflate)

            def onClipping():
                inflate.setChecked(False)
                cube.setChecked(False)
                mesh.projectionMethod = mesh.ProjectionMethod.Clipping

            clipping.clicked.connect(onClipping)

            def onCube():
                inflate.setChecked(False)
                clipping.setChecked(False)
                mesh.projectionMethod = mesh.ProjectionMethod.Cube

            cube.clicked.connect(onCube)

            colorDialog.currentColorChanged.connect(onColorSelect)
            colorBt.clicked.connect(onColorButton)

            meshBox = QtWidgets.QGroupBox()
            meshLayout = QtWidgets.QHBoxLayout()
            meshBox.setLayout(meshLayout)
            projectBox = QtWidgets.QGroupBox()
            projectLayout = QtWidgets.QHBoxLayout()
            projectBox.setLayout(projectLayout)

            meshLayout.addWidget(label)
            meshLayout.addWidget(colorBt)
            meshLayout.addWidget(opacity)
            meshLayout.minimumSize()
            projectLayout.addWidget(inflate)
            projectLayout.addWidget(clipping)
            projectLayout.addWidget(cube)
            projectLayout.minimumSize()

            #meshGroupLayout.addWidget(meshBox)
            #meshGroupLayout.addWidget(projectBox)
            return meshBox, projectBox

    def isfloat(self, value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    '''
    class MouseInteractorHighlightCell(vtk.vtkInteractorStyleTrackballCamera):
        
        Previously handled the manual flattening of the projection meshes, now obsolete.
        
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
    '''
class SimpleView(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.ren = vtk.vtkRenderer()
        org = organizer.Organizer(self.ren)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self,org,self.ren)
        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        #style = self.ui.MouseInteractorHighlightCell(org, self.ui)
        #style.SetDefaultRenderer(self.ren)
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        #self.iren.SetInteractorStyle(style)
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
        msgBox.setText(str(value))
        msgBox.exec()

        # Call the normal Exception hook after
        sys._excepthook(exctype, value, traceback)

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