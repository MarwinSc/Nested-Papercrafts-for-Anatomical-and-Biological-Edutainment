import vtkmodules.all as vtk
import sys
from PyQt5 import QtWidgets, QtCore
import organizer
import os
from mesh import Mesh
import random

from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

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
        ##For flattening.
        self.pickedIds = [[[]],[[]],[[]]]
        self.countOfPickedRegions = 0

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
            #dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
            #dlg.setFilter("*.stl")
            filenames = []

            if dlg.exec_():
                self.filenames.append(dlg.selectedFiles()[0])
                layout = addMesh(dlg.selectedFiles()[0], newGroup=True)
                for i in range(1,len(dlg.selectedFiles())):
                    self.filenames.append(dlg.selectedFiles()[i])
                    addMesh(dlg.selectedFiles()[i], layout=layout,newGroup=False)
                # delete the button now
                addFileButton.deleteLater()

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

        def onUnfoldPaperMesh():
            try:
                iterations = int(unfoldIterationsTextfield.text())
            except:
                iterations = 100
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
            org.projectPass(resolution = [width, width])
            #org.projectPassTemp()
            self.vtkWidget.update()
#            self.centralWidget.update()

        projectPerTriangle = QtWidgets.QPushButton("Project")
        projectPerTriangle.clicked.connect(onProjectPerTriangle)

        importTextfield = QtWidgets.QLineEdit()
        importTextfield.setText("model")

        def onImport():
            org.importUnfoldedMeshPass(importTextfield.text())
            self.vtkWidget.GetRenderWindow().Render()
        importButton = QtWidgets.QPushButton("Import .obj Papermesh")
        importButton.clicked.connect(onImport)

        '''
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

        def onFinish():
            org.finish()
        finishButton = QtWidgets.QPushButton("Finish")
        finishButton.clicked.connect(onFinish)
        '''

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
        importGroup = QtWidgets.QGroupBox()
        importLayout = QtWidgets.QVBoxLayout()
        importGroup.setLayout(importLayout)

        importLayout.addWidget(importButton)
        importLayout.addWidget(importTextfield)

        paperCreationBox = QtWidgets.QGroupBox()
        paperCreationLayout = QtWidgets.QVBoxLayout()
        paperCreationBox.setLayout(paperCreationLayout)
        paperCreationLayout.addWidget(unfoldPaperMeshButton)
        paperCreationLayout.addWidget(unfoldIterationsTextfield)
        paperCreationLayout.addWidget(importGroup)
        paperCreationLayout.addWidget(projectPerTriangle)
        paperCreationLayout.addWidget(resolutionWidth)
        paperCreationLayout.addWidget(brightMultiplicationButton)

        #editBox = QtWidgets.QGroupBox()
        #editLayout = QtWidgets.QVBoxLayout()
        #editBox.setLayout(editLayout)
        #editLayout.addWidget(flattenButton)
        #editLayout.addWidget(selectRegionButton)
        #editLayout.addWidget(finishButton)

        #Insert Groups
        layoutLeft.addWidget(camGroup)
        layoutLeft.addWidget(previewGroup)
        layoutLeft.addWidget(paperCreationBox)
        #layoutLeft.addWidget(editBox)

        #debug
        debugBox = QtWidgets.QGroupBox()
        debugBox_Layout = QtWidgets.QVBoxLayout()
        debugBox.setLayout(debugBox_Layout)

        debugBox_Layout.addWidget(booleanButton)
        debugBox_Layout.addWidget(hierarchySlider)
        debugBox_Layout.addWidget(hierarchical_difference_button)
        debugBox_Layout.addWidget(testButton)
        debugBox_Layout.addWidget(treeToStringButton)

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

        def addMesh(name, layout=None, level=0, parent = None, childId = -1, newGroup = True):
            '''
            Method forwarding the loaded structure as mesh to the organizer class.
            Handling the creation of UI elements and event calls for the loaded files.
            Gets called recursively, either to create a new hierarchical mesh child, or to append a mesh to an existing hierarchical mesh.
            :param name: The filename.
            :param layout: The layout in which the new UI elements are inserted, if newGroup = False.
            :param level: The level of the loaded structures.
            :param parent: The parent Hierarchical Mesh of the new loaded structures.
            :param childId: The ID of the child to which the loaded meshes should be appended to.
            :param newGroup: Boolean indicating if a new Box should be created, to hold the created UI elements or if the new elements are inserted in an existing Box.
            :return: If newGroup is True, the created Box is returned.
            '''
            idx = self.filenames.index(name)
            r = lambda: random.randint(0, 255)
            levelColor = '#%02X%02X%02X' % (r(), r(), r())

            mesh = Mesh(name,idx)

            if level > 0:
                hierarchicalMesh = org.addMesh(mesh,parent,childId)
            else:
                hierarchicalMesh = org.addMesh(mesh, None, childId)
            mesh.hierarchicalMesh = hierarchicalMesh

            ren.RemoveAllViewProps()
            if renderStructuresButton.isChecked():
                ren.RemoveAllViewProps()
                org.hierarchical_mesh_anchor.renderStructures(ren)
            if renderPaperMeshesButton.isChecked():
                org.hierarchical_mesh_anchor.renderPaperMeshes(ren)

            label = QtWidgets.QLabel(name.split("/")[-1])

            colorBt = QtWidgets.QPushButton("Color")
            colorDialog = QtWidgets.QColorDialog()
            opacity = QtWidgets.QLineEdit("0.5")

            def onOpacityChange():
                if self.isfloat(opacity.text()):
                    mesh.setOpacity(float(opacity.text()))
                    self.vtkWidget.update()

            opacity.textChanged.connect(onOpacityChange)

            def onColorSelect():
                red = colorDialog.currentColor().red()/255
                green = colorDialog.currentColor().green()/255
                blue = colorDialog.currentColor().blue()/255
                mesh.setColor([red,green,blue])
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

            if newGroup:
                meshGroupBox = QtWidgets.QGroupBox()
                meshGroupBox.setObjectName("meshGroup{}".format(level))
                meshGroupBox.setStyleSheet("QWidget#meshGroup{}".format(level)+" {border: 1px solid " + levelColor + ";}")
                meshGroupLayout = QtWidgets.QVBoxLayout()
                meshGroupBox.setLayout(meshGroupLayout)

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

            addSubmeshButton = QtWidgets.QPushButton("Add Submesh")
            addMeshButton = QtWidgets.QPushButton("Add Mesh")

            def onAddSubMesh():
                '''
                Recursive call of the addMesh method for a new sub mesh.
                :return:
                '''
                dlg = QtWidgets.QFileDialog()
                dlg.setFileMode(QtWidgets.QFileDialog.ExistingFiles)

                filename = os.path.join(self.dirname, "../Meshes")
                dlg.setDirectory(filename)

                if dlg.exec_():
                    self.filenames.append(dlg.selectedFiles()[0])
                    layout = addMesh(dlg.selectedFiles()[0], layout=meshGroupLayout, level=level + 1, parent=hierarchicalMesh, childId=childId + 1,newGroup=True)

                    for i in range(1, len(dlg.selectedFiles())):
                        self.filenames.append(dlg.selectedFiles()[i])
                        addMesh(dlg.selectedFiles()[i],layout=layout,level = level+1,parent=hierarchicalMesh,childId=childId+1,newGroup=False)


                self.vtkWidget.GetRenderWindow().Render()

            addSubmeshButton.clicked.connect(onAddSubMesh)

            def onAddMesh():
                '''
                Recursive call of the addMesh method to append a new mesh to an hierarchical mesh.
                :return:
                '''
                dlg = QtWidgets.QFileDialog()
                dlg.setFileMode(QtWidgets.QFileDialog.ExistingFiles)

                filename = os.path.join(self.dirname, "../Meshes")
                dlg.setDirectory(filename)

                if dlg.exec_():
                    for i in dlg.selectedFiles():
                        self.filenames.append(i)
                        addMesh(i, layout=meshGroupLayout, level=level, parent=parent, childId=childId, newGroup=False)

                self.vtkWidget.GetRenderWindow().Render()

            addMeshButton.clicked.connect(onAddMesh)

            if newGroup:
                meshGroupLayout.addWidget(addMeshButton)
                meshGroupLayout.addWidget(addSubmeshButton)
                meshGroupLayout.addWidget(meshBox)
                meshGroupLayout.addWidget(projectBox)
                if layout:
                    layout.addWidget(meshGroupBox)
                else:
                    layoutRight.addWidget(meshGroupBox)
                return meshGroupLayout
            else:
                layout.addWidget(meshBox)
                layout.addWidget(projectBox)

    def isfloat(self, value):
        try:
            float(value)
            return True
        except ValueError:
            return False


    class MouseInteractorHighlightCell(vtk.vtkInteractorStyleTrackballCamera):
        '''
        Previously handled the manual flattening of the projection meshes, now obsolete.
        '''
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