import vtkmodules.all as vtk
import os
from meshProcessing import MeshProcessing
from projector import Projector
from imageProcessing import ImageProcessor
import util
from mu3d.mu3dpy.mu3d import Graph
from src.hierarchicalMesh import HierarchicalMesh
from PyQt5 import QtWidgets
import time

class Organizer():
    '''
    Class responsible for Rendering Tasks, I/O and forwarding method calls to ImageProcessing and MeshProcessing.
    '''
    def __init__(self,ren):
        self.ren = ren

    dirname = os.path.dirname(__file__)

    actorList = []
    renList = []

    height = 2000
    width = 2000

    meshProcessor = MeshProcessing()
    projector = Projector()
    imageProcessor = ImageProcessor(height,width)


    camera = vtk.vtkCamera()

    depthPeeling = True
    occlusion = 0.1
    numberOfPeels = 10
    filter = False
    sessionMultiplySaves = 0

    fullViewport = [0.0, 0.0, 1.0, 1.0]
    noViewport = [0.0, 0.0, 0.0, 0.0]

    resultRen = vtk.vtkRenderer()
    resultRen.SetBackground(255.0, 255.0, 255.0)
    resultRen.ResetCamera()
    resultRen.InteractiveOff()

    ctf = vtk.vtkColorTransferFunction()
    ctf.SetColorSpaceToHSV()
    ctf.HSVWrapOn()
    ctf.AddHSVPoint(0.0, 0.0, 1.0, 1.0)
    ctf.AddHSVPoint(25.0, 0.25, 1.0, 1.0)
    ctf.AddHSVPoint(50.0, 0.5, 1.0, 1.0)
    ctf.AddHSVPoint(75.0, 0.75, 1.0, 1.0)
    ctf.AddHSVPoint(100.0, 1.0, 1.0, 1.0)

    hierarchical_mesh_anchor = HierarchicalMesh(None,None,meshProcessor)

    def setUp(self):
        '''
        initial settings for the renderer and the camera.
        :return:
        '''
        # Camera
        self.camera.SetViewUp(0, 0, 1)
        self.camera.SetPosition(0, -600, 0)
        self.camera.SetFocalPoint(0, 0, 0)
        self.camera.ComputeViewPlaneNormal()

        self.ren.SetBackground(255.0, 255.0, 255.0)
        self.ren.SetActiveCamera(self.camera)
        self.camera.Dolly(1.5)

        self.ren.ResetCameraClippingRange()
        self.ren.ResetCamera()

    def resetCamera(self):
        '''
        reset the camera.
        :return:
        '''
        self.camera.SetViewUp(0, 0, 1)
        self.camera.SetPosition(0, -600, 0)
        self.camera.SetFocalPoint(0, 0, 0)
        self.camera.ComputeViewPlaneNormal()
        self.camera.SetScreenBottomLeft(0.1,0.1,0.1)
        self.camera.SetScreenTopRight(0.2,0.2,0.2)
        self.camera.SetScreenBottomRight(0.3,0.3,0.3)
        self.ren.ResetCameraClippingRange()
        self.ren.ResetCamera()

    def onMultiply(self,depthPeeling,filter,brighten):
        '''
        Forwarding the rendering and color multiplication of the loaded structures.
        :param depthPeeling:
        :param filter:
        :param brighten:
        :return:
        '''
        imageActor = vtk.vtkImageActor()

        #firsttime multiplication
        if not self.window.HasRenderer(self.resultRen):
            self.window.AddRenderer(self.resultRen)

        self.resultRen.RemoveActor(imageActor)
        self.ren.SetViewport(self.noViewport)
        self.resultRen.SetViewport(self.fullViewport)
        result = self.imageProcessor.multiplyingActors(depthPeeling,filter,brighten,self.hierarchical_mesh_anchor.getAllMeshes(),self.camera,self.height,self.width,self.occlusion,self.numberOfPeels)

        imageActor = vtk.vtkImageActor()
        imageActor.GetMapper().SetInputData(result)
        imageActor.Update()

        self.resultRen.AddActor(imageActor)
        self.resultRen.Render()

        filename = os.path.join(self.dirname, "../out/2D/multiply{}.png".format(self.sessionMultiplySaves))
        util.writeImage(result,filename)

        self.sessionMultiplySaves += 1

    def brightenMultiplication(self):
        '''
        Multiplies the created unfolding images of the structures into a single unfolded texture.
        :return:
        '''
        reader = vtk.vtkPNGReader()
        imgList=[]

        for i in range(len(self.hierarchical_mesh_anchor.getAllMeshes())):

            filename = os.path.join(self.dirname, "../out/2D/unfolding{}.png".format(i))

            reader.SetFileName(filename)

            castFilter = vtk.vtkImageCast()
            castFilter.SetInputConnection(reader.GetOutputPort())
            castFilter.SetOutputScalarTypeToUnsignedChar()
            castFilter.Update()

            imageDim = castFilter.GetOutput().GetDimensions()

            width = imageDim[0]
            height = imageDim[1]

            imgList.append(self.imageProcessor.optimizedBrighten(castFilter.GetOutput(), width, height, str(i)))

            if i == 1:
                result = self.imageProcessor.normalizeMultiplication(imgList[i-1], imgList[i], width, height).GetOutput()
            elif(i > 1):
                resultCast = self.imageProcessor.normalizeMultiplication(result, imgList[i], width, height)

        writer = vtk.vtkPNGWriter()
        filename = os.path.join(self.dirname, "../out/2D/texture.png")

        writer.SetFileName(filename)
        writer.SetInputConnection(resultCast.GetOutputPort())
        writer.Write()
        self.finish()
    '''
    def addMesh(self, mesh, parent, childId):
        
        Method that manages adding a mesh to the hierarchy by either:
        Create a hierarchical mesh and initialize the tree,
        create a hierarchical mesh and add as child,
        append a mesh to level 0,
        append a mesh to a child.

        :param mesh: the new mesh.
        :param parent: if given, the parent mesh.
        :param childId: the list id of the child to which the mesh is appended.
        :return: the hierarchical mesh that was created/updated.
        
        self.ren.RemoveAllViewProps()

        if hasattr(self,"hierarchical_mesh_anchor"):
            if parent:
                # append as new child
                if childId < 0 or len(parent.children) == 0:
                    hierarchicalMesh = HierarchicalMesh(parent,mesh,self.meshProcessor)
                    parent.children.append(hierarchicalMesh)
                    self.hierarchical_mesh_anchor.renderStructures(self.ren)
                    self.hierarchical_mesh_anchor.renderPaperMeshes(self.ren)
                    return hierarchicalMesh
                # append to the meshes of a child
                else:
                    parent.children[childId].appendMesh(mesh)
                    self.hierarchical_mesh_anchor.renderStructures(self.ren)
                    self.hierarchical_mesh_anchor.renderPaperMeshes(self.ren)
                    return parent.children[childId]
            # add at level 0
            else:
                self.hierarchical_mesh_anchor.appendMesh(mesh)
                self.hierarchical_mesh_anchor.renderStructures(self.ren)
                self.hierarchical_mesh_anchor.renderPaperMeshes(self.ren)
                return self.hierarchical_mesh_anchor
        # initialize Tree
        else:
            self.hierarchical_mesh_anchor = HierarchicalMesh(None,mesh,self.meshProcessor)
            self.hierarchical_mesh_anchor.renderStructures(self.ren)
            self.hierarchical_mesh_anchor.renderPaperMeshes(self.ren)
            return self.hierarchical_mesh_anchor
    '''
    def addMesh(self, meshes):
        newHierarchicalMesh = HierarchicalMesh(None,meshes,self.meshProcessor)
        self.hierarchical_mesh_anchor.add(newHierarchicalMesh)
        return newHierarchicalMesh

    def directImportPapermesh(self, mesh):
        hm = HierarchicalMesh(None,None,self.meshProcessor)
        hm.papermesh = mesh.getActor().GetMapper().GetInput()
        hm.mesh = mesh.getActor()
        hm.meshes.append(mesh)
        hm.setName(hm.writePapermeshStlAndOff("Temp"))
        self.hierarchical_mesh_anchor.add(hm)
        return hm

    def draw_level(self, level):
        self.hierarchical_mesh_anchor.render(level, self.ren)

    def hierarchical_difference(self):
        self.hierarchical_mesh_anchor.recursive_difference()

    def colorFilterImage(self,color):
        self.imageProcessor.canvas_source.SetDrawColor(color[0],color[1],color[2],255)
        self.imageProcessor.canvas_source.FillBox(0, self.width-1, 0, self.height-1)
        self.imageProcessor.canvas_source.Update()

    def changeCameraForFilter(self,pos,fop,clr,vup,dis):
        self.camera.SetPosition(pos)
        self.camera.SetFocalPoint(fop)
        self.camera.SetClippingRange(clr)
        self.camera.SetViewUp(vup)
        self.camera.SetDistance(dis)

    def swapViewports(self, resultRenMaxView):
        '''
        Swaps between the normal renderer self. ren and self.resultRen.
        :param resultRenMaxView:
        :return:
        '''
        if resultRenMaxView:
            self.ren.SetViewport(self.noViewport)
            self.resultRen.SetViewport(self.fullViewport)
        else:
            self.ren.SetViewport(self.fullViewport)
            self.resultRen.SetViewport(self.noViewport)

    def unfoldPaperMeshPass(self, iterations):
        '''
        Forwards the unfold call to the hierarchical tree.
        :param iterations: Iterations of the mu3d unfolding.
        :return:
        '''
        self.hierarchical_mesh_anchor.unfoldWholeHierarchy(iterations)

    '''
    def importUnfoldedMeshPass(self, name):
        #deprecated
        self.hierarchical_mesh_anchor.unfoldedActor = self.meshProcessor.importUnfoldedMesh(name)
        self.meshProcessor.createDedicatedMeshes(self.hierarchical_mesh_anchor)

    def importPapermeshAnchor(self):
        #deprecated
        filename = os.path.join(self.dirname, "../out/3D/papermesh.obj")
        self.hierarchical_mesh_anchor.papermesh = util.readObj(filename)
    '''

    def project(self, hierarchy, resolution):
        '''
        Calls the projectPerTriangle() and createUnfoldedPaperMesh() from the projector class.
        :param hierarchy:
        :param resolution:
        :return:
        '''
        resultMeshes = []
        idx = 0
        meshes = hierarchy.getAllMeshes(asActor=False)
        for a in meshes:
            try:
                #todo projection for whole hierarchy not just level one child one
                if a.hierarchicalMesh.getLevel() > 1: raise Exception("Projection for nested meshes not implemented")
                mesh = (self.projector.projectPerTriangle(a.projectionActor, a.getActor(), idx, resolution))
                resultMeshes.append(mesh)
                self.projector.createUnfoldedPaperMesh(mesh, a.hierarchicalMesh.unfoldedActor, idx)
            except Exception as e:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setText(str(e))
                msgBox.exec()
            idx += 1
        return resultMeshes

    def projectPass(self,resolution = [500,500]):
        #todo projection for whole hierarchy not just level one child one
        self.ren.SetViewport([0.0, 0.0, 0.0, 0.0])
        hm = self.hierarchical_mesh_anchor.children[0]

        actors = self.project(hm, resolution)
        count = 0

        filename = os.path.join(self.dirname, "../out/2D/unfolding{}.png".format(count))
        readerFac = vtk.vtkImageReader2Factory()
        imageReader = readerFac.CreateImageReader2(filename)
        imageReader.SetFileName(filename)
        texture = vtk.vtkTexture()
        texture.SetInputConnection(imageReader.GetOutputPort())



        hm.unfoldedActor.SetTexture(texture)
        hm.unfoldedActor.Modified()

        util.writeObj(hm.unfoldedActor.GetMapper().GetInput(), "tempPaper")

        self.renderers = self.setUpResultRenderers(self.camera,len(actors)+1)
        self.ren.GetRenderWindow().AddRenderer(self.renderers[0])
        self.renderers[0].AddActor(hm.unfoldedActor)

        for actor in actors:
            filename = os.path.join(self.dirname, "../out/2D/texture/texture{}.png".format(count))
            readerFac = vtk.vtkImageReader2Factory()
            imageReader = readerFac.CreateImageReader2(filename)
            imageReader.SetFileName(filename)
            texture = vtk.vtkTexture()
            texture.SetInputConnection(imageReader.GetOutputPort())

            self.ren.GetRenderWindow().AddRenderer(self.renderers[count+1])
            actors[count].GetProperty().SetColor([1.0,1.0,1.0])
            actors[count].SetTexture(texture)
            self.renderers[count+1].AddActor(actors[count])
            count += 1

        actors.insert(0,hm.unfoldedActor)
        self.dedicatedPaperMeshes = actors
        return self.renderers

    def finish(self):
        '''
        #todo whole hierarchy not just level 1 child 1
        Loads the final multiplied texture and assigns it to the papermesh.
        :return:
        '''
        if hasattr(self, "renderers"):
            for ren in self.renderers:
                ren.SetViewport(self.noViewport)
        self.ren.SetViewport(self.fullViewport)

        readerFac = vtk.vtkImageReader2Factory()
        filename = os.path.join(self.dirname,"../out/2D/texture.png")
        imageReader = readerFac.CreateImageReader2(filename)
        imageReader.SetFileName(filename)
        texture = vtk.vtkTexture()
        texture.SetInputConnection(imageReader.GetOutputPort())

        hm = self.hierarchical_mesh_anchor.children[0]
        hm.unfoldedActor.SetTexture(texture)

        self.ren.RemoveAllViewProps()
        hm.unfoldedActor.GetProperty().SetOpacity(1.0)
        self.ren.AddActor(hm.unfoldedActor)

    def setUpResultRenderers(self, camera, maxNofChilds):
        '''
        sets up the horizontally concatenated renderers for each projected structure.
        #todo, set up renderers for nested structures below level zero.
        :param camera:
        :param maxNofChilds:
        :return:
        '''
        renderers = []
        width = 1.0 / maxNofChilds
        right = 0.0
        left = width

        for i in range(maxNofChilds):

            render = vtk.vtkRenderer()
            render.SetActiveCamera(camera)
            render.SetViewport([right, 0.0, left, 1.0])
            render.SetBackground(255.0,255.0,255.0)
            renderers.append(render)
            right += width
            left += width

        return renderers

    '''
    def onFlatten(self,pickerIds, inflateStruc):
        for ren in self.renderers:
            self.ren.GetRenderWindow().RemoveRenderer(ren)

        self.meshProcessor.meshInteractor.improveCells(pickerIds)
        self.projectPass()

    def onGenerateColorMesh(self, renId, ids):
        actor = self.meshProcessor.meshInteractor.generateColorMesh(renId,ids)
        if self.renderers[renId].GetActors().GetNumberOfItems() > 1:
            self.renderers[renId].RemoveActor(self.renderers[renId].GetActors().GetLastActor())
        self.renderers[renId].AddActor(actor)
    '''

    def clearOutputDirectory(self):
        '''
        Deletes all files in the "out/3D" folder, the "unfolded" subdirectory remains untouched.
        :return:
        '''
        directory = os.path.join(self.dirname, "../out/3D")
        for entry in os.scandir(directory):
            if entry.is_file():
                os.remove(entry.path)

        #dirName = os.path.join(directory,"/unfolded")
        #os.mkdir(dirName)

    def boolean(self):
        self.hierarchical_mesh_anchor.recursive_difference()

    def onUnfoldTest(self):
        self.meshProcessor.unfoldTest(name="lower")

    def writeObjMeshes(self):
        self.meshProcessor.writeObjMeshes(self.hierarchical_mesh_anchor)

    def renderCutPlanes(self, viewpoints, bds):
        for i in range(len(viewpoints)):
            plane = vtk.vtkPlaneSource()
            x0, y0, z0 = (bds[i][1]-abs(bds[i][0]))/2., (bds[i][3]-abs(bds[i][2]))/2., (bds[i][5]-abs(bds[i][4]))/2.
            plane.SetCenter(x0,y0,z0)
            plane.SetNormal(viewpoints[i][0],viewpoints[i][1],viewpoints[i][2])
            plane.Update()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(plane.GetOutput())
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.SetScale(50.0,50.0,50.0)
            actor.GetProperty().SetColor(1.0,0.0,0.0)

            self.ren.AddActor(actor)

    def cutHM(self, viewpoints, bds):
        self.hierarchical_mesh_anchor.addCutPlanes(viewpoints, bds)
        #self.hierarchical_mesh_anchor.cut()
        self.hierarchical_mesh_anchor.recursive_difference()

    def onCutTest(self):
        self.hierarchical_mesh_anchor.recursive_difference()

