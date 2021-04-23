import vtkmodules.all as vtk
import os
import meshProcessing
import imageProcessing
import util
from mu3d.mu3dpy.mu3d import Graph
from src.hierarchicalMesh import HierarchicalMesh

# Class responsible for Rendering Tasks, I/O and forwarding method calls to ImageProcessing and MeshProcessing.
class Organizer():

    def __init__(self,ren):
        self.ren = ren
        # todo move this somewhere meaningful
        self.hierarchical_mesh_anchor = HierarchicalMesh(None, None, "Empty")

    dirname = os.path.dirname(__file__)

    actorList = []
    #TODO change when improving nested system
    nestedList = []
    renList = []

    height = 2000
    width = 2000

    meshProcessor = meshProcessing.MeshProcessing()
    imageProcessor = imageProcessing.ImageProcessor(height,width)

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

    def setUp(self):
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

        imageActor = vtk.vtkImageActor()

        #firsttime multiplication
        if not self.window.HasRenderer(self.resultRen):
            self.window.AddRenderer(self.resultRen)

        self.resultRen.RemoveActor(imageActor)
        self.ren.SetViewport(self.noViewport)
        self.resultRen.SetViewport(self.fullViewport)

        result = self.imageProcessor.multiplyingActors(depthPeeling,filter,brighten,self.actorList,self.camera,self.height,self.width,self.occlusion,self.numberOfPeels)

        imageActor = vtk.vtkImageActor()
        imageActor.GetMapper().SetInputData(result)
        imageActor.Update()

        self.resultRen.AddActor(imageActor)
        self.resultRen.Render()

        filename = os.path.join(self.dirname, "../out/2D/multiply{}.png".format(self.sessionMultiplySaves))
        writer = vtk.vtkPNGWriter()
        writer.SetFileName(filename)
        writer.SetInputData(result)
        writer.Write()

        self.sessionMultiplySaves += 1

    #multiplys the created unfoldings of the structures into a single unfolded texture.
    def brightenMultiplication(self):
        reader = vtk.vtkPNGReader()
        imgList=[]

        for i in range(len(self.actorList)):

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

    def changeColor(self,value,idx):
        self.actorList[idx].GetProperty().SetColor(value)

    def changeOpacity(self, value, idx):
        self.actorList[idx].GetProperty().SetOpacity(value)

    def addActor(self, name, boolNested=False):
        stlReader = vtk.vtkSTLReader()
        stlReader.SetFileName(name)# + ".stl")

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(stlReader.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(self.ctf.GetColor(0.0))
        actor.GetProperty().SetOpacity(0.75)
        if not boolNested:
            self.actorList.append(actor)
        else:
            #TODO real obj based system...
            self.nestedList.append(actor)

        #couldn't import multiple structures with that method call enabled
        #self.hierarchical_mesh_anchor.add(actor, name)

        self.ren.AddActor(actor)

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
        if resultRenMaxView:
            self.ren.SetViewport(self.noViewport)
            self.resultRen.SetViewport(self.fullViewport)
        else:
            self.ren.SetViewport(self.fullViewport)
            self.resultRen.SetViewport(self.noViewport)

    def createPaperMeshPass(self, iterations ,inflateStruc):
        mainGraph = Graph()
        actor = self.meshProcessor.createPapermesh(self.actorList)
        #self.ren.AddActor(self.meshProcessor.tempPaperActor)
        success = self.meshProcessor.mu3dUnfoldPaperMesh(actor, mainGraph, iterations)
        if success:
            #should be removed again when a new paperactor is created todo
            self.ren.AddActor(self.meshProcessor.tempPaperActor)
            self.meshProcessor.createDedicatedMeshes(inflateStruc, self.actorList)

    def createDedicatedPaperMeshesPass(self,inflateStrucList):
        self.meshProcessor.createDedicatedMeshes(inflateStrucList,self.actorList)

    def importUnfoldedMeshPass(self, inflateStruc, name):
        #self.ren.RemoveActor(self.ren.GetActors().GetLastActor())

        actor = self.meshProcessor.importUnfoldedMesh(name)
        actor.GetProperty().SetOpacity(0.5)
        self.ren.AddActor(actor)
        self.meshProcessor.createDedicatedMeshes(inflateStruc,self.actorList)

    def projectPass(self, inflateStruc = None, resolution = [500,500]):

        self.ren.SetViewport([0.0, 0.0, 0.0, 0.0])
        actors = self.meshProcessor.project(inflateStruc,self.actorList, resolution)
        count = 0

        filename = os.path.join(self.dirname, "../out/2D/unfolding{}.png".format(count))
        readerFac = vtk.vtkImageReader2Factory()
        imageReader = readerFac.CreateImageReader2(filename)
        imageReader.SetFileName(filename)
        texture = vtk.vtkTexture()
        texture.SetInputConnection(imageReader.GetOutputPort())

        self.meshProcessor.tempPaperActor.SetTexture(texture)
        self.meshProcessor.tempPaperActor.Modified()

        util.writeObj(self.meshProcessor.tempPaperActor.GetMapper().GetInput(), "tempPaper")

        self.renderers = self.setUpResultRenderers(self.camera,len(actors)+1)
        self.ren.GetRenderWindow().AddRenderer(self.renderers[0])
        self.renderers[0].AddActor(self.meshProcessor.tempPaperActor)

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

        actors.insert(0,self.meshProcessor.tempPaperActor)
        self.dedicatedPaperMeshes = actors
        return self.renderers

    def finish(self):
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

        self.meshProcessor.tempPaperActor.SetTexture(texture)

        self.ren.RemoveAllViewProps()
        self.meshProcessor.tempPaperActor.GetProperty().SetOpacity(1.0)
        self.ren.AddActor(self.meshProcessor.tempPaperActor)

    def setUpResultRenderers(self,camera,count):

        renderers = []

        width = 1.0 / count
        right = 0.0
        left = width

        for i in range(count):

            render = vtk.vtkRenderer()
            render.SetActiveCamera(camera)
            render.SetViewport([right, 0.0, left, 1.0])
            render.SetBackground(255.0,255.0,255.0)
            renderers.append(render)
            right += width
            left += width

        return renderers

    def onFlatten(self,pickerIds, inflateStruc):
        for ren in self.renderers:
            self.ren.GetRenderWindow().RemoveRenderer(ren)

        self.meshProcessor.meshInteractor.improveCells(pickerIds)
        self.projectPass(inflateStruc)

    def onGenerateColorMesh(self, renId, ids):
        actor = self.meshProcessor.meshInteractor.generateColorMesh(renId,ids)
        if self.renderers[renId].GetActors().GetNumberOfItems() > 1:
            self.renderers[renId].RemoveActor(self.renderers[renId].GetActors().GetLastActor())
        self.renderers[renId].AddActor(actor)

    def boolean(self):
        self.meshProcessor.booleanCGAL(self.nestedList)

    def onUnfoldTest(self):
        self.meshProcessor.unfoldTest(name="lower")
        #self.binder.unfoldTest()
