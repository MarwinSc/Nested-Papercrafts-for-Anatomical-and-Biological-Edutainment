import vtkmodules.all as vtk
import os
from meshProcessing import MeshProcessing
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
    imageProcessor = ImageProcessor()


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


    def addMesh(self, meshes):
        '''
        Adds a new hierarchicalMesh object to the hierarchy.
        '''
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


    def project(self,resolution=[500,500]):
        self.hierarchical_mesh_anchor.project(resolution)

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

    def writeUnfolded(self):
        self.hierarchical_mesh_anchor.writeAllUnfoldedMeshesInHierarchy()

    def importUnfolded(self):
        self.hierarchical_mesh_anchor.importAllUnfoldedMeshesInHierarchy()