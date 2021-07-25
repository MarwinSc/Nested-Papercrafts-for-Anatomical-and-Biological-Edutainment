import sys, os
import vtkmodules.all as vtk
import itertools, math, random
from PyQt5 import QtCore, QtWidgets
from palettable.colorbrewer.qualitative import Set3_12 as cmap
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtkmodules.util.numpy_support import vtk_to_numpy
from math import sin,cos,pi,sqrt,floor
from viewpointselection import entropy

# import viewpointselection.stability
# import entropy
# import stability

class ViewPointComputation(QtWidgets.QMainWindow):

    def __init__(self, actor_list = None, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)

        self.frame = QtWidgets.QFrame()
        self.setWindowTitle("Entropy")
        self.setGeometry(0, 0, 800, 500)

        self.vl = QtWidgets.QVBoxLayout()
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget)

        self.ren = vtk.vtkRenderer()
        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.renWin = self.vtkWidget.GetRenderWindow()

        if actor_list is None:
            hierarchical_actor_list = self.surface_rendering(r"../viewpointselection/spheres")
        else:
            hierarchical_actor_list = actor_list

        flat_actor_list = list(itertools.chain.from_iterable(hierarchical_actor_list))

        for actor in flat_actor_list:
            self.ren.AddActor(actor)

        self.ren.ResetCamera()

        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)

        self.show()
        self.iren.Initialize()

        self.max_viewpoints_list, self.bds_list = self.init_hierarchical_entropy(hierarchical_actor_list)
        self.init_view_with_cutplanes(flat_actor_list, self.max_viewpoints_list, self.bds_list)

        #self.stability(actor_list)

    def getMaxViewPoints_AndBounds(self):
        return self.max_viewpoints_list, self.bds_list

    def init_hierarchical_entropy(self, hierarchical_actor_list):
        max_viewpoints_list = []
        bds_list = []
        lvl=0
        for level in hierarchical_actor_list:
            self.ren.RemoveAllViewProps()
            for actor in level:
                self.ren.AddActor(actor)
            self.ren.Render()
            bds_list.append(self.ren.ComputeVisiblePropBounds())

            N_views = 100  # change this
            print("This is level %i" % lvl)
            print("Getting %i views, this will take a while" % N_views)
            array_of_images, array_of_camera_positions = self.get_views(N_views, self.renWin, self.ren)
            print("Calculating entropy")
            min_entropy, max_entropy = entropy.entropy_calculation(array_of_images)
            min_entropy_view, max_entropy_view = array_of_camera_positions[min_entropy], array_of_camera_positions[
                max_entropy]
            # print("Min entropy view is: ", min_entropy_view, "for level %i" % lvl)
            print("Max entropy view is: ", max_entropy_view, "for level %i" % lvl)
            self.show_max_entropy_view(max_entropy_view)
            lvl+=1
            max_viewpoints_list.append(max_entropy_view)

        return max_viewpoints_list,bds_list

    def init_view_with_cutplanes(self, flat_actor_list, max_viewpoints_list, bds_list):
        for actor in flat_actor_list:
            self.ren.AddActor(actor)
        self.ren.Render()

        max_viewpoints_list_new = max_viewpoints_list[:len(max_viewpoints_list) - 1]
        bds_list_new = bds_list[1:len(max_viewpoints_list)]

        for (max_entropy_view, bds, bds_out) in zip (max_viewpoints_list_new, bds_list_new, bds_list[:len(max_viewpoints_list) - 1]):
            x, y, z = max_entropy_view[0], max_entropy_view[1], max_entropy_view[2]
            x0, y0, z0 = (bds[1]-bds[0])/2., (bds[3]-bds[2])/2., (bds[5]-bds[4])/2.
            camera = self.ren.GetActiveCamera()
            camera.SetPosition(x, y, z)
            self.ren.Render()

            normal = self.ren.GetActiveCamera().GetViewPlaneNormal()

            plane = vtk.vtkPlaneSource()
            plane.SetCenter(0,0,0)
            plane.SetNormal(normal[0], normal[1], normal[2])
            plane.Update()

            dist1 = math.sqrt(vtk.vtkMath().Distance2BetweenPoints(plane.GetPoint1(), plane.GetOrigin()))
            dist2 = math.sqrt(vtk.vtkMath().Distance2BetweenPoints(plane.GetPoint2(), plane.GetOrigin()))
            dist = max(dist1, dist2)
            scale = max((abs(bds_out[1]-bds_out[0])/dist, abs(bds_out[3]-bds_out[2])/dist))

            transform = vtk.vtkTransform()
            transform.Scale((scale+10, scale+10, scale+10)) #10 is just a margin here

            transformFilter = vtk.vtkTransformPolyDataFilter()
            transformFilter.SetInputConnection(plane.GetOutputPort())
            transformFilter.SetTransform(transform)
            transformFilter.Update()

            plane_mapper = vtk.vtkPolyDataMapper()
            plane_mapper.SetInputData(transformFilter.GetOutput())
            plane_actor = vtk.vtkActor()
            plane_actor.GetProperty().SetColor(random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1))
            plane_actor.SetMapper(plane_mapper)
            self.ren.AddActor(plane_actor)
            plane_actor.SetPosition(x0,y0,z0)
        self.ren.Render()

    # when you need to grab the planes, for each plane do: plane_actor.GetPosition() (this is the origin) and plane.GetNormal() (this is the normal) and plane_actor.GetBounds() will give you the boundaries of the plane (xmin,xmax,...)

    #def stability(self, actor_list):
    #    return stability.Stability(actor_list).stability_calculation()

    def surface_rendering(self, mypath):
        hierarchical_list = []

        for level in [f.path for f in os.scandir(mypath) if f.is_dir()]:
            level_list = []
            for file in os.listdir(level):
                if file.endswith(".obj"):
                    level_list.append(os.path.join(level, file))
            hierarchical_list.append(level_list)

        flat_list = list(itertools.chain.from_iterable(hierarchical_list))
        r, g, b, alpha, lut = self.visual_properties(len(flat_list))

        hierarchical_actor_list = []
        j=0
        for hierarchical_level in hierarchical_list:
            actor_list = []
            for i in range(len(hierarchical_level)):
                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInputData(self.loadIVRData(hierarchical_level[i]))

                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                actor.GetProperty().SetColor(float(r[i+j]), float(g[i+j]), float(b[i+j]))
                actor.GetProperty().SetOpacity(alpha[i+j])
                actor_list.append(actor)
                j+=1
            hierarchical_actor_list.append(actor_list)
        return hierarchical_actor_list

    def loadIVRData(self, name):
        reader = vtk.vtkOBJReader()
        reader.SetFileName(name)
        reader.Update()

        return reader.GetOutput()

    def visual_properties(self, N):
        r, g, b, alpha, lut = [], [], [], [], []

        for i in range(N):
            if i >= 12:
                i = i - (int(i / 12) * 12)
            else:
                i = i
            rgb = cmap.mpl_colors[i]
            r_i, g_i, b_i = rgb[0], rgb[1], rgb[2]
            alpha_i = 0.5

            h, s, v = self.rgb2hsv(rgb[0], rgb[1], rgb[2])

            lut_i = vtk.vtkLookupTable()
            lut_i.SetRange(0, 1)
            lut_i.SetHueRange(h, h)
            lut_i.SetValueRange(v, v)
            lut_i.SetSaturationRange(s, s)
            lut_i.SetAlphaRange(0, float(alpha_i))
            lut_i.SetRampToLinear()
            lut_i.Modified()
            lut_i.Build()

            r.append(r_i)
            g.append(g_i)
            b.append(b_i)
            alpha.append(alpha_i)
            lut.append(lut_i)

        return r, g, b, alpha, lut

    def rgb2hsv(self, r, g, b):
        r, g, b = r / 255.0, g / 255.0, b / 255.0

        mx = max(r, g, b)
        mn = min(r, g, b)
        df = mx - mn

        if mx == mn:
            h = 0
        elif mx == r:
            h = (60 * ((g - b) / df) + 360) % 360
        elif mx == g:
            h = (60 * ((b - r) / df) + 120) % 360
        elif mx == b:
            h = (60 * ((r - g) / df) + 240) % 360
        if mx == 0:
            s = 0
        else:
            s = df / mx
        v = mx

        return h / 360.0, s, v

    def get_views(self, N, renWin, ren):

        array_of_images = []
        array_of_camera_positions = []

        camera = ren.GetActiveCamera()
        position = camera.GetPosition()
        focal = camera.GetFocalPoint()
        r = sqrt((position[0] - focal[0])*(position[0] - focal[0]) + (position[1] - focal[1])*(position[1] - focal[1]) + (position[2] - focal[2])*(position[2] - focal[2]))

        hN = floor(sqrt(N))
        for i in range(0, hN):
            for j in range(0, hN):
                theta = i * pi / hN
                phi = j * 2 * pi / hN
                x, y, z = focal[0] + r * sin(theta) * cos(phi), focal[1] + r * sin(theta) * sin(phi), focal[2] + r * cos(theta)
                camera.SetPosition(x, y, z)
                ren.Render()

                w2if = vtk.vtkWindowToImageFilter()
                w2if.SetInput(renWin)
                w2if.Update()
                vtk_image = w2if.GetOutput()

                width, height, _ = vtk_image.GetDimensions()
                vtk_array = vtk_image.GetPointData().GetScalars()
                components = vtk_array.GetNumberOfComponents()

                array_of_images.append(vtk_to_numpy(vtk_array).reshape(height, width, components))
                array_of_camera_positions.append([x, y, z])

        return array_of_images, array_of_camera_positions

    def show_max_entropy_view(self, max_entropy_view):
        x, y, z = max_entropy_view[0], max_entropy_view[1], max_entropy_view[2]
        camera = self.ren.GetActiveCamera()
        camera.SetPosition(x, y, z)
        self.ren.Render()

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    window = ViewPointComputation()

    sys.exit(app.exec_())
