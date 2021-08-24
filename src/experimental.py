import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy


class UnfoldingPrinter:

    def __init__(self, unfolded_model_path, unfolded_gt_path, unfolded_mirror_gt_path):
        self.unfolded_model = unfolded_model_path
        self.unfolded_gt = unfolded_gt_path
        self.unfolded_mirror_gt = unfolded_mirror_gt_path

        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(vtk.vtkNamedColors().GetColor3d('White'))
        self.window = vtk.vtkRenderWindow()
        self.window.AddRenderer(self.renderer)
        self.window.SetSize(720, 480)

        self.i_renderer = vtk.vtkRenderWindowInteractor()
        self.i_renderer.SetRenderWindow(self.window)
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.i_renderer.SetInteractorStyle(style)

        model_reader = self.__load_from_obj(self.unfolded_model)
        self.__render_triangle_obj(model_reader, self.renderer, color=[0.5, 0.5, 0.5])

        gt_reader = self.__load_from_obj(self.unfolded_gt)
        self.__render_triangle_obj(gt_reader, self.renderer, color=[0.3, 0.3, 0.3])

        mirror_gt_reader = self.__load_from_obj(self.unfolded_mirror_gt)
        self.__render_triangle_obj(mirror_gt_reader, self.renderer, color=[0.1, 0.1, 0.1])

        self.__label_gt(gt_reader, self.renderer, 0.1, [250, 125, 0])
        self.__label_gt(mirror_gt_reader, self.renderer, 0.1, [125, 125, 125])

    def __render_triangle_obj(self, reader, renderer, color=None):
        if color is None:
            color = [0.1, 0.1, 0.1]

        model_mapper = vtk.vtkPolyDataMapper()
        model_mapper.SetInputConnection(reader.GetOutputPort())
        model_mapper.SetColorModeToDefault()
        model_actor = vtk.vtkActor()
        model_actor.SetMapper(model_mapper)

        model_actor.GetProperty().SetColor(color)
        model_actor.GetProperty().SetLineWidth(2)
        model_actor.GetProperty().SetAmbient(0.4)
        model_actor.GetProperty().EdgeVisibilityOn()
        renderer.AddActor(model_actor)

    def __load_from_obj(self, path):
        reader = vtk.vtkOBJReader()
        reader.SetFileName(path)
        reader.Update()
        return reader

    def __label_gt(self, reader, renderer, scale=0.1, color=None):
        if color is None:
            color = [255, 0, 0]
        polydata = reader.GetOutput()
        points = vtk_to_numpy(polydata.GetPoints().GetData()).astype(float)
        indices = np.unique(points, axis=0, return_index=True)[1]
        points = [points[index] for index in sorted(indices)]

        id = 1
        for i in range(0, len(points), 4):
            label = vtk.vtkVectorText()
            label.SetText(str(id))
            id += 1
            pos = (points[i] + points[i + 1] + points[i + 2] + points[i + 3]) / 4
            lblMapper = vtk.vtkPolyDataMapper()
            lblMapper.SetInputConnection(label.GetOutputPort())

            # Set up an actor for the node label
            lblActor = vtk.vtkFollower()
            lblActor.SetMapper(lblMapper)
            lblActor.SetScale(scale, scale, scale)
            lblActor.SetPosition(pos[0] - scale / 2, pos[1] - scale / 2, 0.0)
            lblActor.GetProperty().SetColor(color[0], color[1], color[2])
            renderer.AddActor(lblActor)

    def display(self):
        self.i_renderer.Initialize()
        self.window.Render()
        self.i_renderer.Start()


def main():
    up = UnfoldingPrinter(r"C:\repositories\AnatomicalEdutainer\meshes\model.obj",
                          r"C:\repositories\AnatomicalEdutainer\meshes\gluetabs.obj",
                          r"C:\repositories\AnatomicalEdutainer\meshes\mirrorgt.obj")
    up.display()


if __name__ == "__main__":
    main()
