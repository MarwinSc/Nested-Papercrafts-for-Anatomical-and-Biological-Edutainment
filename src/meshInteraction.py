import vtkmodules.all as vtk
import numpy as np

class MeshInteraction():

    def __init__(self,dedicatedPaperMeshes):
        self.dedicatedPaperMeshes = dedicatedPaperMeshes

    def extractMeanNormal(self, mesh, ids):

        normalsFilter = vtk.vtkPolyDataNormals()
        normalsFilter.ComputeCellNormalsOn()

        normalsFilter.SetInputData(mesh)
        normalsFilter.Update()
        poly = normalsFilter.GetOutput()

        normals = poly.GetCellData().GetArray("Normals")

        mean = np.array([0.0 ,0.0 ,0.0])
        for i in ids:
            normal = normals.GetTuple3(i)
            mean += normal
        mean = mean / len(ids)
        return mean

    def extractMeanOrigin(self, mesh, ids):
        mean = np.array([0.0 ,0.0 ,0.0])
        for i in ids:
            p = [0.0 ,0.0 ,0.0]
            mesh.GetPoint(i ,p)
            mean += np.array(p)
        mean = mean / len(ids)
        np.round(mean ,2)
        return mean

    def cellIdsToPointIDs(self, mesh, ids):
        result = []
        for i in ids:
            pntIDs = mesh.GetCell(i).GetPointIds()
            result.append(pntIDs.GetId(0))
            result.append(pntIDs.GetId(1))
            result.append(pntIDs.GetId(2))
        return np.unique(np.array(result))

    def improveCells(self, pickerIds):

        for j in range(len(pickerIds)):

            self.dedicatedPaperMeshes[j].GetProperty().SetOpacity(0.8)
            clean = vtk.vtkCleanPolyData()
            clean.SetInputData(self.dedicatedPaperMeshes[j].GetMapper().GetInput())
            clean.Update()
            mesh = clean.GetOutput()

            for i in range(len(pickerIds[j])):

                newGeometry = vtk.vtkPolyData()
                newGeometry.DeepCopy(mesh)

                pointIds = self.cellIdsToPointIDs(mesh ,pickerIds[j][i])

                normal = self.extractMeanNormal(newGeometry ,pickerIds[j][i])
                origin = self.extractMeanOrigin(newGeometry ,pointIds)

                plane = vtk.vtkPlane()
                plane.SetOrigin(origin)
                plane.SetNormal(normal)

                for i in pointIds:
                    point1 = newGeometry.GetPoint(i)

                    projected1 = [0.0 ,0.0 ,0.0]
                    plane.ProjectPoint(point1 ,projected1)
                    newGeometry.GetPoints().SetPoint(i ,projected1)

                mesh = newGeometry

            self.dedicatedPaperMeshes[j].GetMapper().SetInputData(mesh)

    def generateColorMesh(self, renId, ids):

        newGeometry = vtk.vtkPolyData()
        newGeometry.DeepCopy(self.dedicatedPaperMeshes[renId].GetMapper().GetInput())

        cellData = vtk.vtkUnsignedCharArray()
        cellData.SetNumberOfComponents(4)

        for i in range(newGeometry.GetNumberOfCells()):
            if i in ids[-1]:
                cellColor = [0.0, 150.0, 0.0, 155.0]
                cellData.InsertTuple(i, cellColor)
            elif (any(i in list for list in ids[:-1])):
                cellColor = [100.0, 0.0, 0.0, 155.0]
                cellData.InsertTuple(i, cellColor)
            else:
                cellColor = [0.0, 0.0, 0.0, 0.0]
                cellData.InsertTuple(i, cellColor)

        offsetted = vtk.vtkWarpVector()
        offsetted.SetInputData(newGeometry)
        offsetted.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,
                                         vtk.vtkDataSetAttributes.NORMALS)
        offsetted.SetScaleFactor(1.0)
        offsetted.Update()
        newGeometry = offsetted.GetOutput()

        newGeometry.GetCellData().SetScalars(cellData)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(newGeometry)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        return actor


