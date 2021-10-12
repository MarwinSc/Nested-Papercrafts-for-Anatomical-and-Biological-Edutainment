import vtkmodules.all as vtk
from vtkmodules.numpy_interface.dataset_adapter import numpy_support
import numpy as np
import util

#Class responsible for 2D image related processing steps.
class ImageProcessor():

    width = 2000
    height = 2000

    renList = []

    canvas_source = vtk.vtkImageCanvasSource2D()
    canvas_source.SetExtent(0, width-1, 0, height-1, 0, 0)
    canvas_source.SetScalarTypeToUnsignedChar()
    canvas_source.SetNumberOfScalarComponents(3)
    canvas_source.SetDrawColor(255,255,255)
    canvas_source.FillBox(0, width-1, 0, height-1)
    canvas_source.Update()

    #Main method to muliply structures
    def multiplyingActors(self,dethPeeling,filter,brightBool,actorList,camera,height,width,occlusion,numberOfPeels):

        self.renList.clear()

        #render every actor in an own renderer window
        for a in actorList:

            ren, iren, renWin, wti = util.getbufferRenIntWin(camera,height,width)
            renWin.SetOffScreenRendering(True)
            self.renList.append([ren, iren, renWin, wti])

            if dethPeeling:
                ren.SetUseDepthPeeling(True)
                ren.SetOcclusionRatio(occlusion)
                ren.SetMaximumNumberOfPeels(numberOfPeels)
            else:
                ren.SetUseDepthPeeling(False)

            ren.AddActor(a)
            renWin.Render()
            wti.Update()

        #if only one actor
        if len(actorList) == 1:
            result = self.renList[0][3].GetOutput()
            if brightBool:
                result = self.optimizedBrighten(result,self.width,self.height,"0")

        elif len(actorList) >= 2:
            #first two actors
            image = self.renList[0][3].GetOutput()
            image2 = self.renList[1][3].GetOutput()
            if brightBool:
                image = self.optimizedBrighten(image,self.width,self.height,"0")
                image2 = self.optimizedBrighten(image2,self.width,self.height,"1")

            result = self.normalizeMultiplication(image,image2,self.width,self.height).GetOutput()

            #the rest
            i = 0
            for a in actorList[2:]:
                img = result
                img2 = self.renList[2+i][3].GetOutput()
                if brightBool:
                    img2 = self.optimizedBrighten(img2,self.width,self.height,"2")
                result = self.normalizeMultiplication(img, img2,self.width,self.height).GetOutput()
                i = i + 1

        if filter:
            img = result
            img2 = self.canvas_source.GetOutput()
            result = self.normalizeMultiplication(img, img2,self.width,self.height).GetOutput()

        return result

    #method carrying out the normalization and multiplication of the structurs
    def normalizeMultiplication(self, image, image2):

        img1 = image / 255
        img2 = image2 / 255

        result = (img1 * img2)

        result = result * 255

        return result

    def optimizedBrighten(self,image):

        img1 = image / 255

        averagePixelRed = np.mean(img1[:,:,0])
        averagePixelGreen = np.mean(img1[:,:,1])
        averagePixelBlue = np.mean(img1[:,:,2])

        x = img1[:,:,0] + img1[:,:,1] + img1[:,:,2]

        white = np.where(x>2.0)
        black = np.where(x<=1.1)
        #colored = np.where(np.logical_and(x>1.0,x<2.0))
        colored = np.where(x<=2.0)

        if averagePixelBlue > averagePixelGreen and averagePixelRed > averagePixelGreen:
            img1[colored[0],colored[1],1] = img1[colored[0],colored[1],1] + (((img1[colored[0],colored[1],0]+img1[colored[0],colored[1],2])*0.5))
            img1[colored[0],colored[1],0] = 1.0
            img1[colored[0],colored[1],2] = 1.0

            img1[black[0],black[1],0] = 1.0
            img1[black[0],black[1],1] = 0.0
            img1[black[0],black[1],2] = 1.0
        if averagePixelGreen > averagePixelBlue and averagePixelRed > averagePixelBlue:
            img1[colored[0],colored[1],2] = img1[colored[0],colored[1],2] + (((img1[colored[0],colored[1],0]+img1[colored[0],colored[1],1])*0.5))
            img1[colored[0],colored[1],1] = 1.0
            img1[colored[0],colored[1],0] = 1.0

            img1[black[0],black[1],0] = 1.0
            img1[black[0],black[1],1] = 1.0
            img1[black[0],black[1],2] = 0.0
        if averagePixelGreen > averagePixelRed and averagePixelBlue > averagePixelRed:
            img1[colored[0],colored[1],0] = img1[colored[0],colored[1],0] + (((img1[colored[0],colored[1],2]+img1[colored[0],colored[1],1])*0.5))
            img1[colored[0],colored[1],1] = 1.0
            img1[colored[0],colored[1],2] = 1.0

            img1[black[0],black[1],0] = 0.0
            img1[black[0],black[1],1] = 1.0
            img1[black[0],black[1],2] = 1.0
        img1[white[0], white[1], :] = 1.0

        img1 = np.clip(img1 * 255,0.0,255.0)

        return img1


