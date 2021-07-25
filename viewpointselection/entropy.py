from collections import Counter
from math import log
import numpy as np
#from matplotlib import pyplot as plt
import os


def entropy_calculation(image_list):
    entropy_list = []
    N = 1
    for img in image_list:
        img_resized = img[::10,::10]
        size = img_resized.shape[0] * img_resized.shape[1]

        # Make a single 24-bit number for each pixel
        res = np.zeros((img_resized.shape[0], img_resized.shape[1], 0))
        img_flat = np.dot(img_resized.astype(np.uint32),[1,256,65536]).flatten()

        color_count = Counter(img_flat)
        print("For image %i " % N, "out of ", len(image_list), "images:")
        print("found %i different colors" % len(color_count))

        samples_probability = [float(count) / size for count in color_count.values()]

        norm = np.linalg.norm(np.array(samples_probability))
        if norm != 0:
            norm_samples_probability = np.array(samples_probability)/norm
        else:
            norm_samples_probability = np.array(samples_probability)

        e = -sum([p * log(p, 2) for p in norm_samples_probability])
        entropy_list.append(e)
        print (N, "out of ", len(image_list), "views")
        N+=1
        print("image entropy is: ", e)

    print("min image entropy is: ", np.array(entropy_list).min(), "for image: ", np.array(entropy_list).argmin())
    print("max image entropy is: ", np.array(entropy_list).max(), "for image: ", np.array(entropy_list).argmax())

    return np.array(entropy_list).argmin(), np.array(entropy_list).argmax()