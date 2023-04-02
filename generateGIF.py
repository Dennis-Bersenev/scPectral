from PIL import Image, ImageDraw
import glob
import cv2

images = []
n_cells = 456
n_genes = 100


"""
Using cv2 to make it into a video
"""
for i in range(n_cells):
    path = "./out/cell{i}.png".format(i = i + 1)
    images.append(cv2.imread(path))



height,width,layers=images[0].shape

video=cv2.VideoWriter('video.avi',cv2.VideoWriter_fourcc(*'XVID'), 6, (width,height))

for j in range(n_cells):
    video.write(images[j])

cv2.destroyAllWindows()
video.release()

"""
Using PIL to make it into a gif
"""
# for i in range(n_cells):
#     path = "./out/cell{i}.png".format(i = i + 1)
#     images.append(Image.open(path))


# images[0].save("heatmap_animated.gif", save_all=True, append_images=images[1:], duration=500)
