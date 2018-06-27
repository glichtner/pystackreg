import os
import pystackreg
import numpy as np
from matplotlib import pyplot as plt
from skimage import io
data_path = os.path.join(os.getenv('DAVID_DATA_PATH'), "180614_Rut1_Papain", "Use These Files")

ref = io.imread(os.path.join(data_path, "pygreg_test_src.tif"))
mov = io.imread(os.path.join(data_path, "pygreg_test_rot.tif"))

print(ref.shape)
print(mov.shape)

out = pystackreg.register(ref,mov)

f, ax = plt.subplots(3,1)
ax[0].imshow(ref)
ax[1].imshow(mov)
ax[2].imshow(out)

plt.show()


