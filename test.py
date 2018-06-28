import os
from pystackreg import StackReg
import numpy as np
from matplotlib import pyplot as plt
from skimage import io
data_path = os.path.join(os.getenv('DAVID_DATA_PATH'), "180614_Rut1_Papain", "Use These Files")

ref = io.imread(os.path.join(data_path, "pygreg_test_src.tif"))
mov = io.imread(os.path.join(data_path, "pygreg_test_rot.tif"))

print(ref.shape)
print(mov.shape)



#m,refpts,movpts = pystackreg._register(ref,mov)
#print(m)
#print(refpts)
#print(movpts)
#
#out = pystackreg._transform(mov, m)

#out = sr.register_transform(ref, mov)

sr = StackReg(StackReg.RIGID_BODY)
sr.register(ref,mov)
out = sr.transform(mov)


f, ax = plt.subplots(3,1)
ax[0].imshow(ref)
ax[1].imshow(mov)
ax[2].imshow(out)

plt.show()


