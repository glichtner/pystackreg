import os
from pystackreg import StackReg
import numpy as np

from skimage import io

def test():
	data_path = os.path.join(os.getenv('DAVID_DATA_PATH'), "180614_Rut1_Papain", "Use These Files")

	ref = io.imread(os.path.join(data_path, "pygreg_test_src.tif"))
	mov = io.imread(os.path.join(data_path, "pygreg_test_rot.tif"))

	print(ref.shape)
	print(mov.shape)

	#out = pystackreg.register(ref,mov)

	sr = StackReg(StackReg.RIGID_BODY)
	sr.register(ref,mov)
	out = sr.transform(mov)




