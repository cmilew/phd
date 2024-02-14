import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir("C:/Users/candi/Desktop/QDC8/QDC8_soft/QDC_V1/")
data_file = "ZDATA.bin"

dt = np.dtype('uint16')

f = open(data_file,"rb")
data = f.read()
#print(data)
print(dt)
zdata = np.frombuffer(data, dt)
print(zdata)

CPT = []

CPT1 = []
CPT2 = []
CPT3 = []
CPT4 = []
CPT5 = []
CPT6 = []
CPT7 = []
CPT8 = []

plt.figure()

nb_event = np.size(zdata)//27
print(nb_event)


for num_piste in range(8):
    for num_event in range(nb_event):
        CPT.append(np.uint32(((zdata[3 + num_piste * 2 + num_event * 27]) << 16) + (zdata[4 + num_piste * 2 + num_event * 27])))
    plt.plot([i for i in range(nb_event)],CPT)
    CPT.clear()
    

for num_event in range(nb_event):
    CPT1.append(np.uint32(((zdata[3 + 0 * 2 + num_event * 27]) << 16) + (zdata[4 + 0 * 2 + num_event * 27])))
    CPT2.append(np.uint32(((zdata[3 + 1 * 2 + num_event * 27]) << 16) + (zdata[4 + 1 * 2 + num_event * 27])))
    CPT3.append(np.uint32(((zdata[3 + 2 * 2 + num_event * 27]) << 16) + (zdata[4 + 2 * 2 + num_event * 27])))
    CPT4.append(np.uint32(((zdata[3 + 3 * 2 + num_event * 27]) << 16) + (zdata[4 + 3 * 2 + num_event * 27])))
    CPT5.append(np.uint32(((zdata[3 + 4 * 2 + num_event * 27]) << 16) + (zdata[4 + 4 * 2 + num_event * 27])))
    CPT6.append(np.uint32(((zdata[3 + 5 * 2 + num_event * 27]) << 16) + (zdata[4 + 5 * 2 + num_event * 27])))
    CPT7.append(np.uint32(((zdata[3 + 6 * 2 + num_event * 27]) << 16) + (zdata[4 + 6 * 2 + num_event * 27])))
    CPT8.append(np.uint32(((zdata[3 + 7 * 2 + num_event * 27]) << 16) + (zdata[4 + 7 * 2 + num_event * 27])))


plt.figure("CPTx")
plt.subplot(421)
plt.title("CPT1")
plt.plot([i for i in range(nb_event)],CPT1)

plt.subplot(422)
plt.title("CPT2")
plt.plot([i for i in range(nb_event)],CPT2)

plt.subplot(423)
plt.title("CPT3")
plt.plot([i for i in range(nb_event)],CPT3)

plt.subplot(424)
plt.title("CPT4")
plt.plot([i for i in range(nb_event)],CPT4)

plt.subplot(425)
plt.title("CPT5")
plt.plot([i for i in range(nb_event)],CPT5)

plt.subplot(426)
plt.title("CPT6")
plt.plot([i for i in range(nb_event)],CPT6)

plt.subplot(427)
plt.title("CPT7")
plt.plot([i for i in range(nb_event)],CPT7)

plt.subplot(428)
plt.title("CPT8")
plt.plot([i for i in range(nb_event)],CPT8)

plt.show()
