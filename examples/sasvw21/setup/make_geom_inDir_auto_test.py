import sys
sys.path.insert(0, '../../../python/')

import numpy as np
import matplotlib.pyplot as plt
# import solps
import gitr
import math

solps_rz = 'assets/geom-SASV6/solps_rz.txt'
gitr_rz = 'assets/geom-SASV6/gitr_rz.txt'


def vectors(lines1, lines2, inDir, plt):

   x1 = lines1[:, 0]
   z1 = lines1[:, 1]
   x2 = lines1[:, 2]
   z2 = lines1[:, 3]
   slope = lines2[:, 4]
   
   for i in range(0,len(x1)):
       if slope[i]==0: 
           perpSlope = 1.0e12
       else:
           perpSlope = -np.sign(slope[i])/np.abs(slope[i]);

       rPerp = -inDir[i]/np.sqrt(perpSlope*perpSlope+1);
       zPerp = -inDir[i]*np.sign(perpSlope)*np.sqrt(1-rPerp*rPerp);
       plt.quiver([x1[i] + (x2[i]-x1[i])/2], [z1[i] + (z2[i]-z1[i])/2], [rPerp/10], [zPerp/10], width=0.0015, scale=5, headwidth=4)



# read data in ogr r,z wall geometry from original gemotry file
with open(gitr_rz) as f: solps_geom = f.readlines()[1:]
solps_geom[-1] += '\n' # need this for index counting in r,z extraction

# creates empty matricies for r_ogr and z_ogr
r_ogr = z_ogr = np.empty(0)
for row in solps_geom:
    rvalue = float(row.split(' ')[0])
    r_ogr = np.append(r_ogr, rvalue)

    zvalue = float(row.split(' ')[1][0:-1])
    z_ogr = np.append(z_ogr, zvalue)

# creates a closed geometery by putting the first data point at the end
r_ogr = np.append(r_ogr, r_ogr[0])
z_ogr = np.append(z_ogr, z_ogr[0])


plt.plot(r_ogr,z_ogr)

# plt.axis('scaled')

angle = 2*math.pi/(len(r_ogr)-1)
radius = .5

# Creates index same length as r and z
index = np.ones(len(r_ogr))

index = np.ones(100)
angle = 2*math.pi/(len(index)-1)

for i in range(len(index)):
    index[i] = i+1

# calibration r and z
cal_unit = .99
r_cal = r_ogr*cal_unit + min(r_ogr)*(1-cal_unit+.001)
z_cal = z_ogr*cal_unit

# making Calibration circle
x = np.ones(len(index))
y = np.ones(len(index))
for i in range(len(index)):
    x[i] = radius*math.cos(angle*index[i])+(min(r_ogr)+max(r_ogr))/2
    y[i] = radius*math.sin(angle*index[i])
   
plt.close()
# plt.plot(x,y)

# plt.axis('scaled')


plt.plot(r_cal,z_cal)

# plt.axis('scaled')

# make inDir and inDir limits
inDir1 = np.ones(len(index))
inDir2 = np.ones(len(r_ogr))

# # V1.0
# if (math.pi/2/angle/len(index)) > .25:
#     inDir1[round(math.pi/2/angle)-1] = -1
# else:
#     inDir1[round(math.pi/2/angle)-2] = -1

# if (3*math.pi/2/angle/len(index)) < .75:
#     inDir1[round(3*math.pi/2/angle)-1] = -1
# else:
#     inDir1[round(3*math.pi/2/angle)-2] = -1

# V2.0
# for i in range(len(index)):
#     if x[i] <= (max(x)+min(x))/2:
#         inDir1[int(index[i])-2] = -1
#         inDir1[int(index[i])-2] = -1

# V3.0 - success
# sets bounds for typical switch conditions
lower = 2
upper = 4
# creates list f1 of index values under 100
if len(index) >= 100:
    f1 = len(index) - round(len(index)/100)*100
else:
    print('index must have a length of 100 or larger')
# scrolls through f1 in mulitples of 4 and 8 from 3 and 2 respectively to set 
# bounds for special conditions
i=3
q = 2
while i < 100:
    if f1 == i:
        upper = 3
    if f1 == q:
        upper = 3
        lower = 3
    i = i+4
    q = q+8
# applies bounds
for i in range(len(index)):
    if index[i] >= index[(round(len(index)*.25)-lower)] and  index[i] <= index[(round(len(index)*.75)-upper)]:
        inDir1[int(index[i])] = -1

lines1 = gitr.gitr_lines_from_points(x,y)
lines2 = gitr.gitr_lines_from_points(r_ogr,z_ogr)
lines3 = gitr.gitr_lines_from_points(r_cal,z_cal)
# vectors(lines1, lines1, inDir1, plt)
# vectors(lines2, lines2, inDir2, plt)

def ex1(lines1, lines2, plt):

   x1 = lines1[:, 0]
   z1 = lines1[:, 1]
   x2 = lines1[:, 2]
   z2 = lines1[:, 3]
   y1 = lines2[:, 0]
   r1 = lines2[:, 1]
   y2 = lines2[:, 2]
   r2 = lines2[:, 3]
   slope = lines1[:, 4]
   inDir = np.ones(len(x1))

   for i in range(0,len(x1)):
        midx = x1[i] + (x2[i]-x1[i])/2
        midz = z1[i] + (z2[i]-z1[i])/2
        midy = y1[i] + (y2[i]-y1[i])/2
        midr = r1[i] + (r2[i]-r1[i])/2
        
        if slope[i]==0: 
            perpSlope = 1.0e12
        else:
            perpSlope = -np.sign(slope[i])/np.abs(slope[i]);
            rPerp = -inDir[i]/np.sqrt(perpSlope*perpSlope+1)
            zPerp = -inDir[i]*np.sign(perpSlope)*np.sqrt(1-rPerp*rPerp)
        
        
        plt.quiver([x1[i] + (x2[i]-x1[i])/2], [z1[i] + (z2[i]-z1[i])/2], [rPerp/10], [zPerp/10], width=0.0015, scale=5, headwidth=4)


ex1(lines2, lines3, plt)



# perpSlope = -np.sign(slope[i])/np.abs(slope[i])
# rPerp = -inDir[i]/np.sqrt(perpSlope*perpSlope+1)
# zPerp = -inDir[i]*np.sign(perpSlope)*np.sqrt(1-rPerp*rPerp)



# q = np.linspace(-1,1,len(index))
# w = np.zeros(len(index))
# plt.plot(q,w)
# inDir0 = np.ones(len(index))
# lines0 = gitr.gitr_lines_from_points(q,w)
# vectors(lines0, inDir0, plt)


















