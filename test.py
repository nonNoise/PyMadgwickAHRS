# -*- coding: utf-8 -*-
import time
import math
import MadgwickAHRS
import pandas as pd

q = MadgwickAHRS.MadgwickAHRS()

data = pd.read_csv('sample.csv', names=('gx', 'gy', 'gz','ax','ay','az','mx','my','mz'))

for i in range(len(data.index)) :
    start = time.time()
    acc = (data['ax'][i],data['ay'][i],data['az'][i])
    gyro = (data['gx'][i],data['gy'][i],data['gz'][i])
    meg = (data['mx'][i],data['my'][i],data['mz'][i])
    time.sleep(0.1)
    print(acc)
    print(gyro)
    print(meg)
    
    
    q.MadgwickAHRSupdate(
        gx = gyro[0]*math.pi/180.0,
        gy = gyro[1]*math.pi/180.0,
        gz = gyro[2]*math.pi/180.0,
        ax = acc[0],
        ay = acc[1],
        az = acc[2],
        mx = meg[0]+40,
        my = meg[1]+40,
        mz = meg[2]+40  
    )
    """
    q.MadgwickAHRSupdateIMU(
        gx = gyro[0]*math.pi/180.0,
        gy = gyro[1]*math.pi/180.0,
        gz = gyro[2]*math.pi/180.0,
        ax = acc[0],
        ay = acc[1],
        az = acc[2],
    )  
    """  
    print("-"*10)
    print(u"Acceleration [g]")
    print(acc)
    print(u"Gyro [deg/s]")
    print(gyro)
    print(u"Magnetism[mT]")
    print(meg)
    print(u"Roll")
    print(q.GetRoll())
    print(u"Pitch")
    print(q.GetPitch())
    print(u"Yaw")
    print(q.GetYaw())
    print(u"SampleFreq [Hz]")
    print(1/(time.time() - start))
    q.sampleFreq = 1/(time.time() - start)