
# Implementation of Madgwick's IMU and AHRS algorithms.
# See: http:#www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
#
# Date			Author          Notes
# 29/09/2011	SOH Madgwick    Initial release
# 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
# Open-source resources available on this website are provided under the 
# GNU General Public Licence unless an alternative licence is provided in source.
#
#
# Thank you for such a wonderful source. 
# This source forked objects written in C and ported it to Python.
# By Y.Kitagami


import math

def invSqrt(x):
    import struct
    i = struct.unpack('>i', struct.pack('>f', x))[0]
    i = 0x5f3759df - (i >> 1)
    y = struct.unpack('>f', struct.pack('>i', i))[0]
    return y * (1.5 - 0.5 * x * y * y)

class MadgwickAHRS :

    def __init__(self,quaternion=None, beta=None):
        self.sampleFreq = 3.0 # sample frequency in Hz
        self.beta = 0.1     # 2 * proportional gain
        if quaternion is not None:
            self.quaternion = quaternion
        if beta is not None:
            self.beta = beta
        self.q0 = 1.0
        self.q1 = 0.0
        self.q2 = 0.0
        self.q3 = 0.0
    def MadgwickAHRSupdate(self,gx,gy,gz,ax,ay,az,mx,my,mz):
        # Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
        #if((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
        #    MadgwickAHRSupdateIMU(gx, gy, gz, ax, ay, az)
        #    return
        #}
        q0 = self.q0
        q1 = self.q1
        q2 = self.q2
        q3 = self.q3
        sampleFreq = self.sampleFreq
        beta = self.beta
        # Rate of change of quaternion from gyroscope
        qDot1 = 0.5 * (-q1 * gx - q2 * gy - q3 * gz)
        qDot2 = 0.5 * (q0 * gx + q2 * gz - q3 * gy)
        qDot3 = 0.5 * (q0 * gy - q1 * gz + q3 * gx)
        qDot4 = 0.5 * (q0 * gz + q1 * gy - q2 * gx)

        # Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
        if( (ax != 0.0) and (ay != 0.0) and (az != 0.0)) :

            # Normalise accelerometer measurement
            recipNorm = invSqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm   

            # Normalise magnetometer measurement
            recipNorm = invSqrt(mx * mx + my * my + mz * mz)
            mx *= recipNorm
            my *= recipNorm
            mz *= recipNorm

            # Auxiliary variables to avoid repeated arithmetic
            _2q0mx = 2.0 * q0 * mx
            _2q0my = 2.0 * q0 * my
            _2q0mz = 2.0 * q0 * mz
            _2q1mx = 2.0 * q1 * mx
            _2q0 = 2.0 * q0
            _2q1 = 2.0 * q1
            _2q2 = 2.0 * q2
            _2q3 = 2.0 * q3
            _2q0q2 = 2.0 * q0 * q2
            _2q2q3 = 2.0 * q2 * q3
            q0q0 = q0 * q0
            q0q1 = q0 * q1
            q0q2 = q0 * q2
            q0q3 = q0 * q3
            q1q1 = q1 * q1
            q1q2 = q1 * q2
            q1q3 = q1 * q3
            q2q2 = q2 * q2
            q2q3 = q2 * q3
            q3q3 = q3 * q3

            # Reference direction of Earth's magnetic field
            hx = mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 + _2q1 * my * q2 + _2q1 * mz * q3 - mx * q2q2 - mx * q3q3
            hy = _2q0mx * q3 + my * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my * q1q1 + my * q2q2 + _2q2 * mz * q3 - my * q3q3
            _2bx = math.sqrt(hx * hx + hy * hy)
            _2bz = -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 + _2q1mx * q3 - mz * q1q1 + _2q2 * my * q3 - mz * q2q2 + mz * q3q3
            _4bx = 2.0 * _2bx
            _4bz = 2.0 * _2bz

            # Gradient decent algorithm corrective step
            s0 = -_2q2 * (2.0 * q1q3 - _2q0q2 - ax) + _2q1 * (2.0 * q0q1 + _2q2q3 - ay) - _2bz * q2 * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            s1 = _2q3 * (2.0 * q1q3 - _2q0q2 - ax) + _2q0 * (2.0 * q0q1 + _2q2q3 - ay) - 4.0 * q1 * (1 - 2.0 * q1q1 - 2.0 * q2q2 - az) + _2bz * q3 * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            s2 = -_2q0 * (2.0 * q1q3 - _2q0q2 - ax) + _2q3 * (2.0 * q0q1 + _2q2q3 - ay) - 4.0 * q2 * (1 - 2.0 * q1q1 - 2.0 * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            s3 = _2q1 * (2.0 * q1q3 - _2q0q2 - ax) + _2q2 * (2.0 * q0q1 + _2q2q3 - ay) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5 - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5 - q1q1 - q2q2) - mz)
            recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3) # normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            # Apply feedback step
            qDot1 -= beta * s0
            qDot2 -= beta * s1
            qDot3 -= beta * s2
            qDot4 -= beta * s3


        # Integrate rate of change of quaternion to yield quaternion
        q0 += qDot1 * (1.0 / sampleFreq)
        q1 += qDot2 * (1.0 / sampleFreq)
        q2 += qDot3 * (1.0 / sampleFreq)
        q3 += qDot4 * (1.0 / sampleFreq)

        # Normalise quaternion
        recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
        q0 *= recipNorm
        q1 *= recipNorm
        q2 *= recipNorm
        q3 *= recipNorm

        self.q0 = q0
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3

    def MadgwickAHRSupdateIMU(self,gx,gy,gz,ax,ay,az):
        q0 = self.q0
        q1 = self.q1
        q2 = self.q2
        q3 = self.q3
        sampleFreq = self.sampleFreq
        beta = self.beta
        # Rate of change of quaternion from gyroscope
        qDot1 = 0.5 * (-q1 * gx - q2 * gy - q3 * gz)
        qDot2 = 0.5 * (q0 * gx + q2 * gz - q3 * gy)
        qDot3 = 0.5 * (q0 * gy - q1 * gz + q3 * gx)
        qDot4 = 0.5 * (q0 * gz + q1 * gy - q2 * gx) 
        if( (ax != 0.0) and (ay != 0.0) and (az != 0.0)) :
            # Normalise accelerometer measurement
            recipNorm = invSqrt(ax * ax + ay * ay + az * az)
            ax *= recipNorm
            ay *= recipNorm
            az *= recipNorm   

            # Auxiliary variables to avoid repeated arithmetic
            _2q0 = 2.0 * q0
            _2q1 = 2.0 * q1
            _2q2 = 2.0 * q2
            _2q3 = 2.0 * q3
            _4q0 = 4.0 * q0
            _4q1 = 4.0 * q1
            _4q2 = 4.0 * q2
            _8q1 = 8.0 * q1
            _8q2 = 8.0 * q2
            q0q0 = q0 * q0
            q1q1 = q1 * q1
            q2q2 = q2 * q2
            q3q3 = q3 * q3

            # Gradient decent algorithm corrective step
            s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay
            s1 = _4q1 * q3q3 - _2q3 * ax + 4.0 * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az
            s2 = 4.0 * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az
            s3 = 4.0 * q1q1 * q3 - _2q1 * ax + 4.0 * q2q2 * q3 - _2q2 * ay
            recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3) # normalise step magnitude
            s0 *= recipNorm
            s1 *= recipNorm
            s2 *= recipNorm
            s3 *= recipNorm

            # Apply feedback step
            qDot1 -= beta * s0
            qDot2 -= beta * s1
            qDot3 -= beta * s2
            qDot4 -= beta * s3
        
        # Integrate rate of change of quaternion to yield quaternion
        q0 += qDot1 * (1.0 / sampleFreq)
        q1 += qDot2 * (1.0 / sampleFreq)
        q2 += qDot3 * (1.0 / sampleFreq)
        q3 += qDot4 * (1.0 / sampleFreq)

        # Normalise quaternion
        recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3)
        q0 *= recipNorm
        q1 *= recipNorm
        q2 *= recipNorm
        q3 *= recipNorm

        self.q0 = q0
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3

    def GetRoll(self):
        q0 = self.q0
        q1 = self.q1
        q2 = self.q2
        q3 = self.q3
        roll = math.atan2( q0*q1 + q2*q3 ,  0.5 - q1*q1 - q2*q2) 
        roll = roll * 57.29578
        return roll
    def GetPitch(self):
        q0 = self.q0
        q1 = self.q1
        q2 = self.q2
        q3 = self.q3
        pitch = math.asin(-2.0 * (q1*q3 - q0*q2)) 
        pitch = pitch * 57.29578
        return pitch
    def GetYaw(self):
        q0 = self.q0
        q1 = self.q1
        q2 = self.q2
        q3 = self.q3
        yaw = math.atan2(q1*q2 + q0*q3, 0.5 - q2*q2 - q3*q3)
        yaw = yaw * 57.29578 + 180.0
        return yaw