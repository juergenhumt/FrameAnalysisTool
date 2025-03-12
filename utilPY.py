import re, sys
import numpy as np
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
from math import cos, sin, asin, acos, atan2, radians, degrees, sqrt, pi
from typing import Tuple

"""

 This file is part of FrameAnalysisTool
 

     FrameAnalysisTool, is free  software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by the 
     Free Software Foundation, either version 3 of the License or any later 
     version.
 
     FrameAnalysisTool is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
 
     You should have received a copy of the GNU General Public License along 
     with FrameAnalysisTool.  If not, see <http://www.gnu.org/licenses/>.

 Copyright 2025 Juergen Humt

"""
# App.getDocument('PietFusSpMkXII').getObject('Sketch007').movePoint(0,2,App.Vector(0.000646,282.052562,0),0)


def sinc(x):
    """Sinc function as used in MATLAB (normalized by pi)."""
    return np.sinc(x / np.pi)

def biarc_plot(x0, y0, l0, theta0, kappa0,
               x1, y1, l1, theta1, kappa1,
               fmt1=None, fmt2=None):
    """
    Plot a biarc.

    Arguments:
    x0, y0    -- initial point of the first arc
    l0        -- length of the first arc
    theta0    -- initial angle of the first arc
    kappa0    -- curvature of the first arc
    x1, y1    -- final point of the second arc
    l1        -- length of the second arc
    theta1    -- final angle of the second arc
    kappa1    -- curvature of the second arc
    fmt1      -- format dictionary for the first arc (e.g., {'color': 'blue', 'linewidth': 3})
    fmt2      -- format dictionary for the second arc (e.g., {'color': 'red', 'linewidth': 3})
    """

    # Default plot styles if not provided
    if fmt1 is None:
        fmt1 = {'color': 'blue', 'linewidth': 3}
    if fmt2 is None:
        fmt2 = {'color': 'red', 'linewidth': 3}

    # First arc
    ell = np.linspace(0, l0, 100)  # Equivalent to MATLAB 0:l0/100:l0
    tmp = (kappa0 / 2) * ell
    S = sinc(tmp)
    x = x0 + ell * S * np.cos(theta0 + tmp)
    y = y0 + ell * S * np.sin(theta0 + tmp)
    
    plt.plot(x, y, **fmt1)  # Use Python unpacking for dictionary arguments

    # Second arc
    ell = np.linspace(0, l1, 100)
    tmp = (kappa1 / 2) * ell
    S = sinc(tmp)
    x = x1 - ell * S * np.cos(theta1 - tmp)
    y = y1 - ell * S * np.sin(theta1 - tmp)

    plt.plot(x, y, **fmt2)

    plt.axis('equal')  # Keep aspect ratio square for correct arc visualization
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.title("Biarc Plot")
    plt.grid()

'''
# Example usage:
plt.figure()  # Create a new figure
biarc_plot(0, 0, 2, np.pi/4, 0.5, 2, 2, 2, np.pi/2, -0.5, {'color': 'blue', 'linewidth': 3}, {'color': 'red', 'linewidth': 3})
plt.show()  # Show the plot
atan2(-0.1, 0.1) + pi
'''

def circlefit(x, y):
# Circle fitting function
# By Izhak Bucher (2025).
# See legal notes in the description
    x = np.array(x)
    y = np.array(y)
    A = np.c_[2 * x, 2 * y, np.ones(len(x))]
    b = x**2 + y**2
    c = np.linalg.lstsq(A, b, rcond=None)[0]
    x_center, z_center = c[0], c[1]
    radius = np.sqrt(c[2] + x_center**2 + z_center**2)
    return x_center, z_center, radius


# Helper function to convert degrees to radians
def deg_to_rad(deg):
    return np.deg2rad(deg)

def atan36(z,x):
  a = atan2(z,x) % (2*pi)
  return a

   
def rotMat(alf, kDir=1):
  fhi = np.radians(alf)
  c, s = np.cos(fhi), np.sin(fhi)
  rM = np.array([[c, kDir*s],[-kDir*s, c]])
  return rM

def gFnc(alf, Rq, vC1, vC3, R3, d, bta):
    """
    Computes rotation angle necessary to make two consecutive
    inverse curve elements tangent to the end circles.

    Parameters:
    alf : float
        Angle alpha in radians.
    Rq : float
        Radius Rq.
    vC1 : numpy.ndarray
        2x1 vector.
    vC3 : numpy.ndarray
        2x1 vector.
    R3 : float
        Radius R3.
    d : float
        Distance parameter.
    bta : float
        Angle beta in radians.

    Returns:
    rsX : float
        Computed result.
    """
    vD = vC1 + d * np.array([-np.cos(alf), np.sin(alf)])
    vR2 = vD + Rq * np.array([np.cos(bta - alf), np.sin(bta - alf)])

    vR2b = np.linalg.norm(vC3 - vR2)

    rsX = vR2b - (Rq + R3)

    return abs(rsX)


def r2als2(R2, R1, R3, d):
  # equation #2
  try:
    xh1 = ((d**2 - (R1 - R2)**2 - (R2 + R3)**2)/(-2*abs(R2 - R1)*(R2 + R3))) # eq. #1
    xh1 = acos(xh1)
  except ValueError:
    sys.exit('Math Error in r2als2, try incresing R2')
  return xh1

def r2fnc2(R2, R1, R3, d):
  
  aLs = r2als2(R2, R1, R3, d)
  aLsD = degrees(aLs)
  l_ = R2*sin(0.5*aLs)
  
  # eq #4 in #5
  xRet = ((l_**2 + R3**2 - 2*l_*R3*cos(0.5*(pi+aLs))) - (R2**2 + (R2+R3)**2 - 2*R2*(R2+R3)*cos(aLs)))
  return abs(xRet)


def r2fnc(R2, R1, R3, d, bta):
  btaR = radians(bta)
  
  aLs = r2als2(R2, R1, R3, d)
  
  aL = acos(((R3-R2)**2 - (R2+R1)**2 - d**2)/(-2*(R1+R2)*d)) # eq. #1
  lv2 = ((d**2 + R3**2 - 2*R3*d*cos(btaR)))  # eq.#7
  l = sqrt(lv2)
  gA = acos((R3**2 - l**2 - d**2)/(-2*l*d))
  
  # eq #8a in #7
  xRet = (l**2 + R1**2 - 2*l*R1*(cos(gA)*cos(aL) + sin(gA)*sin(aL)) - (2*R2**2*(1 + cos(aLs))))
  return abs(xRet)


def gAs(R1, R2, R3, d, cX= 1.0):
  # the origins of R1 and R3 are lying at a distance 
  # d from each other and a desired value of R2 is
  # entered. From these values the angle alfa is 
  # calculated which the arc of R2 spans from the
  # center of R1 to the center of R3
  # if the circle R2 is to touch R1/R3 concave
  # cX has to be -1.0
  # 
  R23 = (R2-cX*R3)
  R21 = (R2-cX*R1)
  xNr = 4*((R23)**2 + cX*(R3-R1)*R23)
  xEq = (d**2 - (R3-R1)**2)/xNr
  siAlf2 = sqrt(xEq)
  fhiRad = 2*asin(siAlf2)

  fhi31 = acos((R21**2 - d**2 - R23**2)/(-2*d*R23))
  fhi21 = pi - fhiRad - fhi31
  fhiRadD = degrees(fhiRad)
  return  [fhiRad, fhiRadD, fhi21, fhi31]


def gAs__(R1, R2, R3, d, cX1= 1.0, cX3= 1.0):
  # the origins of R1 and R3 are lying at a distance 
  # d from each other and a desired value of R2 is
  # entered. From these values the angle alfa is 
  # calculated which the arc of R2 spans from the
  # center of R1 to the center of R3
  # if the circle R2 is to touch R1/R3 concave
  # cX has to be -1.0
  # 
  R23 = (R2-cX3*R3)
  R21 = (R2-cX1*R1)
  xNr = 4*((R23)**2 + cX3*(R3-R1)*R23)
  xEq = (d**2 - (R3-R1)**2)/xNr
  siAlf2 = sqrt(xEq)
  fhiRad = 2*asin(siAlf2)

  fhi31 = acos((R21**2 - d**2 - R23**2)/(-2*d*R23))
  fhi21 = pi - fhiRad - fhi31
  fhiRadD = degrees(fhiRad)
  return  [fhiRad, fhiRadD, fhi21, fhi31]




def gAs_(R1, R2, R3, d):
  # the origins of R1 and R3 are lying at a distance 
  # d from each other and a desired value of R2 is
  # entered. From these values the angle alfa is 
  # calculated which the arc of R2 spans from the
  # center of R1 to the center of R2
  # 
  R23 = (R2-R3)
  R21 = (R2-R1)
  xNr = 4*((R23)**2 + (R3-R1)*R23)
  xEq = (d**2 - (R3-R1)**2)/xNr
  siAlf2 = sqrt(xEq)
  fhiRad = 2*asin(siAlf2)

  fhi31 = acos((R21**2 - d**2 - R23**2)/(-2*d*R23))
  fhi21 = pi - fhiRad - fhi31
  return  [fhiRad, degrees(fhiRad), fhi21, fhi31]


def fndTriaVrtx(d, R1, R2):
    # Calculate x-coordinate of the intersection point
    x = (R1**2 - R2**2 + d**2) / (2 * d)
    
    # Calculate y-coordinate (considering both positive and negative roots)
    y_squared = R1**2 - x**2
    if y_squared < 0:
        raise ValueError("Invalid triangle configuration: no solution exists.")
    
    y = sqrt(y_squared)
    
    # Return both possible coordinates for point P
    return (x, y), (x, -y)

def circle_segment( center, radius, alfR, start_angle, end_angle, ds= 4, cX=1.0):
    aLen = abs(radius*(end_angle - start_angle))
    nPts = int(np.floor(aLen/float(ds)))
    
    angles = np.linspace(np.radians(start_angle), np.radians(end_angle), nPts)
    if (cX > 0) and (angles[0] > angles[-1]):
        dIr = -1
    else:
        dIr = 1
    #    angles = angles[::-1]  # Reverse the points
   
    angles = angles + np.radians(alfR)
    x = center[0] + radius * np.cos(angles)
    y = center[1] + radius * np.sin(angles)

    return x, y, nPts-1, dIr



def getSubShapeIndex(shp, sub):
    '''
    shp : a shape
    sub : a subshape (Vertex, Edge or Face)
    returns the index i of sub in shp 
    e.g if sub is a Vertex in shp
    sub = shp.Vertexes[i]
    '''
    dict = {'Vertex' : 'Vertexes', 'Edge' : 'Edges', 'Face' : 'Faces'}

    try:
        subs = getattr(shp, dict[sub.ShapeType])
    except KeyError:
        return -1
    
    for i, f in enumerate(subs):
        if f.isSame(sub): 
            return i
    return -1

def quaternion_to_axis_angle(quaternion: Tuple[float, float, float, float]) -> Tuple[Tuple[float, float, float], float]:
    """
    Convert quaternion to axis-angle form.

    Axis-angle is a two-element tuple where
    the first element is the axis vector (x, y, z),
    and the second element is the angle in degrees.

    See:
        https://github.com/FreeCAD/FreeCAD/blob/0.19.2/src/Base/Rotation.cpp#L119-L140
        https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/index.htm
    """
    qx, qy, qz, qw = quaternion

    s = sqrt(1 - qw**2)
    normalization_factor = 1 if s < 0.001 else s
    x = qx / normalization_factor
    y = qy / normalization_factor
    z = qz / normalization_factor
    axis = (x, y, z)

    angle = degrees(2 * acos(qw))

    return (axis, angle)


def euler_to_quaternion(yaw: float,
                        pitch: float,
                        roll: float) -> Tuple[float, float, float, float]:
    """
    Convert Euler angles (in degrees) to quaternion form:
        q0 = x, q1 = y, q2 = z and q3 = w
    where the quaternion is specified by q = w + xi + yj + zk.

    See:
        https://github.com/FreeCAD/FreeCAD/blob/0.19.2/src/Base/Rotation.cpp#L632-L658
        https://en.wikipedia.org/wiki/Quaternion
    """
    y = radians(yaw)
    p = radians(pitch)
    r = radians(roll)

    c1 = cos(y / 2.0)
    s1 = sin(y / 2.0)
    c2 = cos(p / 2.0)
    s2 = sin(p / 2.0)
    c3 = cos(r / 2.0)
    s3 = sin(r / 2.0)

    qx = (c1 * c2 * s3) - (s1 * s2 * c3)
    qy = (c1 * s2 * c3) + (s1 * c2 * s3)
    qz = (s1 * c2 * c3) - (c1 * s2 * s3)
    qw = (c1 * c2 * c3) + (s1 * s2 * s3)

    return (qx, qy, qz, qw)

def v3toqtr(p1, p2, p3):
   # Define the points on the target plane pln1
   v1 = np.array(p1)
   v2 = np.array(p2)
   v3 = np.array(p3)
   
   # Step 1: Normal of the initial x-y plane
   normal_xy = np.array([0, 0, 1])  # Normal of x-y plane (z-axis)
   
   # Step 2: Calculate the normal of the target plane pln1 using cross product
   # Vector from v1 to v2 and v1 to v3
   v12 = v2 - v1
   v13 = v3 - v1
   
   # Cross product of v12 and v13 gives the normal of the plane
   normal_pln1 = np.cross(v12, v13)
   normal_pln1 = normal_pln1 / np.linalg.norm(normal_pln1)  # Normalize
   
   # Step 3: Calculate the quaternion to rotate the normal of the x-y plane to the normal of pln1
   # Align the z-axis to the normal of pln1
   rotation_align = R.align_vectors([normal_pln1], [normal_xy])[0]
   
   # Step 4: Calculate the translation vector to move the x-y plane to pln1's position
   translation_vector = v1  # Offset the x-y plane to start from point v1
   
   # Step 5: Retrieve the quaternion components from the rotation
   qtrX = rotation_align.as_quat()  # Format is [x, y, z, w]
   
   # Output the results
   print(f"Translation vector to move x-y plane to pln1: {translation_vector}")
   print(f"Quaternion for rotation from x-y plane to pln1: [x = {qtrX[0]}, y = {qtrX[1]}, z = {qtrX[2]}, w = {qtrX[3]}]")
    
   return [v1, qtrX]



def pts2rot(p1_,p2_,alf):
   # Define the vectors
   # v1 = np.array([0, 10, -20])
   # v2 = np.array([100, 20, -50])
   # v3 = np.array([100, 50, -30])  # Point to rotate
   # angle_deg = 10  # Rotation angle
   v1 = np.array([p1_[0],p1_[1],p1_[2]])
   v2 = np.array([p2_[0],p2_[1],p2_[2]]) 
 
   # Step 1: Define the rotation axis as the vector from v1 to v2 and normalize it
   v12 = v2 - v1
   d12 = np.linalg.norm(v12)
   vN12 = v12 / d12 

   v3 = v1 + 0.5*d12*np.array([vN12[0], vN12[1]*(1.0 + 0.25), vN12[2]])
   # Step 2: Translate v3 so that the rotation axis passes through the origin
   v3_translated = v3 - v1
   
   # Step 3: Create the quaternion representing a rotation about the axis
   rotation = R.from_rotvec(np.radians(alf) * vN12)
   
   # Step 4: Apply the rotation to the translated v3
   v3_rotated_translated = rotation.apply(v3_translated)
   
   # Step 5: Translate v3 back to its original position
   v3_rotated = v3_rotated_translated + v1
   
   # Output the result
   print(f"New coordinates of v3 after rotation: {v3_rotated}")
   return v3_rotated



def pts2qtr(p1_, p2_, alf):
   # Define the initial points and rotation
   # p1 = np.array([0, 0, -100])   # Point p1
   # p2 = np.array([100, 50, -150])  # Point p2

   p1 = np.array([p1_[0],p1_[1],p1_[2]])
   p2 = np.array([p2_[0],p2_[1],p2_[2]]) 

   fhi = np.radians(alf)  # Convert to radians
   
   # Step 1: Calculate the normal of the initial plane (x-y plane)
   nIniPln = np.array([0, 0, 1])  # Normal of x-y plane is along z-axis
   
   # Step 2: Calculate the normal of the new plane defined by points p1 and p2
   v12 = p2 - p1  # Vector from p1 to p2
   nNewPln = np.cross(v12, np.array([1, 0, 0]))  # Cross with x-axis to find normal
   nNewPln = nNewPln / np.linalg.norm(nNewPln)  # Normalize the new plane normal
   
   # Step 3: Calculate the quaternion that aligns the initial plane's normal with the new plane's normal
   rotation_align = R.align_vectors([nNewPln], [nIniPln])[0]  # Align z-axis to new normal
   
   # Step 4: Apply additional rotation of -5 degrees around the vector from p1 to p2
   rotation_about_vector = R.from_rotvec(fhi * v12 / np.linalg.norm(v12))  # Rotation around vector
   
   # Combine the two rotations
   total_rotation = rotation_align * rotation_about_vector  # Multiply the quaternions
   
   # Step 5: Calculate the offset vector between the initial plane and new plane (difference between p1 and original origin)
   offset_vector = p1  # Since the plane was originally in the x-y plane, the offset is simply p1
   
   # Step 6: Extract quaternion from the total rotation
   quaternion = total_rotation.as_quat()  # Quaternion in the form [x, y, z, w]
   
   # Output the results
   print(f"Offset vector: {offset_vector}")
   print(f"Quaternion representing the total transformation: [x = {quaternion[0]}, y = {quaternion[1]}, z = {quaternion[2]}, w = {quaternion[3]}]")
   # ============

def pts2qtr__(p1_, p2_, alf):
   # Define points p1 and p2
   # p1 = np.array([0, 0, -100])
   # p2 = np.array([100, 50, -150])
   p1 = np.array([p1_[0],p1_[1],p1_[2]])
   p2 = np.array([p2_[0],p2_[1],p2_[2]]) 

   # Step 1: Calculate the vector from p1 to p2 (direction vector)
   v_p1p2 = p2 - p1
   v_p1p2_norm = v_p1p2 / np.linalg.norm(v_p1p2)  # Normalize the direction vector
   
   # Step 2: Calculate the normal of the new plane
   # The new plane normal is orthogonal to the direction vector v_p1p2
   # Since the plane is defined to pass through these two points, the normal to the plane
   # should be perpendicular to v_p1p2 and form a plane with the xy-plane initially
   
   # Step 3: The initial normal of the xy-plane is along the z-axis [0, 0, 1]
   initial_normal = np.array([0, 0, 1])
   
   # Find the quaternion that rotates the initial normal [0, 0, 1] to align with the vector v_p1p2
   rotation_to_align = R.align_vectors([v_p1p2_norm], [initial_normal])[0]  # Align vectors returns (rotation, RMSD)
   
   # Step 4: Additional rotation of -5 degrees about the axis v_p1p2 (rotation vector between p1 and p2)
   additional_rotation = R.from_rotvec(np.deg2rad(alf) * v_p1p2_norm)
   
   # Step 5: Combine the two rotations (first alignment, then the additional -5 degrees)
   combined_rotation = additional_rotation * rotation_to_align  # Combine rotations
   
   # Step 6: Compute the quaternion for the total transformation
   qtr = combined_rotation.as_quat()  # Returns (x, y, z, w)
   
   # Step 7: Calculate the offset vector
   # The offset vector is the projection of p1 onto the plane defined by the normal
   # The equation of the plane is n · (r - p1) = 0, where n is the new normal
   # The offset vector moves the plane so that it passes through p1
   
   # Get the new normal after transformation
   new_normal = combined_rotation.apply(initial_normal)
   
   # The offset is the distance along the new normal to bring the plane to pass through p1
   offset_magnitude = np.dot(p1, new_normal)  # Project p1 onto the normal
   oVec = offset_magnitude * new_normal
   
   # Output results
   print(f"Quaternion for the complete transformation: {qtr}")
   print(f"Offset vector between the initial and new plane: {oVec}")
   
   qtrV = [qtr[0],qtr[1],qtr[2]]
   qtrA = qtr[3]

   return [oVec, qtrV, qtrA]

def pts2rot__(p1_, p2_, alf):
   # pts2rot calculates a third point which is
   # rotated about the vector from p1 to p2 by
   # the angle alfa with respect to the global 
   # z direction
   p1 = np.array([p1_[0],p1_[1],p1_[2]])
   p2 = np.array([p2_[0],p2_[1],p2_[2]])
   p12 = p2 - p1
   d12 = np.linalg.norm(p12)
   pZ  = np.array([0, 0, 1]) 
   vPZ = np.cross(p12,pZ)
   vPY = np.cross(p12,vPZ)
   
   # Define the points and rotation parameters
   p1 = np.array([0, 0, -100])   # Point p1
   p2 = np.array([100, 50, -150])  # Point p2
   p3 = np.array([50, 20, -120])   # Example point p3 lying 50 mm from v12 (you can adjust as needed)
   angle_deg = 5  # Rotation angle in degrees
   
   # Step 1: Define the rotation axis vector v12 and normalize it
   v12 = p2 - p1
   v12_unit = v12 / np.linalg.norm(v12)
   
   # Step 2: Translate p3 so the rotation axis passes through the origin
   p3_translated = p3 - p1
   
   # Step 3: Create the quaternion representing a 5-degree rotation around v12
   rotation = R.from_rotvec(np.radians(angle_deg) * v12_unit)
   
   # Step 4: Apply the rotation to p3 in its translated position
   p3_rotated_translated = rotation.apply(p3_translated)
   
   # Step 5: Translate p3 back to its original position
   p3_rotated = p3_rotated_translated + p1
   
   # Output the result
   print(f"New position of p3 after rotation: {p3_rotated}")
   return p3_rotated



def pts2qtr_(p1_, p2_, alf):
   # Define the points p1 and p2
   # p1 = np.array([0, 0, -100])
   # p2 = np.array([100, 50, -150])
   p1 = np.array([p1_[0],p1_[1],p1_[2]])
   p2 = np.array([p2_[0],p2_[1],p2_[2]]) 
   # Step 1: Calculate the offset vector (vector from p1 to p2)
   offset_vector = p2 - p1
   
   # Step 2: Normalize the offset vector (axis of rotation)
   rotation_axis = offset_vector / np.linalg.norm(offset_vector)
   
   # Step 3: Define the rotation angle (in degrees) and convert to radians
   rotation_angle_deg = alf
   rotation_angle_rad = np.radians(rotation_angle_deg)
   
   # Step 4: Calculate the quaternion representing the -5 degree rotation about the axis
   # scipy.spatial.transform.Rotation.from_rotvec uses the angle-axis representation to create a rotation
   rotation_quaternion = R.from_rotvec(rotation_angle_rad * rotation_axis)
   
   # Get the quaternion as (x, y, z, w)
   qtr = rotation_quaternion.as_quat()  # This returns (x, y, z, w)
   
   # Output the results
   print(f"Offset vector (p2 - p1): {offset_vector}")
   print(f"Normalized rotation axis: {rotation_axis}")
   print(f"Rotation quaternion (x, y, z, w): {qtr}")
    
   qtrV = [qtr[0],qtr[1],qtr[2]]
   qtrA = qtr[3]

   return [qtrV, qtrA]



def rot2qtr(alf1, alf2):
   # Step 1: Create the quaternion for the first rotation (-10 degrees around the x-axis)
   rotation1 = R.from_euler('x', alf1, degrees=True)
   
   # Step 2: Create the quaternion for the second rotation (-20 degrees around the new y-axis)
   rotation2 = R.from_euler('y', alf2, degrees=True)
   
   # Step 3: Combine the two rotations by quaternion multiplication
   combined_rotation = rotation1 * rotation2  # Quaternion multiplication in scipy
   
   # Step 4: Retrieve the quaternion components (w, x, y, z)
   quaternion = combined_rotation.as_quat()
   
   # Print the result
   print(f"Combined quaternion: w = {quaternion[3]}, x = {quaternion[0]}, y = {quaternion[1]}, z = {quaternion[2]}")
   
   # Step 5: Optional: Retrieve the combined axis and angle from the quaternion
   axis_angle = combined_rotation.as_rotvec()  # Rotation vector (axis * angle in radians)
   axis = axis_angle / np.linalg.norm(axis_angle)  # Normalize to get axis
   angle_rad = np.linalg.norm(axis_angle)  # Angle in radians
   angle_deg = np.rad2deg(angle_rad)  # Convert angle to degrees
   
   print(f"Combined rotation axis: {axis}")
   print(f"Combined rotation angle: {angle_deg} degrees")
   return [axis, angle_deg]
   


def splitStr(inpStr, spltChrs="[,;:]"):
  # split up string after ,;: and blank
  inpStr = inpStr.strip()

  # replace any of ,;: with a blank
  str1 = re.sub(spltChrs,' ',inpStr)
  # replace multiple blanks by just one
  str2 = " ".join(str1.split())
  return str2.split()


def splitPrStr(inpStr):
#split up string after [{(,;:)}] and blank
#
  inpStr=inpStr.strip()
#replace any of ,;: with a blank
  str1= re.sub('[{(,;:)}]',' ',inpStr)
#replace multiple blanks by just one
  str2= " ".join(str1.split())
  return str2.split()

def modWngDat(flName):
  inF = open(flName,'r')

  # skip first line  
  inpLine = inF.readline()

  wDat = []
  while True:
    inpLine = inF.readline()
    if not inpLine:
       break
    tk=splitStr(inpLine)
    wD=[]
    wD.append(float(tk[0]))
    wD.append(float(tk[1]))
    wDat.append(wD)

  return wDat

def rdWngDat(flName, jB):
# jB =0 -> inner      jB =1 ->outer
  inF = open(flName,'r')
  
  wDat = []
  while True:
    inpLine = inF.readline()
    if not inpLine:
       break
    tk=splitStr(inpLine)
    wDat.append(float(tk[0]))

  if (jB > 0.5):
    kOfs= 4
    wDat[2] = wDat[3]  # inner end is tip of inner segment
    wDat[3] = wDat[1]  # outer end is wing span

    for k in range(4,8):
      wDat[k]=wDat[k+kOfs]

  wDat[8] = wDat[3] - wDat[2]
  inF.close()
  return wDat


def rdObjDat2(flName):
  print('reading ' + flName)
# jB =0 -> inner      jB =1 ->outer
  inF = open(flName,'r')

  wDat = []
  while True:
    inpLine = inF.readline()
    if not inpLine:
       break
    tk=splitStr(inpLine)
    if len(tk) > 2:
        wD = []
        wD.append(int(tk[2]))
        wD.append(int(tk[0]))
        wD.append(tk[1])
        wDat.append(wD)


  inF.close()
  wDatS = sorted(wDat) 

  return wDatS
