# This calculation program uses a Cartesian coordinate system with the y-axis along the long axis and the x-axis along the short axis of the first tape. Light intensity was numerically calculated by superimposing arbitrary tapes with angles up to three variables.

# Importing Modules
import numpy as np
import matplotlib.pyplot as plt
import time
import csv
from decimal import Decimal, ROUND_HALF_UP
import colorsys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
import math



#Tape number:x
xx = 3
#Tape number:p
pp = 3
#Tape number:z
zz = 3

#The number which determines the states of polarizer (closed/open/45°)
nn = 4
# Amount of change in optical path difference in calculation
step = 5
# K=100/∫_380^750 y(λ)S(λ)dλ  y(λ)：Color-matching function(CIE1931) / S(λ)：optical spectrum
K = 5.279886888
#sRGB（standard RGB）
tosRGB = [
        [3.2406, -1.5372, -0.4986],
        [-0.9689, 1.8758, 0.0415],
        [0.0557, -0.2040, 1.0570]
    ]


#_______________________________________________________________________________________________________________________________________________
## data file (optical retardation: Value at 380 nm normalized to 100)
file_name = "data/data_d_100.csv"
# data file (Light source: White LED light source)
file_name_lED = "data/data_light.csv"
# data file (Light source: White LED light source ②)
file_name_lED2 = "data/data_light2.csv"
# data file (Light Source: Smartphone-iphone_13_mini)
file_name_lED3 = "data/data_light_iphone13mini.csv"
# data file (Light source: White LED light source (surface-pro-7))
file_name_lED4 = "data/data_light_surface.csv"
# data file (Correction data: C2)
file_name_correct = "data/data_correct.csv"
# data file (color-matching function(CIE1931) row[１]=x(λ),row[2]=y(λ),rao[3]=z(λ))
file_name_color_f = "data/data_e_color.csv"
#________________________________________________________________________________________________________________________________________________

# Optical path difference (with open function, installation of optical path difference around 100nm)
with open(file_name, encoding='utf8') as f1:
    reader = csv.reader(f1)
    header = next(reader)

    λ = 380
    dd = []
    λλ = []
    for row in reader:
        λ = int(row[0])
        d = float(row[1])
        dd.append(d)
        λλ.append(λ)
        
        if int(row[0]) != λ:
            λ += 1
            # Reset arrays for repeated enforcement
            if λ > 750:
                break 

dd.append(d)
λλ.append(λ)

#Change the type of installed filedata from list to array to describe the polarization state
array_dd = np.array(dd)
array_λλ = np.array(λλ)

#＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿＿

# Light source (with open function to install white LED light source data)
with open(file_name_lED4, encoding='utf8') as f2:
    reader = csv.reader(f2)
    header = next(reader)

    λ = 380
    ll = []
    λλ = []
    for row in reader:
        λ = int(row[0])
        l = float(row[1])
        ll.append(l)
        
        if int(row[0]) != λ:
            λ += 1
            # Reset arrays for repeated enforcement
            if λ > 750:
                break 

ll.append(l)

#Change the type of installed filedata from list to array to describe the polarization state
array_ll = np.array(ll)

#_______________________________________________________________________________________________________________________________________________

# Correction data (with open function to install correction data)
with open(file_name_correct, encoding='utf8') as f3:
    reader = csv.reader(f3)
    header = next(reader)

    λ = 380
    cc = []
    λλ = []
    for row in reader:
        λ = int(row[0])
        c = float(row[1])
        cc.append(c)
        
        if int(row[0]) != λ:
            λ += 1
            # Reset arrays for repeated enforcement
            if λ > 750:
                break 
            
cc.append(c)

#Change the type of installed filedata from list to array to describe the polarization state
array_cc = np.array(cc)

#_______________________________________________________________________________________________________________________________________________

# Color-matching function (with open function, install color-matching function)
with open(file_name_color_f, encoding='utf8') as f4:
    reader = csv.reader(f4)
    header = next(reader)

    λ = 380
    cfx = []
    cfy = []
    cfz = []
    λλ = []
    for row in reader:
        λ = int(row[0])
        fx = float(row[1])
        fy = float(row[2])
        fz = float(row[3])
        cfx.append(fx)
        cfy.append(fy)
        cfz.append(fz)
        
        if int(row[0]) != λ:
            λ += 1
            # Reset arrays for repeated enforcement
            if λ > 750:
                break 
        
cfx.append(fx)
cfy.append(fy)
cfz.append(fz)

# Change the type of installed filedata from list to array to describe the polarization state
array_cfx = np.array(cfx)
array_cfy = np.array(cfy)
array_cfz = np.array(cfz)

#_______________________________________________________________________________________________________________________________________________

for count in range(230,235,5):
    #Optical path difference (normalized value multiplied by real number)
    array_ddc = [num*(count/100) for num in np.array(dd)]
    

    #List of initialized arrays for each condition : Each combination of retardation, tape angle, polarizer, and number of tapes

    # Tristimulus values : Array list ___________________________________________________________________________________________
    arr_spcX= np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_spcY = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_spcZ = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))

    # Array to store color mixing ratios for each condition
    arr_spcXx = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_spcYy = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_spcZz = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))

    # Array to store R"G"B"(standard RGB） for each condition
    arr_RR = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_GG = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_BB = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))

    # Array to store R'G'B'(standard RGB） for each condition 
    arr_R = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_G = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_B = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))

    # Array to store RGB for each condition
    arr_Rr = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_Gr = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_Br = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))

    # Array to store HSV for each condition
    arr_H = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_S = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_V = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))

    # Array of xyz, 3D coordinates in HSV color space
    arr_x = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_y = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))
    arr_z = np.zeros(((nn),(zz),(pp),(xx),int(180/step)+1, int(180/step)+1,int(180/step)+1))

    #_____________________________________________________________________________________________________________________________________

    def r_theta(theta):  # Rotation matrix R(θ)
        return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    def mai_r_theta(theta):  #Rotation matrix R(-θ)
        return np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])

    def jhons(theta):  # Jones Matrix
        return np.array([[np.sin(theta)**2, -np.sin(theta)*np.cos(theta)], [-np.sin(theta)*np.cos(theta), np.cos(theta)**2]])





    def main():  # Main function (mechanism for calculating an array of RGB values under all conditions)
        start_time = time.time()
        
        for nn in range(2,3): # open gate // ※ closed gate :(1,2), 45° gate :(3,4)
            for zz in range(1,3): # tape number z(from 1 to 2)
                for pp in range(1,3): # tape number p(from 1 to 2)
                    for xx in range(1,3): # tape number x(from 1 to 2)
                        for ε in range(0,181,step): #Third angle 
                            for β in range(0, 181, step):#Second angle 
                                for w in range(0, 181, step):#First angle or vertical polarizer angle 
                                    #Variable assigning the sum of XYZ (color mixing ratio)
                                    spcX = 0.0
                                    spcY = 0.0
                                    spcZ = 0.0
                                    # Rotation angle of the first polarizer (α in Fig. 2.)
                                    a = np.deg2rad(-w)
                                    E_1 = np.array([[-np.sin(a)], [np.cos(a)]])
                                    # Rotation angle of the second cellophane sheet (β in Fig. 2.)
                                    b = np.deg2rad(β-w)
                                    # Rotation angle of the third cellophane sheet (γ in Fig. 2.)
                                    e = np.deg2rad(ε-w)
                                    # Rotation angle of the second polarizer (ε in Fig. 2.)
                                    if nn == 1: # closed gate
                                        c = np.deg2rad(-w-90)
                                    elif nn ==2: # open gate
                                        c = np.deg2rad(-w)
                                    elif nn == 3: # 45° gate
                                        c = np.deg2rad(-w-45)
                                    
                                    #Jones Matrix calculation using optical path difference d at each λ
                                    for i in range(380,751): 
                                        # λ
                                        l = i
                                        ## Phase difference (for x number of tapes)
                                        delta = 2*array_ddc[l-380]*np.pi/l
                                        ##  Phase difference (for p number of tapes)
                                        deltap = 2*array_ddc[l-380]*np.pi/l
                                        ##  Phase difference (for z number of tapes)
                                        deltaz = 2*array_ddc[l-380]*np.pi/l
                                        # Cellophane's Jones Matrix formula (explicitly representing imaginary i in 1j): x pieces
                                        cellox = np.array([[1, 0], [0, np.exp(-1j*delta*xx)]])
                                        # Cellophane's Jones Matrix formula (explicitly representing imaginary i in 1j): p pieces
                                        cellop = np.array([[1, 0], [0, np.exp(-1j*deltap*pp)]])
                                        # Cellophane's Jones Matrix formula (explicitly representing imaginary i in 1j): z pieces
                                        celloz = np.array([[1, 0], [0, np.exp(-1j*deltaz*zz)]])
                                        # Conversion to E_2 (Polarizer and x pieces of tape)
                                        E_2 = np.dot(cellox, E_1)
                                        # Conversion to E_3 (Polarizer and, x number of tapes and p number of tapes)
                                        E_3 = np.dot(r_theta(b), np.dot(
                                            cellop, np.dot(mai_r_theta(b), E_2)))
                                        # Conversion to E_4 (Polarizer and, x number of tapes, p number of tapes, z number of tapes)
                                        E_4 = np.dot(r_theta(e), np.dot(
                                            celloz, np.dot(mai_r_theta(e), E_3)))
                                        # Conversion to E_5 (Finally, the effect of the second polarizer)
                                        E_5 = np.dot(jhons(c), E_4)
                                        # Calculate I (light intensity) at each wavelength
                                        I = (np.abs(np.abs(E_5[0]**2) + np.abs(E_5[1]**2)))
                                        # Calculate the measured spectrum by multiplying the intensity array by the light source data array　
                                        sp = I*(array_ll[l-380])
                                        # Corrected measurement configuration spectrums are calculated by integrating the correction data.
                                        spc = sp*array_cc[l-380]
                                        # Find the sum of the products of each of the three types of Color-matching function.
                                        spcX += float((spc*array_cfx[l-380])*K)
                                        spcY += float((spc*array_cfy[l-380])*K)
                                        spcZ += float((spc*array_cfz[l-380])*K)
                                    print("step:"+str(1),"nn:"+str(nn),"z:"+str(zz),"p:"+str(pp),"x:"+str(xx),"Third angle:"+str(ε),"Second angle:"+str(β), "First angle:"+str(w))
                                    # Assign the sum of the products of the three isochromatic functions to the array
                                    arr_spcX[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = spcX 
                                    arr_spcY[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = spcY
                                    arr_spcZ[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = spcZ
            
        # Calculate color mixing ratio xyz from arr_spc(X,Y,Z) for all conditions. Convert color mixing ratio xyz to R" G" B" values for all conditions
        for nn in range(2,3):
            for zz in range(1,3):
                for pp in range(1,3):
                    for xx in range(1,3):
                        for ε in range(0,181,step):
                            for β in range(0, 181, step):
                                for w in range(0, 181, step):
                                    arr_spcXx[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = (arr_spcX[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])/((arr_spcX[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])+(arr_spcY[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])+(arr_spcZ[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]))
                                    arr_spcYy[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = (arr_spcY[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])/((arr_spcX[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])+(arr_spcY[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])+(arr_spcZ[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]))
                                    arr_spcZz[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = (arr_spcZ[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])/((arr_spcX[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])+(arr_spcY[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])+(arr_spcZ[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]))
                                    #Create a matrix of color mixing ratios (3 rows by 1 column) as E5 and calculate the product with the sRGB conversion matrix as a [R"G"B"] array
                                    E_5 = np.array([[arr_spcXx[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]], [arr_spcYy[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]],
                                                    [arr_spcZz[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]]])
                                    arr_RGB = np.dot(tosRGB, E_5)
                                    arr_RR[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]= arr_RGB[0]
                                    arr_GG[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = arr_RGB[1]
                                    arr_BB[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = arr_RGB[2]
                                    print("step:"+str(2),"nn:"+str(nn),"z:"+str(zz),"p:"+str(pp),"x:"+str(xx),"Third angle:"+str(ε),"Second angle:"+str(β), "First angle:"+str(w))
                                    #R'G'B' values are calculated from the R"G"B" array.
                                    # This is called γ correction. Normally, sRGB γ correction is approximated by γ=2.2, but strictly speaking, it has a formula that results in a straight line at low luminance and a curve with γ=2.4 at high luminance. In this case, γ=2.4 is used for the correction.
                                    for color in [arr_RR, arr_GG, arr_BB]:
                                        value = color[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]
                                        if value <= 0.0031308:
                                            value *= 12.92 * 255
                                        else:
                                            value = (1.055 * (value ** (1 / 2.4)) - 0.055) * 255
                                        color[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = min(255, max(0, value))  
        
    
        # RGB values are calculated from R'G'B' values, rounded to string type by Decimal and converted to number type by .to_integral_value
        for nn in range(2,3):
            for zz in range(1,3):
                for pp in range(1,3):
                    for xx in range(1,3):
                        for ε in range(0,181,step):
                            for β in range(0, 181, step):
                                for w in range(0, 181, step):
                                    arr_Rr[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]= (Decimal(str(arr_RR[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])).quantize(Decimal('0'), rounding=ROUND_HALF_UP).to_integral_value())
                                    arr_Gr[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]= (Decimal(str(arr_GG[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])).quantize(Decimal('0'), rounding=ROUND_HALF_UP).to_integral_value())
                                    arr_Br[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]= (Decimal(str(arr_BB[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])).quantize(Decimal('0'), rounding=ROUND_HALF_UP).to_integral_value())
                                    print("step:"+str(3),"nn:"+str(nn),"z:"+str(zz),"p:"+str(pp),"x:"+str(xx),"Third angle:"+str(ε),"Second angle:"+str(β), "First angle:"+str(w))
                            
        # Converts RGB values to HSV values.
        for nn in range(2,3):
            for zz in range(1,3):
                for pp in range(1,3):
                    for xx in range(1,3):
                        for ε in range(0,181,step):
                            for β in range(0, 181, step):
                                for w in range(0, 181, step):
                                    HSV = colorsys.rgb_to_hsv((arr_Rr[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])/255,(arr_Gr[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])/255,
                                                            (arr_Br[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])/255)
                                    arr_H[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = HSV[0]
                                    arr_S[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = HSV[1]
                                    arr_V[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]= HSV[2]
                                    print("step:"+str(4),"nn:"+str(nn),"z:"+str(zz),"p:"+str(pp),"x:"+str(xx),"Third angle:"+str(ε),"Second angle:"+str(β), "First angle:"+str(w))
        
        # Add Figure
        fig = plt.figure(figsize=[24,12])

        # Add 3DAxes
        ax1 = fig.add_subplot(2, 4, 1)
        ax2 = fig.add_subplot(2, 4, 2)
        ax3 = fig.add_subplot(2, 4, 3)
        ax4 = fig.add_subplot(2, 4, 4)
        ax5 = fig.add_subplot(2, 4, 5)
        ax6 = fig.add_subplot(2, 4, 6)
        ax7 = fig.add_subplot(2, 4, 7)
        ax8 = fig.add_subplot(2, 4, 8)

        # Set the title of Axes
        ax1.grid(True)  # grid display ON
        ax1.set_xlim(-1,1)
        ax1.set_ylim(-1,1)
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")
        ax1.set_title("open gate z1-p1-x1/R(retardation)"+str(count))  # Title
        # Describes a circle on a graph (ec = edge color)
        c = patches.Circle(xy=(0, 0), radius=1.0,fill=False, ec='r')
        # Add a circle to the graph.
        ax1.add_patch(c)
        
        ax2.grid(True)  # grid display ON
        ax2.set_xlim(-1,1)
        ax2.set_ylim(-1,1)
        ax2.set_xlabel("x")
        ax2.set_ylabel("y")
        ax2.set_title("open gate z1-p1-x2/R(retardation)"+str(count))  # Title
        #Describes a circle on a graph (ec = edge color)
        c = patches.Circle(xy=(0, 0), radius=1.0,fill=False, ec='r')
        # Add a circle to the graph.
        ax2.add_patch(c)
        
        ax3.grid(True)  # grid display ON
        ax3.set_xlim(-1,1)
        ax3.set_ylim(-1,1)
        ax3.set_xlabel("x")
        ax3.set_ylabel("y")
        ax3.set_title("open gate z1-p2-x1/R(retardation)"+str(count))  # Title
        # Describes a circle on a graph (ec = edge color)
        c = patches.Circle(xy=(0, 0), radius=1.0,fill=False, ec='r')
        # Add a circle to the graph.
        ax3.add_patch(c)
        
        ax4.grid(True)  # grid display ON
        ax4.set_xlim(-1,1)
        ax4.set_ylim(-1,1)
        ax4.set_xlabel("x")
        ax4.set_ylabel("y")
        ax4.set_title("open gate z1-p2-x2/R(retardation)"+str(count))  # Title
        # Describes a circle on a graph (ec = edge color)
        c = patches.Circle(xy=(0, 0), radius=1.0,fill=False, ec='r')
        # Add a circle to the graph.
        ax4.add_patch(c)
        
        ax5.grid(True)  # grid display ON
        ax5.set_xlim(-1,1)
        ax5.set_ylim(-1,1)
        ax5.set_xlabel("x")
        ax5.set_ylabel("y")
        ax5.set_title("open gate z2-p1-x1/R(retardation)"+str(count))  # Title
        # Describes a circle on a graph (ec = edge color)
        c = patches.Circle(xy=(0, 0), radius=1.0,fill=False, ec='r')
        # Add a circle to the graph.
        ax5.add_patch(c)
        
        ax6.grid(True)  # grid display ON
        ax6.set_xlim(-1,1)
        ax6.set_ylim(-1,1)
        ax6.set_xlabel("x")
        ax6.set_ylabel("y")
        ax6.set_title("open gate z2-p1-x2/R(retardation)"+str(count))  # Title
        # Describes a circle on a graph (ec = edge color)
        c = patches.Circle(xy=(0, 0), radius=1.0,fill=False, ec='r')
        # Add a circle to the graph.
        ax6.add_patch(c)

        ax7.grid(True)  # grid display ON
        ax7.set_xlim(-1,1)
        ax7.set_ylim(-1,1)
        ax7.set_xlabel("x")
        ax7.set_ylabel("y")
        ax7.set_title("open gate z2-p2-x1/R(retardation)"+str(count))  # Title
        # Describes a circle on a graph (ec = edge color)
        c = patches.Circle(xy=(0, 0), radius=1.0,fill=False, ec='r')
        # Add a circle to the graph.
        ax7.add_patch(c)
        
        ax8.grid(True)  # grid display ON
        ax8.set_xlim(-1,1)
        ax8.set_ylim(-1,1)
        ax8.set_xlabel("x")
        ax8.set_ylabel("y")
        ax8.set_title("open gate z2-p2-x2/R(retardation)"+str(count))  # Title
        # Describes a circle on a graph (ec = edge color)
        c = patches.Circle(xy=(0, 0), radius=1.0,fill=False, ec='r')
        # Add a circle to the graph.
        ax8.add_patch(c)
        
        # Display HSV values in polar coordinates
        
        for nn in range(2,3):
            for zz in range(1,3):
                for pp in range(1,3):
                    for xx in range(1,3):
                        X = []
                        Y = []
                        for ε in range(0,181,step):
                            for β in range(0, 181, step):
                                for w in range(0, 181, step):
                                    arr_x[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = (arr_S[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])*np.cos(2*np.pi*(arr_H[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]))
                                    arr_y[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = (arr_S[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)])*np.sin(2*np.pi*(arr_H[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]))
                                    arr_z[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)] = arr_V[(nn)][(zz)][(pp)][(xx)][int(ε/step)][int(β/step)][int(w/step)]
                        X.append(arr_x)
                        Y.append(arr_y)
                        if xx == 1 and pp ==1 and zz==1 and nn==2:
                            color = "b"
                            ax1.scatter(X,Y,s=0.1,c=color)
                        elif xx == 2 and pp ==1 and zz==1 and nn == 2:
                            color = "r"
                            ax2.scatter(X,Y,s=0.1,c=color)
                        elif xx == 1 and pp ==2 and zz==1 and nn == 2:
                            color = "g"
                            ax3.scatter(X,Y,s=0.1,c=color)
                        elif xx == 2 and pp ==2 and zz == 1 and nn==2:
                            color = "y"
                            ax4.scatter(X,Y,s=0.1,c=color)
                        elif xx == 1 and pp ==1 and zz == 2 and nn==2:
                            color = "c"
                            ax5.scatter(X,Y,s=0.1,c=color)
                        elif xx == 2 and pp ==1 and zz == 2 and nn == 2:
                            color = "0.6"
                            ax6.scatter(X,Y,s=0.1,c=color)
                        elif xx == 1 and pp ==2 and zz == 2 and nn == 2:
                            color = "k"
                            ax7.scatter(X,Y,s=0.1,c=color)
                        elif xx == 2 and pp ==2 and zz == 2 and  nn==2:
                            color = "m"
                            ax8.scatter(X,Y,s=0.1,c=color)
                        


        end_time = time.time()
        print("実行時間:"+str(end_time-start_time)+"秒")
        # Save the graph to a fig file
        plt.savefig("fig/"+str(count)+".png")
            
    if __name__:
        main()     
