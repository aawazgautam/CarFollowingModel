# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 19:52:11 2022

@author: Awaz
"""
""" TEAM GREEN"""

"""*** NOTICE: This code is entirely based on the 
Zhao, X.M. and Gao, Z.Y., 2005. A new car-following model: full veloc-ity and acceleration difference model. 
The European Physical Journal B-Condensed Matter and Complex Systems, 47(1), pp.145-150.   ***"""

###It is an effort to solve the given question depending upon Optimal Velocity Method.

# Conversion factor

vcon= 0.44704  #Factor to convert the velocity from mph to m/s.
dcon= 0.3048   #Factor to convert the distance ft to m
pi=22/7

# Constants that remain unchanged
"""The sensitivity of the driver."""
kai = 0.85  # /s

# step 1 : Mentioning the boundary conditions given in Question
# Initial velocity of truck
# Velocities are given in mph, which are converted to m/s
vo_tr= 55*vcon                   # m/s   It is the initial velocity of the truck.
vo_suv= 57*vcon                  # m/s   It is the initial velocity of the SUV.
vo_sed= 54*vcon                  # m/s   It is the initial velocity of the sedan.
vo_car= 55*vcon                  # m/s   It is the initial velocity of the car.



"""We are considering the position of Sedan car that is at the last as (0,0)"""
# Distance are given in ft, which are converted to m
xo_tr= (10+15+25+17+30+80)*dcon   # m    It is the initial position of the truck.
xo_suv= (10+15+25+17)*dcon        # m    It is the initial position of the SUV.      
xo_sed= (10+15)*dcon              # m    It is the initial position of the sedan.
xo_car= 0*dcon                    # m    It is the initial position of the car.

# lengths of the vehicles given in ft, which are converted to m
l_tr= (80)*dcon                   # m    It is the initial position of the truck.
l_suv= (17)*dcon                  # m    It is the initial position of the truck.
l_sed= 15*dcon                    # m    It is the initial position of the truck. 
l_car= 15*dcon                    # m    It is the initial position of the truck.

"""The Truck accelerates to a speed of 60 mph for 5 minutes and then reaches to a speed of 
40mph in next 40 seconds."""
"""The total time frame considered is 5 min 40 seconds i.e 340 seconds"""
T= 340 #seconds  Total time.

""""For this simulation we calcualte values at every 1 seconds."""
del_T = 1 #sec

# step 2: Import libraries.
import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt


# step 2: Creating a function to calculate OV(optimal velocity). 
def Vs(s):
    v1= 6.75  #m/s
    v2= 7.91  #m/s
    c1= 0.13  #per m
    c2= 1.57  #dimensionless
    vs = v1+v2*math.tanh(math.radians((c1*s - c2)))
    
    return (vs)

# step 3: For loop to calculate position, velocity, acceleration, 
N=int(T/del_T)
lst=[1]*N
time=np.array(lst,dtype='int')
"""Creating an array of the size(=no of points were kinematics of the system is calculated."""

# Arrays of acceleration for truck, suv, sedan and car
a_trk=np.array(lst,dtype='float64')
a_suv=np.array(lst,dtype='float64')
a_sed=np.array(lst,dtype='float64')
a_car=np.array(lst,dtype='float64')

# Arrays of velocity for truck, suv, sedan and car
v_trk=np.array(lst,dtype='float64')
v_suv=np.array(lst,dtype='float64')
v_sed=np.array(lst,dtype='float64')
v_car=np.array(lst,dtype='float64')

# Arrays of position for truck, suv, sedan and car
x_trk=np.array(lst,dtype='float64')
x_suv=np.array(lst,dtype='float64')
x_sed=np.array(lst,dtype='float64')
x_car=np.array(lst,dtype='float64')

# Arrays of netto distance for truck, suv, sedan and car
s_suv=np.array(lst,dtype='float64')
s_sed=np.array(lst,dtype='float64')
s_car=np.array(lst,dtype='float64')

i=0      # Initializing a value to a variable to be used in for loop.
r= range(0, N-1)

"""Assigning given boundary conditions as the first element of respective arrays""" 
v_trk[0] = (vo_tr)               #Initial velocity of truck.
x_trk[0] = (xo_tr)               #Initial position of truck.

v_suv[0] = vo_suv                          #Initial velocity of SUV.
x_suv[0] = xo_suv                          #Initial position of SUV.
s_suv[0] = x_trk[0] - x_suv[0] - l_tr      #Initial net to distance of SUV to truck.

v_sed[0] = vo_sed                          #Initial velocity of sedan.
x_sed[0] = xo_sed                          #Initial position of sedan.
s_sed[0] = x_suv[0] - x_sed[0] - l_suv     #Initial net to distance of sedan to SUV.

v_car[0] = vo_car                          #Initial velocity of car.
x_car[0] = xo_car                          #Initial position of car.
s_car[0] = x_sed[0] - x_car[0] - l_car     #Initial net to distance of car to sedan.

lst[0]=1
for i in r:
    if i<300:
        v_trk[i+1] = v_trk[0] + 5*vcon*(i+1)/300  # Given condition for initial acceleration upto 5 min
        
    else:
        v_trk[i+1] = v_trk[300] + 0.5*vcon*(300-i) # Given condition for deceleration of 40 seconds after initial 5 min
    
    x_trk[i+1] = x_trk[i] + (v_trk[i]+v_trk[i+1])/2*del_T             # Calculating the position of truck as the sum of previous position plus average velocity times time.
    
    a_trk[i] =(v_trk[i+1]-v_trk[i]) / del_T
    
    v_suv[i+1] = v_suv[i] + kai*del_T*(Vs(s_suv[i]) - v_suv[i])       # Calculating the velocity of SUV using OV function.
    
    x_suv[i+1] = x_suv[i] + (v_suv[i]+v_suv[i+1])/2*del_T             # Calculating the position of SUV.
    
    s_suv[i+1] = x_trk[i+1] - x_suv[i+1] - l_tr                       # Calculating netto distance from SUV to truck.   
    
    a_suv[i] =(v_suv[i+1]-v_suv[i]) / del_T
    
    v_sed[i+1] = v_sed[i] + kai*del_T*(Vs(s_sed[i]) - v_sed[i])       # Calculating the velocity of sedan.
    
    x_sed[i+1] = x_sed[i] + (v_sed[i]+v_sed[i+1])/2*del_T             # Calculating position of sedan.
    
    s_sed[i+1] = x_suv[i+1] - x_sed[i+1] - l_suv                      # Calculating netto distance from sedan to SUV
     
    a_sed[i] =(v_sed[i+1]-v_sed[i]) / del_T
    
    v_car[i+1] = v_car[i] + kai*del_T*(Vs(s_car[i]) - v_car[i])       # Calculating the velocity of car.
    
    x_car[i+1] = x_car[i] + (v_car[i]+v_car[i+1])/2*del_T             # Calculating position of car.
    
    s_car[i+1] = x_sed[i+1] - x_car[i+1] - l_sed                      # Calculating netto distance from car to sedan.
     
    a_car[i] =(v_car[i+1]-v_car[i]) / del_T
          
    lst[i+1]=i+2    # This variable is stored for the purpose of using it as time axis in graph.

v_suv[0]

"""Plotting the values calculated in the form of Graph"""
## Creating pandas series of the position variables of truck, suv, sedan and car.
a = pd.Series(x_trk)  
b = pd.Series(x_suv)
c = pd.Series(x_sed)
d = pd.Series(x_car)
## Creating pandas series of the time variable
t = pd.Series(lst)

"""INPUT FOR PLOT OF POSITION VS TIME"""
plt.plot(t,a, label = "Truck") 
plt.plot(t,b, label = "SUV") 
plt.plot(t,c, label = "Sedan") 
plt.plot(t,d, label = "Car") 
plt.xlabel('Time (seconds)')
plt.ylabel('Position (meter)')
plt.title('Position - Time Relationship')
plt.legend() # legend upper left
plt.show() # show the plot 

"""INPUT FOR PLOT OF VELOCITY VS TIME"""
plt.plot(t,(pd.Series(v_trk) ), label = "Truck") 
plt.plot(t,pd.Series(v_suv) , label = "SUV") 
plt.plot(t,(pd.Series(v_sed) ), label = "Sedan") 
plt.plot(t,(pd.Series(v_car) ), label = "Car") 
plt.xlabel('Time (seconds)')
plt.ylabel('Velocity (meter)')
plt.title('Velocity - Time Relationship')
plt.legend() # legend upper left
plt.show() # show the plot 


"""INPUT FOR PLOT OF ACCELERATION VS TIME"""
plt.plot(t,(pd.Series(a_trk) ), label = "Truck") 
plt.plot(t,pd.Series(a_suv) , label = "SUV") 
plt.plot(t,(pd.Series(a_sed) ), label = "Sedan") 
plt.plot(t,(pd.Series(a_car) ), label = "Car") 
plt.xlabel('Time (seconds)')
plt.ylabel('Position (meter)')
plt.title('Acceleration - Time Relationship')
plt.legend() # legend upper left
plt.show() # show the plot 

""" Creating the CSV files to export the datas of acceleration, velocity, time for the vehicles in the system"""
df = pd.DataFrame({"Time" : lst,"acc_truck" : a_trk, "acc_suv" : a_suv, "acc_sedan" : a_sed, "acc_car" : a_car,"Time" : lst,"velo_truck" : v_trk, "velo_suv" : v_suv, "velo_sedan" : v_sed, "velo_car" : v_car,"Time" : lst,"Position_truck" : x_trk, "Position_suv" : x_suv, "Position_sedan" : x_sed, "Position_car" : x_car})
df.to_csv("submission.csv", index=False)

### Graphical User Interface###

import PySimpleGUI as sg

sg.theme('LightBlue7') #keep things interesting for your user

layout = [[sg.Text('This program gives the position of vehicles  in above Car Following Model.')],
          [sg.Text('Enter the time in seconds(0 - 340).'),sg.Input(key='-IN-',size=(10,1))],
          [sg.Text('Position of car:', key='-Out-')],
          [sg.Button('Show'), sg.Exit()]]

window = sg.Window('F2C Converter',layout )
while True:
    event, values = window.read()
    if event == 'Show':
        F = int(values["-IN-"])
        C = x_car[F]
        C = round(C,2)
        strc = "The position of car is:" + str(C)
        ""
                
        window["-OUT-"].update(strc)
        print(event.values)
    if event == sg.WIN_CLOSED or event == 'Exit':
        break

window.close()


