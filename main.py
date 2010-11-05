"""
Author: Isaac Goodin
Date Created: 1/12/2010
Last Updated: 8/31/2010

Program to interface with the Neuron.py class and allow simple configuration and execution of the simulation. 
"""

from neuron import *

file  = open("neuron.txt",'w')

steps = 100000 #Total Time(Ms)
skip  = 10000  #Time to Skip(Ms)


#Creating Neuron Objects
N1 = Neuron("Neuron 1")
N2 = Neuron("Neuron 2")
N3 = Neuron("Neuron 3")

#Making a list of all the Neurons
Neuron_list = [N1,N2,N3]

#Setting Initial values
N1.X =[-21.0,0.11,0.41,0.11]
N2.X =[-22.0,0.12,0.42,0.12]
N3.X =[-23.0,0.13,0.43,0.13]

#Changing Neuron Parameters
N1.gsr = 0.42
N2.gsr = 0.40
N3.gsr = 0.38

#Setting Neuron Inputs and correlation strength
N2.Input([[N1,0.0035]])
N3.Input([[N2,0.003]])

#Run the simulation and output data 
Run_rk4(file,steps,skip,Neuron_list)
