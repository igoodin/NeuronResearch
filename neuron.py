"""
Author: Isaac Goodin
Date Created: 1/12/2010
Last Updated: 8/31/2010

Simulation of Neuron firing and synchronization using the Hodgkin-Huxley model of Neurons.
"""

import math 
class Neuron:

	#Neuron object based on the Hodgkin Huxley model of a neuron

    def __init__(self,type="Empty"):
        #Ignoring type for now
        self.type = type
        self.input = None

		#Neuron parameters from experimental values
        self.phi   = 0.124
        self.rho   = 0.607
        self.theta = 0.17
        self.nu    = 0.012
        self.cm    = 1.0
        self.gleak = 0.1
        self.gna   = 1.5
        self.gk    = 2.0
        self.gsd   = 0.25
        self.gsr   = 0.20
        self.tk    = 2.0
        self.tsd   = 10.0
        self.tsr   = 20.0
        self.sna   = 0.25
        self.sk    = 0.25
        self.ssd   = 0.09
        self.vona  = -25.0
        self.vok   = -25.0
        self.vosd  = -40.0
        self.vleak = -60.0
        self.vna   = 50.0
        self.vk    = -90.0
        self.vsd   = 50.0
        self.vsr   = -90.0
                    
    def Input(self,input_neurons):
        self.input = input_neurons
            
    def Calc_input(self):
		#calculates the total input for the neuron
        input_neurons = self.input
        if self.input != None:
                i=0
                val=0.0
                for element in input_neurons:
                        val +=element[1]*(element[0].X[0] - self.X[0])
                        i+=1
                return val
        else:
                return 0.0
                    
    def DE(self,Z,Eq,Ic):
        if(Eq==0):
                aNai = 1.0/(1.0 + math.exp(-self.sna * (Z[0] - self.vona)))
        
                Ileak = self.gleak*(Z[0]-self.vleak);
                INa   = self.rho * self.gna * aNai * (Z[0] - self.vna)
                IK    = self.rho * self.gk  * Z[1] * (Z[0] - self.vk)
                Isd   = self.rho * self.gsd * Z[2] * (Z[0] - self.vsd)
                Isr   = self.rho * self.gsr * Z[3] * (Z[0] - self.vsr)
        
                val = (-Ileak - INa - IK - Isd - Isr + Ic - 1.0)/self.cm
        elif(Eq==1):
                aKi = 1.0/(1.0 + math.exp(-self.sk * (Z[0] - self.vok)))
                val = self.phi * (aKi - Z[1]) /self.tk
        elif(Eq==2):
                asdi = 1.0/(1.0 + math.exp(-self.ssd * (Z[0] - self.vosd)))
                val = self.phi * (asdi - Z[2]) / self.tsd
        elif(Eq==3):
                Isd   = self.rho*self.gsd*Z[2]*(Z[0]-self.vsd)
                val = -self.phi * (self.nu * Isd + self.theta * Z[3]) / self.tsr
        return val

    def rk4(self):
		#Fourth order Runge Kutta Method
	
        #Get the input sigals from the other Neurons
        Ic = self.Calc_input()

		#step size
        h=0.1
        
        k1=[0.0,0.0,0.0,0.0]
        k2=[0.0,0.0,0.0,0.0]
        k3=[0.0,0.0,0.0,0.0]
        k4=[0.0,0.0,0.0,0.0]
        Xn=[0.0,0.0,0.0,0.0]

        #These loops are split apart because Xn must be updated before the next K calculation
        for i in range(4):
                # Calculates K1 and updates variable array
                k1[i] = h*self.DE(self.X,i,Ic)
                Xn[i] = self.X[i] + k1[i]/2.

        for i in range(4):
                # Calculates K2 and updates variable array
                k2[i] = h*self.DE(Xn,i,Ic)
                Xn[i] = self.X[i] + k2[i]/2.

        for i in range(4):
                # Calculates K3 and updates variable array
                k3[i] = h*self.DE(Xn,i,Ic)
                Xn[i] = self.X[i] + k3[i]

        for i in range(4):
                # Calculates K4 and updates variable array
                k4[i] = h*self.DE(Xn,i,Ic)
                self.X[i] = self.X[i] + (k1[i]+2.*k2[i]+2.*k3[i]+k4[i])/6.

def Run_rk4(file,steps,skip,Nlist):
	#runs the Rk4 method and outputs data for plotting
    t=0
    for i in range(steps):         
        #skip initial time
        while(t<skip):
            for neuron in Nlist:
                    neuron.rk4()
            t+=1
        file.write(str(i*0.1)+ ' ')
        for neuron in Nlist:
            neuron.rk4()
            file.write(str(neuron.X[0])+' ')
        file.write('\n')
    file.close()
