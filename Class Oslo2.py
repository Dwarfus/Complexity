# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:05:38 2017
@author: pk3014
CN project
"""
# 2,4,8,16,32,64,128,256,512,1024,2048,4096
import numpy as np
import matplotlib.pyplot as plt

class Oslo:
    def __init__(self, L=256, treshold = (1,2), prob=0.5, nruns = 50):
        self.z=np.zeros((L)) # THe values of slopes for each point
        self.ztres = np.ones((L)) # The treshold values for each point

        self.firsth = 0 # THis keeps the height of the first instance
        self.heights = [0]
        self.nruns = nruns # This is the number of times the code should run
        self.s = [] # avalanche size for each grain added
        self.counter = 0 # used to measure the time in units of total grains added. 
        self.lost = 0 # The count of grains leaving the system
        self.steady = 0 # check whether we are in steady state, this might not be needed after while
        self.treshold = treshold # The values of possible treshold slopes
        self.L = L # The size of the system
        self.p = prob # The probability of the first treshold
        self.avalanche = 0 # the counter of size of avalanches for one grain added
        i=0
        while (i<self.L):
            # this sets the initial treshold values
            self.ztres[i]= np.random.choice(treshold, p=(self.p,1-self.p))
            i+=1
        self.call()    
            
    def drive(self):
        """
        THis method is the drive algorithm. E.g. it adds one point to z0
        """
        self.z[0]+=1
        self.firsth +=1
       # print ("added")
        self.avalanche=0
        self.counter +=1
        self.relax=1
        
        
    def relaxation(self):
        """"
        This method is the relaxation one. It checks for any z>ztres and carries out an appropriate action
        It is important to note that always scan through the whole array as by relaxin i, you can now have i-1 
        to topple even though it was previously stable.
        It is important to run the relaxation as long as there is something to relax
        """
        
        while self.relax==1:        
            self.toppling = 0
            self.toppling = np.where(self.z>self.ztres)
        
            if self.toppling[0].size>0: # checking if there is something to relax
                integ = self.toppling[0][0]
                if integ==0:
                    self.z[0]-=2
                    self.z[1]+=1
                    self.firsth-=1
                
                elif integ>0 and integ<self.L-1:
                    self.z[integ]-=2
                    self.z[integ-1]+=1
                    self.z[integ+1]+=1
                
                elif integ ==self.L-1:
                    self.z[integ]-=1
                    self.z[integ-1]+=1
                    self.lost +=1
                
                self.ztres[integ] = np.random.choice(self.treshold, p=(self.p,1-self.p))
                #print "topple"
                self.avalanche+=1
                

            else:
            
                self.s.append(self.avalanche)
                self.heights.append(self.firsth)
                self.relax = 0
          
         
                
                    
    def moveheights(self):
        
        # This method calculates the heights from the value of self.firsth and the self.z        
        self.W= 50
        
        leng = len(self.heights)-self.W
        i = self.W
        self.moh=np.zeros((1,leng))
        while i< (leng):
            self.moh[0][i] = (np.mean(self.heights[i-self.W:i+self.W+1]))
            
            i+=1     

    """
    def task2c(self, t0=50, T=len(a.heights)-1):
        self.averageh = np.sum(self.heights[t0+1:t0+T+1])/T
        squares = [x**2 for x in self.heights[t0+1:t0+T+1]]
        self.sqaverageh = (np.sum(squares))/T
        self.sigma = np.sqrt(self.sqaverageh-self.averageh**2) 
        print(self.averageh, self.sqaverageh, T, self.sigma)"""
                    


        

        
        
    def call(self):
        self.run=0

        while self.run<self.nruns:     
            self.drive()
            self.relaxation()
            self.run+=1
        #self.task2c()
        #self.moveheights()    
        #self.heights()
        #self.avalancheplot()
        
#a = Oslo(16,(1,2), 0.5,500 )
        

class Results:
    def __init__(self, L=[], nruns=[], steady=[]):
        """
        This class contains the different function used to call the Oslo.
        Each method then calculates whatever it needs.         
        """
        self.L = L
        self.nruns = nruns
        self.steady = steady
        self.difruns = len(self.L)
        self.results = np.array([0.0]*self.difruns)
        
        self.call2b()
        
        
    
    
    def call2(self):
        """This is used to call for task 2a
        finds the value of steady state height and 
        the first occurance of that value.
        It plots the height vs t and highlight the above point"""
        i=0
        while(i<self.difruns):
            a = Oslo(self.L[i],(1,2), 0.5,self.nruns[i] )
            res1=a.heights
            legend1=("L= %i" %self.L[i])
            rs1=np.mean(res1[self.steady[i]:])
            pos1 = np.where(res1>rs1)
            print(self.L[i],rs1, pos1[0][0])
            scale = list(range(0,self.nruns[i]+1))
            plt.plot(scale, res1, label=legend1)
            plt.axhline(y=rs1)
            plt.axvline(x=pos1[0][0])

            i+=1
            
        plt.xlabel("Grains added")
        plt.ylabel("Height of the system")
        plt.legend(loc = 'upper left')
        plt.show()


        
    def call2b(self):
        """This is used to call for task 2b
        Implementing data collaps"""
        i=0
        self.w=25
        plt.figure(1)
        plt.title("Data collapse for moving average")
        while(i<self.difruns):
            a = Oslo(self.L[i],(1,2), 0.5,self.nruns[i] )
            l = float(self.L[i])
            self.heights= a.heights
            j = 0
            leng = len(self.heights)-2*self.w
            self.moh = np.zeros((1, leng))
            while j<(leng):
                self.moh[0][j]= np.mean(self.heights[j:j+2*self.w+1])
                j+=1
                
                
            
            self.res = [float(y/l) for y in self.moh[0]]
            legend1=("L= %i" %self.L[i])
            self.scale1 = list(range(self.w,leng+self.w))                        
            self.scale = [float(x) for x in self.scale1]
            self.scale = [float(x/(l**2)) for x in self.scale]            
            plt.loglog(self.scale, self.res,  label=legend1)

            i+=1
        # Now for the highest L I will calculate the gradient
        self.scafit= self.scale[int(self.steady[-1]*0.1):int(self.steady[-1]*0.8)]
        self.p = np.polyfit(np.log10(self.scafit),np.log10(self.res[int(self.steady[-1]*0.1):int(self.steady[-1]*0.8)]),deg=1 )
        print(self.p)
        # Plotting the linear fit - works          
        self.values =[ z**self.p[0]*10**(self.p[1]) for z in self.scafit]
        plt.loglog(self.scafit, self.values, label="line of best fit")
        plt.xlabel("log10(t/(L**2))")
        plt.ylabel("log10(Moving average/L)")
        plt.legend(loc = 0)
        plt.show()
        
  
        
        
        
    def function( L,a,c,w ):
        return a+c(L**(-w))

    def call2c(self):
        i=0
        print("L, average height, std, number of points used,  check normalization")
        self.averages = []
        self.sigmas= []
        while(i<self.difruns):
            
            a = Oslo(self.L[i],(1,2), 0.5,self.nruns[i] )
            
            l = self.L[i]
            t0 = self.steady[i]
            
             
            self.steadyheights = a.heights[t0+1:]
            T = len(self.steadyheights)
            self.averageh = float(np.sum(self.steadyheights)/float(T))
            #squares = [x**2 for x in self.heights[t0+1:t0+T+1]]
            #self.sqaverageh = (np.sum(squares))/T
            #self.sigma = np.sqrt(self.sqaverageh-self.averageh**2) 
            self.sigma2 = np.std(self.steadyheights)
            self.unique = list(set(self.steadyheights))
            self.probabilities = np.zeros((len(self.unique),2))
            j=0
            self.sum=0
            for height in self.unique:
                self.probabilities[j][0] = height
                self.probabilities[j][1] = self.steadyheights.count(height)/len(self.steadyheights)
                self.sum += self.steadyheights.count(height)/len(self.steadyheights)
                j+=1
            
            print(l, self.averageh,  self.sigma2, T, self.sum)
            self.averages.append(self.averageh/l)
            self.sigmas.append(self.sigma2)
            
            i+=1
        
        plt.figure(2)
        plt.plot(np.log10(self.L), np.log10(self.averages), label="Task 2c")
        
        plt.xlabel("log(L)")
        plt.ylabel("log(<h>/l)")
        plt.legend(loc=0)
        plt.show()
        
        """
        self.averageh = np.sum(self.heights[t0+1:t0+T+1])/T
        squares = [x**2 for x in self.heights[t0+1:t0+T+1]]
        self.sqaverageh = (np.sum(squares))/T
        self.sigma = np.sqrt(self.sqaverageh-self.averageh**2) 
        print(self.averageh, self.sqaverageh, T, self.sigma)"""
        

#tresholds in time are[54,226,898,3391,14056,56437]        
b=Results(L=[32], nruns=[6000], steady=[900])  
#b=Results(L=[8,16,32,64,128], nruns=[5100,5300, 5900,8400,20000], steady=[100,300,900, 3400,15000]) 
#b=Results(L=[8,16,32,64,128,256,512], nruns=[5000,5000, 5000,6000,20000,80000,300000], steady=[100,350,1000,3700,15500, 60000])
#b=Results(L=[8,16,32,64,128,256], nruns=[5000,5000, 6000,7000,20000,80000], steady=[100,350,1000,3700,15500, 60000])
