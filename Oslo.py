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
        self.heights = []
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
        
        
    def relaxation(self):
        """"
        This method is the relaxation one. It checks for any z>ztres and carries out an appropriate action
        It is important to note that always scan through the whole array as by relaxin i, you can now have i-1 
        to topple even though it was previously stable.
        It is important to run the relaxation as long as there is something to relax
        """
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
                
        #print self.z  
        #print self.ztres
        #self.toppling = 0
        #self.toppling = np.where(self.z>self.ztres)        
        #if self.toppling[0].size>0:
            self.relaxation()
        else:
            
            self.s.append(self.avalanche)
            self.heights.append(self.firsth)
           # print (self.z, self.firsth)
            #print (self.avalanche)
           

                
    """                
    def heights(self):
        
        # This method calculates the heights from the value of self.firsth and the self.z        
        
        helpH=self.firsth
        i = 0
        for value in self.z:
            self.height[i] = helpH            
            helpH -=self.z[i]            
            i+=1"""
            

    def avalancheplot(self):
        length = len(self.s)+1
        scale = list(range(1, length))
        plt.stem(scale, self.s)   
        plt.show()
        
    def heighplot(self):
        scale = list(range(1, len(self.heights)+1))
        plot1 = plt.plot(scale, self.heights)
        plt.xlabel("Grains added")
        plt.ylabel("Height of the system")
        plt.figlegend(plot1,('label1',), 'upper right')
        plt.show()
        
        
    def call(self):
        self.run=0

        while self.run<self.nruns:     
            self.drive()
            self.relaxation()
            self.run+=1

        self.heighplot()    
        #self.heights()
        #self.avalancheplot()
        
        
a=Oslo(8,(1,2),0.2,1000)

