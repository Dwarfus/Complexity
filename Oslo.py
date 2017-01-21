# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:05:38 2017

@author: pk3014
CN project
"""
# 2,4,8,16,32,64,128,256,512,1024,2048,4096
import numpy as np

class Oslo:
    def __init__(self, L=256, treshold = (1,2), prob=0.5, nruns = 50):
        self.z=np.zeros((L))
        self.ztres = np.ones((L))
        self.height = np.zeros(L)
        self.firsth = 0
        self.nruns = nruns
        self.s = np.zeros((nruns)) # avalanche size
        self.counter = 0 # used to measure the time in units of total grains added. It will be useful when assesing the steady state
        self.treshold = treshold
        self.L = L
        self.p = prob
        self.avalanche = 0 # the counter of size of avalanches
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
        print "added"
        self.avalanche=0
        
        
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
            self.s[self.run] = self.avalanche
            print self.z, self.firsth
            print self.avalanche
            
    def heights(self):
        """
        This method calculates the heights from the value of self.firsth and the self.z        
        """
        helpH=self.firsth
        i = 0
        for value in self.z:
            self.height[i] = helpH            
            helpH -=self.z[i]            
            i+=1
            
            
        
        
    def call(self):
        self.run=0
        while self.run<self.nruns:     
            self.drive()
            self.relaxation()
            self.run+=1
            
        self.heights()
        
        
a=Oslo(16,(1,2),0,200)

