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
                
        #print self.z  
        #print self.ztres
        #self.toppling = 0
        #self.toppling = np.where(self.z>self.ztres)        
        #if self.toppling[0].size>0:

                #self.relaxation()
            else:
            
                self.s.append(self.avalanche)
                self.heights.append(self.firsth)
                self.relax = 0
           # print (self.z, self.firsth)
            #print (self.avalanche)
           
         
                
                    
    def moveheights(self):
        
        # This method calculates the heights from the value of self.firsth and the self.z        
        self.W= 25
        
        leng = len(self.heights)-self.W
        i = 25
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
        #self.task2c()
       # self.moveheights()    
        #self.heights()
        #self.avalancheplot()
        
#a = Oslo(16,(1,2), 0.5,500 )
        

class Results:
    def __init__(self, L=[], nruns=[], steady=[]):
        self.L = L
        self.nruns = nruns
        self.steady = steady
        self.difruns = len(self.L)
        self.results = np.array([0.0]*self.difruns)
        
        self.call2b()
        
        
    
    
    def call2(self):
        """This is used to call for task 2a"""
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
           # print("L = %i, htres=%i, at ttres=%i", %(self.L[i],rs1, pos1 ))
            i+=1
            
        plt.xlabel("Grains added")
        plt.ylabel("Height of the system")
        plt.legend(loc = 'upper left')
        plt.show()
        
    def call2b(self):
        """This is used to call for task 2b
        Implementing data collaps for moving average"""
        self.W= 25
        i=0
        plt.figure(1)
        plt.title("Data collapse for moving average")
        while(i<self.difruns):
            a = Oslo(self.L[i],(1,2), 0.5,self.nruns[i] )
            self.heights = a.heights
        
            leng = len(self.heights)-self.W
            k = 25
            self.moh=np.zeros((1,leng))
            while k< (leng):
                self.moh[0][k] = (np.mean(self.heights[k-self.W:k+self.W+1]))
            
                k+=1 
                
                
            l = self.L[i]
            res1=[self.moh[0]/l]            
            legend1=("L= %i" %self.L[i])
            scale = list(range(0,len(res1[0])))
            scale = [x /(l**2) for x in scale]
            print(np.shape(res1[0]), np.shape(scale))
            
            plt.plot(np.log10(scale), np.log10(res1[0]), label=legend1)

            i+=1
            
        plt.xlabel("log(t(grains)/(L**2))")
        plt.ylabel("log(Moving average/L)")
        plt.legend(loc = 'lower right')
        plt.show()
        
    def call2c(self):
        i=0
        
        while(i<self.difruns):
            
            a = Oslo(self.L[i],(1,2), 0.5,self.nruns[i] )
            
            l = self.L[i]
            t0 = self.steady[i]
            
             
            steadyheights = a.heights[t0+1:]
            T = len(steadyheights)
            self.averageh = np.sum(steadyheights)/T
            #squares = [x**2 for x in self.heights[t0+1:t0+T+1]]
            #self.sqaverageh = (np.sum(squares))/T
            #self.sigma = np.sqrt(self.sqaverageh-self.averageh**2) 
            self.sigma2 = np.std(steadyheights)
            self.unique = list(set(steadyheights))
            self.probabilities = np.zeros((len(self.unique),2))
            j=0
            self.sum=0
            for height in self.unique:
                self.probabilities[i][0]= height
                self.probabilities[i][1] = steadyheights.count(height)/len(steadyheights)
                self.sum += steadyheights.count(height)/len(steadyheights)
                j+=1
            print("L, average height, std, number of points used,  check normalization")
            print(l, self.averageh,  self.sigma2, T, self.sum)
            
            i+=1
            
        """
        self.averageh = np.sum(self.heights[t0+1:t0+T+1])/T
        squares = [x**2 for x in self.heights[t0+1:t0+T+1]]
        self.sqaverageh = (np.sum(squares))/T
        self.sigma = np.sqrt(self.sqaverageh-self.averageh**2) 
        print(self.averageh, self.sqaverageh, T, self.sigma)"""
        
        
b=Results(L=[32], nruns=[5000], steady=[350])  
#b=Results(L=[8,16,32,64], nruns=[5000,5000, 5000,6000], steady=[200,300,1000]) 
#b=Results(L=[8,16,32,64,128,256,512], nruns=[5000,5000, 5000,6000,20000,80000,300000], steady=[100,350,1000,3700,15500, 60000])
#b=Results(L=[8,16,32,64,128,256], nruns=[5000,5000, 6000,7000,20000,80000], steady=[100,350,1000,3700,15500, 60000])
"""        
a=Oslo(8,(1,2),0.5,65000)
results1 = a.heights
legend1= ("L=8")
rs1=np.mean(results1[200:])
pos1 = np.where(results1>rs1)
print(rs1, pos1[0][0])
a=Oslo(16,(1,2),0.5,65000)
results2 = a.heights
legend2= "L=16"
rs2=np.mean(results2[300:])
pos2 = np.where(results2>rs2)
print(rs2, pos2[0][0])
a=Oslo(32,(1,2),0.5,65000)
results3 = a.heights
legend3= "L=32"
rs3=np.mean(results3[1000:])
pos3 = np.where(results3>rs3)
print(rs3, pos3[0][0])
a=Oslo(64,(1,2),0.5,65000)
results4 = a.heights
legend4= "L=64"
rs4=np.mean(results4[4000:])
pos4 = np.where(results4>rs4)
print(rs4, pos4[0][0])
a=Oslo(128,(1,2),0.5,65000)
results5 = a.heights
legend5= "L=128"
rs5=np.mean(results5[40000:])
pos5 = np.where(results5>rs5)
print(rs5, pos5[0][0])
a=Oslo(256,(1,2),0.5,65000)
results6 = a.heights
legend6= "L=256"
rs6=np.mean(results6[58000:])
pos6 = np.where(results6>rs6)
print(rs6, pos6[0][0])
scale = list(range(0, 65001))
plot1 = plt.plot(scale, results1, label = legend1)
plot2 = plt.plot(scale, results2, label = legend2)
plot3 = plt.plot(scale, results3, label = legend3 )
plot4 = plt.plot(scale, results4, label = legend4)
plot5 = plt.plot(scale, results5, label = legend5)
plot6 = plt.plot(scale, results6, label = legend6)
plt.xlabel("Grains added")
plt.ylabel("Height of the system")
plt.legend(loc = 'upper left')
#plt.figlegend((plot1, plot2, plot3, plot4),(legend1,legend2,legend3,legend4), 'upper right')
plt.axhline(y=rs1)
plt.axhline(y=rs2)
plt.axhline(y=rs3)
plt.axhline(y=rs4)
plt.axhline(y=rs5)
plt.axhline(y=rs5)
plt.axhline(y=rs6)
plt.axvline(x=pos1[0][0])
plt.axvline(x=pos2[0][0])
plt.axvline(x=pos3[0][0])
plt.axvline(x=pos4[0][0])
plt.axvline(x=pos5[0][0])
plt.axvline(x=pos6[0][0])
"""
