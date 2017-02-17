# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:05:38 2017
@author: pk3014
CN project
"""
# 2,4,8,16,32,64,128,256,512,1024,2048,4096
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from log_bin_CN_2016 import log_bin
import pickle

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
        self.drop = []
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
        This method is the drive algorithm. E.g. it adds one point to z0
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
                self.drop.append(self.lost)
                self.heights.append(self.firsth)
                self.relax = 0
                self.lost = 0
          
         
                

        
        
    def call(self):
        self.run=0

        while self.run<self.nruns:     
            self.drive()
            self.relaxation()
            self.run+=1
        
#a = Oslo(8,(1,2), 0.5,1000 ) # as threshold is definitely under 100 and I want million runs after it
        
#  Test this first
#with open('L8.pkl', 'wb') as output:    
#    pickle.dump(a, output, pickle.HIGHEST_PROTOCOL) 
        
#with open('company_data.pkl', 'rb') as input:
    #company1 = pickle.load(input)
        
#np.save('outfilename', array)
#np.load()
        
def func( L,a,c,w):
    return a+c*(L**(-w))        

#['L8','L16','L32','L64','L128','L256','L512']
class Analysis:
    def __init__(self, data=['8','16','32','64','128','256','512']): 
            """
            Loads the data from pickled objects.
            The objects contain:
                heights as .heights
                avalanche sizes as .s
                drop sizes as .drop
            The data was collected so that the last 10^6 is always the steady state
            e.g. from previous knowledge of tc the number of grains added was picked 
            so that we can use last 10^6 values as steady state for all systems.
            """           
            self.objects = [] 
            self.data = data
            self.L = [float(x) for x in data]
            for value in data:
                with open('L'+value+'.pkl', 'rb') as input:
                    self.objects.append( pickle.load(input))
                    
            #self.task2a()
            #self.task2b()
            self.task2c()
                    
    def task2a(self):
        """
        Plots the height vs t, either as a normal or loglog plot.
        Also can plot the lines of the steady state height and the t at which
        this height was reached for the first time.        
        """
        
        i=0
        plt.figure(1)
        plt.title("Height as a function of time")    
        #plt.title("Log plot of h(t)")
        for obj in self.objects:
            height = obj.heights
            
            legend1=('L='+self.data[i])
            # find average height of recurrent state
            # The data was collected so that at least last 10^6 points are steady states
            rs=np.mean(height[-1000000:])
            pos1 = np.where(height>rs)
            print(self.data[i],rs, pos1[0][0])
            scale = list(range(0,len(height)))
            plt.plot(scale, height, label=legend1)
            #plt.loglog(scale, height, label=legend1)
            
            #plt.axhline(y=rs) # the average height line
            #plt.axvline(x=pos1[0][0]) # the steady state reached
            plt.grid()

            i+=1
        self.lastpos = pos1[0][0]    
        plt.xlabel("t")
        plt.ylabel("h(t,L)")
        plt.legend(loc = 0)
        plt.show()
        
    def task2b(self):
        """
        The data collapse for moving average.
        Plots h/L vs t/L^2 which results in power low
        The moving average is calculated from t to t+2w, so for j=0 it 
        effectively calculates the moving average around t+w.
        the x axis is calculated appropriately as it starts at w
        """
        i=0
        self.w=25
        plt.figure(2)
        plt.title("Data collapse for moving average")
        for obj in self.objects:
            
            l = float(self.data[i])
            height= obj.heights
            j = 25
            leng = len(height)-2*self.w
            self.moh = np.zeros(( leng))
            #calculates the values of htylda
            while j<(leng):
                # h is here scaled by L to get data collapse
                self.moh[j]= float(np.mean(height[j:j+2*self.w+1])/l)
                j+=1
                
                
            legend1=('L='+self.data[i])
            self.scale1 = list(range(self.w,leng+self.w))                        
            # data collapse of the time
            self.scale = [float(x/(l**2)) for x in self.scale1]   
            
            plt.loglog(self.scale, self.moh,  label=legend1)

            i+=1
        # Now for the highest L I will calculate the gradient
        
        # If running only 2b for L = 512 as highest L uncomment the line below        
        self.lastpos= 225396   
        # Below it calculates the slope from a certain interval on the transient line
        self.scafit= self.scale[int(self.lastpos*0.2):int(self.lastpos*0.9)]
        self.p = np.polyfit(np.log10(self.scafit),np.log10(self.moh[int(self.lastpos*0.2):int(self.lastpos*0.9)]),deg=1 )
        print("The power law coeff. is: %s" %self.p)
        # Plotting the linear fit - in most plots it is not visible, so no need to plot it anyway.
        # It served as good check that the fit is correct, if wanted, uncomment the two lines below.
        #self.values =[ z**self.p[0]*10**(self.p[1]) for z in self.scafit]
        #plt.loglog(self.scafit, self.values)
        plt.xlabel("t/(L^2)")
        plt.ylabel("Moving average/L")
        plt.legend(loc = 0)
        plt.grid()
        plt.show()
        
                

    def task2c(self):
        """
        Also includes 2d                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        """
        i=0
        print("L, average height, std, number of points used,  check normalization, av.slope")
        self.averages = []
        self.sigmas= []
        self.probfinal=[]
        self.heightfinal = []
        self.averh=[]
        for obj in self.objects:
            
            
            l = float(self.data[i])
            
             
            height = obj.heights[-1000000:]
            T = float(len(height))
            self.averageh = float(np.mean(height))
            # It was checked that np.std gives same value as the full expression given in instructions        
            self.sigma = np.std(height)
            self.unique = list(set(height))
            self.probabilities = []
            self.heights=[]
            j=0
            self.sum=0
            for value in self.unique:
                self.heights.append( value)
                self.probabilities.append(float(height.count(value)/T))
                self.sum += height.count(value)/T
                
            
            print(l, self.averageh,  self.sigma, T, self.sum)
            self.probfinal.append(self.probabilities)
            self.heightfinal.append(self.heights)
            self.averh.append(self.averageh)
            self.averages.append(self.averageh/l)
            self.sigmas.append(self.sigma)
            
            i+=1
        
        plt.figure(3)
        plt.title("Scaled height vs L")
        plt.plot(self.L, self.averages, 'o')        
        plt.xlabel("L")
        plt.ylabel("<h>/L ")
        #plt.legend(loc=0)
        plt.grid()
        plt.show()
        
        
        plt.figure(4)
        plt.grid()
        plt.loglog(self.L, self.sigmas, 'o', label="Sigma values")
        
        plt.xlabel("L")
        plt.ylabel("Sigma")
        plt.legend(loc=0)
        plt.show()
        
        sig = np.polyfit(np.log10(self.L), np.log10(self.sigmas), deg=1, full=True)
        print("the sigma fit is is", sig[0], "residuals=", sig[1])
        
        # A simple curve fit gives an idea of a0 and omega
        self.popt, self.pcov= curve_fit(func, self.L, self.averages, p0=[1.7,-0.5,0.5])
        self.perr = np.sqrt(np.diag(self.pcov))
        print("a0,a0*a1, omega1 are (with errors)",self.popt,self.perr)
        
        #this part includes the different method for evaluating low size corrections                
        plt.figure(5)
               
        self.avalues = np.linspace(1.726,1.76,51 )
        plt.title("Plots for different a0")
        self.res = []
        self.residual=[]
        for value in self.avalues:
            plt.loglog( self.L,abs(1-self.averages/value), label=value)
            fit = np.polyfit( np.log10(self.L),np.log10(1-self.averages/value), deg=1, full=True)
            self.res.append(fit[0])
            self.residual.append(fit[1])
         
        minres= np.argmin(self.residual)
        err = (self.avalues[2]-self.avalues[1])/2
        print("The other method gives -omega, a0, a0 error, residual", self.res[minres][0], self.avalues[minres], err, self.residual[minres]) 
        plt.xlabel("L")
        plt.ylabel("1-h/(L*a0)")
        plt.legend(loc= 0)
        plt.grid() 
        plt.show()
        
        plt.figure(3)
        plt.axhline(self.avalues[minres])
        
        # This part is the 2d task        
        j=0
        plt.figure(6)
        plt.grid()
        plt.title("The data collapse of P(h;L)")
        plt.xlabel("(h-<h>)/sigma")
        plt.ylabel("P*sigma")
        while j< len(self.data):
            
            scaledh = [(h-self.averh[j])/self.sigmas[j] for h in self.heightfinal[j]]
            scaledp = [v*self.sigmas[j] for v in self.probfinal[j]]
            plt.plot(scaledh, scaledp, label = 'L='+ self.data[j])
            
            j+=1
        plt.legend(loc=0)
        plt.show()

c= Analysis()


################################################################################################################
#################################################################################################################
#################################################################################################################


                
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
        
        self.call3a()
        
        
    
    
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
            print("%s Done" %(l))
            i+=1
        # Now for the highest L I will calculate the gradient
        self.scafit= self.scale[int(self.steady[-1]*0.15):int(self.steady[-1]*0.75)]
        self.p = np.polyfit(np.log10(self.scafit),np.log10(self.res[int(self.steady[-1]*0.15):int(self.steady[-1]*0.75)]),deg=1 )
        print(self.p)
        # Plotting the linear fit - works          
        self.values =[ z**self.p[0]*10**(self.p[1]) for z in self.scafit]
        plt.loglog(self.scafit, self.values, label="line of best fit")
        plt.xlabel("log10(t/(L**2))")
        plt.ylabel("log10(Moving average/L)")
        plt.legend(loc = 0)
        plt.show()
        
  
        
        
        


    def call2c(self):
        """
        Also includes 2d                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        """
        i=0
        print("L, average height, std, number of points used,  check normalization, av.slope")
        self.averages = []
        self.sigmas= []
        self.probfinal=[]
        self.heightfinal = []
        self.averh=[]
        while(i<self.difruns):
            print("started")
            a = Oslo(self.L[i],(1,2), 0.5,self.nruns[i] )
            self.avslope = np.mean(a.z)
            
            
            l = self.L[i]
            t0 = self.steady[i]
             
            self.steadyheights = a.heights[t0+1:]
            T = float(len(self.steadyheights))
            self.averageh = float(np.mean(self.steadyheights))
            #squares = [x**2 for x in self.heights[t0+1:t0+T+1]]
            #self.sqaverageh = (np.sum(squares))/T
            #self.sigma = np.sqrt(self.sqaverageh-self.averageh**2) 
            self.sigma2 = np.std(self.steadyheights)
            self.unique = list(set(self.steadyheights))
            self.probabilities = []
            self.heights=[]
            j=0
            self.sum=0
            for height in self.unique:
                self.heights.append( height)
                self.probabilities.append(float(self.steadyheights.count(height)/T))
                self.sum += self.steadyheights.count(height)/T
                
            
            print(l, self.averageh,  self.sigma2, T, self.sum, self.avslope)
            self.probfinal.append(self.probabilities)
            self.heightfinal.append(self.heights)
            self.averh.append(self.averageh)
            self.averages.append(self.averageh/l)
            self.sigmas.append(self.sigma2)
            
            i+=1
        
        plt.figure(1)
        plt.plot(self.L, self.averages, 'ro', label="Rescaled average height")
        
        plt.xlabel("L")
        plt.ylabel("<h>/L ")
        plt.legend(loc=0)
        plt.show()
        
        
        plt.figure(3)
        plt.loglog(self.L, self.sigmas, 'ro', label="Sigma values")
        
        plt.xlabel("L")
        plt.ylabel("Sigma ")
        plt.legend(loc=0)
        plt.show()
        sig = np.polyfit(np.log10(self.L), np.log10(self.sigmas), deg=1, full=True)
        print("the sigma fit is is", sig)
        
        self.popt, self.pcov= curve_fit(func, self.L, self.averages, p0=[1.7,-0.5,0.5])
        self.perr = np.sqrt(np.diag(self.pcov))
        print("The curve fit is",self.popt,self.perr)
        #this part includes the different method for evaluating low size corrections        
        
        plt.figure(2)        
        self.avalues = np.linspace(1.7,1.8,21 )
        plt.title("Plots for different a0")
        self.res = []
        self.residual=[]
        for value in self.avalues:
            plt.loglog( self.L,abs(self.averages-value), label=value)
            fit = np.polyfit( np.log10(self.L),np.log10(abs(self.averages-value)), deg=1, full=True)
            self.res.append(fit[0])
            self.residual.append(fit[1])
         
        minres= np.argmin(self.residual)
        
        print("The other method is", self.res[minres], self.avalues[minres]) 
        plt.xlabel("Log10(L)")
        plt.ylabel("Log10(|h/L-a|)")
        plt.legend(loc= 0)
        plt.show()
        
        # This part is the 2d task        
        j=0
        plt.figure(5)
        plt.title("The data collapse of probabilities")
        plt.xlabel("(h-<h>)/sigma")
        plt.ylabel("P*sigma")
        while j< len(self.L):
            
            scaledh = [(h-self.averh[j])/self.sigmas[j] for h in self.heightfinal[j]]
            scaledp = [v*self.sigmas[j] for v in self.probfinal[j]]
            plt.plot(scaledh, scaledp,label=self.L[j])
            
            j+=1
        
        plt.show()
        
        
    def call3a(self):
        """
        This is the function for task 3a and also for 3b if the plot of normal prob is commented out
        This function is ready for use, just need to run it according to the script... Will probably take overnight
        """
        plt.figure(6)
        plt.title("Avalanche size probability")
        plt.xlabel("Avalanche size")
        plt.ylabel("Probability")
        
        i=0
        while(i<self.difruns):
            print("started")
            a = Oslo(self.L[i],(1,2), 0.5,self.nruns[i] )
            t0=self.steady[i]
            self.avalanches = a.s[t0:]
            unique = list(set(self.avalanches))
            unique.sort()
            total = float(len(self.avalanches))
            self.probs = []
            for value in unique:
                self.probs.append(float(self.avalanches.count(value))/total)
                
            print(np.sum(self.probs))
            

            # For 3b just comment out the line below to omit the normal Probs and set different runs parameters            
            plt.loglog(unique, self.probs, 'o', label=len(self.avalanches))
            print("Logbin start")
            d, c = log_bin(self.avalanches, 0., 1., 1.5, 'integer', drop_zeros = False, debug_mode=True)
            plt.loglog(d, c, label=(len(self.avalanches), "Log binned",1.5))
            #d, c = log_bin(self.avalanches, 1., 1.5, 1.75, 'integer', debug_mode=True)
            #plt.loglog(d, c, label=(len(self.avalanches), "Log binned",1.75))
            
            #d, c = log_bin(self.avalanches, 1., 1.5, 1.85, 'integer', debug_mode=True)
            #plt.loglog(d, c, label=(len(self.avalanches), "Log binned",1.85))
            #d, c = log_bin(self.avalanches, 1., 1.5, 2.0, 'integer', debug_mode=True)
            #plt.loglog(d, c, label=(len(self.avalanches), "Log binned",2.0))
            
             
            i+=1
        
        plt.legend()
        plt.show() 
        
        
    def call3c(self):
        """
        Redo for 3c, so far copied 3a/b        
        
        This method is for task 3b. Main difference from 3a is that it only plots the log binned data.
        When setting up the parameters make sure that all have the same steady state number of runs N that is large enough
        """
        plt.figure(7)
        plt.title("Avalanche size probability")
        plt.xlabel("Avalanche size")
        plt.ylabel("Probability")
        
        i=0
        while(i<self.difruns):

            a = Oslo(self.L[i],(1,2), 0.5,self.nruns[i] )
            t0=self.steady[i]
            self.avalanches = a.s[t0:]
            unique = list(set(self.avalanches))
            unique.sort()
            total = float(len(self.avalanches))
            self.probs = []
            for value in unique:
                self.probs.append(float(self.avalanches.count(value))/total)
                
            print(np.sum(self.probs))
                        
            #plt.loglog(unique, self.probs, 'ro', label=len(self.avalanches))
            centres, counts = log_bin(self.avalanches, 0., 1., 1.5, 'integer', drop_zeros = False, debug_mode=True)
            plt.loglog(centres, counts, 'r-', label=(len(self.avalanches)))
             
            i+=1
        
        plt.legend()
        plt.show() 
        return None
        

#tresholds in time are[54,226,898,3391,14056,56437,225745...]    
#b=Results(L=[128.], nruns=[25000], steady=[15000])   
#b=Results(L=[256., 256., 256], nruns=[70000,160000, 1060000], steady=[60000,60000,60000]) # 3a 
#b=Results(L=[16], nruns=[5300], steady=[300])  
