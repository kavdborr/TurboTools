# -*- coding: utf-8 -*-
"""
Created on Sun May  1 14:30:52 2022

@author: Karel Van den Borre
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from xml.dom import minidom
from scipy.interpolate import interp2d
import os

class TurboMap:
    """
    This class contains tools to load, convert and display compressor maps.
    
    Parameters
    ----------
    fileLoc:str
        Location of teh map that needs to be read
        
    width:int
        width of the map
        
    height:int
        height of the map
        
    isCompressor:Bool
        is the map a compressor?
        
    name:string
        opt: provides a custom name for the map. Otherwise the filename is used.
    
    """
    
    name = "" # name of the map for plotting and savinf purposes
    
    isTurbine = False   # is this the map of a turbine?
    isCompressor = False    # is this the map of a compressor?
    
    massFlow = np.zeros([10, 10])   # The map of massflow
    compression = np.zeros([10, 10])    # map of the compression ratio
    efficiency = np.zeros([10, 10]) # map of the adiabatic efficiency
    
    Kp = 1  # compression ratio scaling factor
    Km = 1  # massflow scaling factor
    Keta = 1    # efficiency scaling factor
    
    beta = np.ones([10,1])  # list of beta describing the map
    N = np.ones([10,1]) # list of non dimensional rotational speeds
    
    width = 10  # width of the map (number of columns in ine field)
    height = 10 # height of the map (number of rows in one field)
    
    surgeM = np.array([0,1])
    surgeP = np.array([0,1])
    
    def __init__(self, fileLoc, width, height, isCompressor, name="default"):
        """
        initialiser for the object. Populates all fields needed to display a (non-scaled) map
        the width and height are user provided as the file structure is not always consistent
        
        """
        
        # has a custom name been specified? if so use it!
        if not(name == "default"):
            self.name=name
        else:
            # Take the name of the file without extension as default name
            self.name = fileLoc.split(".")[-2].split("/")[-1]
        
        # set width and height
        self.width=width
        self.height=height
        
        # set wheter compressor or turbine
        if isCompressor:
            self.isCompressor=True
        else:
            self.isTurbine=True
        
        # read the map file in the gasturb format if the extention is right
        if ".Map" in fileLoc:
            self.readMapGT(fileLoc, width, height)
           

    def scaleMap(self, Kp, Km, Keta):
        """function to easily scale the whole map in one line
        
        Parameters
        ----------
        Kp:float
            scaling factor for the pressure ratio
            
        Km:float
            scaling factor for the massflowrate
            
        Keta:float
            scaling factor for the adiabatic efficiency
        """
        self.Kp=Kp
        self.Km=Km
        self.Keta=Keta

    
    def readMapGT(self, fileLoc, width, height):
        """
        Reads the map in the Gasturb format. The Gasturb format can vary slighlty in setup, therefore it is easier to make the user specify the size of the map.
        
        Parameters
        ----------
        fileLoc:str
            Location of teh map that needs to be read
            
        width:int
            width of the map
            
        height:int
            height of the map
            
        """
        mPR = np.zeros([0,1])   # minimum pressure ratio for turbine construction
        MPR = np.zeros([0,1])   # maximum pressure ratio for turbine construction
        fid = open(fileLoc)
        while True:
            tline = fid.readline()
            
            #exit condition is end of file
            if not(tline):
                break
            
            # read the mass flow map
            if "Mass Flow" in tline:
                self.beta, self.N, self.massFlow = self.__readSection(width, height, fid)
                
            # read the efficiency map
            if "Efficiency" in tline:
                self.beta, self.N, self.efficiency = self.__readSection(width, height, fid)  
                
            # read the pressure ratio map (only for compressors)
            if ("Pressure Ratio" in tline) and not("M" in tline):
                self.beta, self.N, self.compression = self.__readSection(width, height, fid)
                
            # read the min and max pressure ratios of the pressure map (only for turbines)
            if "Min Pressure Ratio" in tline:
                mPR = self.__readTurbineSection(fid, width, height)
            if "Max Pressure Ratio" in tline:
                MPR = self.__readTurbineSection(fid, width, height)
                
        # construct the compression ratio map in case of a turbine
        if np.size(mPR,0)!=0:
            self.compression=np.zeros([height, width])
            for i in range(np.size(mPR)):
                self.compression[i,:]=np.linspace(mPR[i], MPR[i], width)
    
    
    def writeMap(self, fileLoc):
        """
        Write the map in the TURBO format, using the xml format.
        For the goal of the current project only gasturb to TURBO conversion is needed.
        In future versions support for the inverse may also become available.
        
        Parameters
        ----------
        fileLoc:str
            fileName for the map, which is saved in convertedMaps.
            
        """
        
        # construct header
        doc = minidom.Document()
        l1 = doc.createElement("table")
        l1.setAttribute("type", "3D")
        l1.setAttribute("name", "Characteristic map " + self.name)
        l1.setAttribute("description", "DESCRIPTION")
        l1.setAttribute("version", "1.0")
        l1.setAttribute("cdate", date.today().strftime("%d/%m/%Y"))
        l1.setAttribute("mdate", date.today().strftime("%d/%m/%Y"))
        l1.setAttribute("revision", "1.0")

        doc.appendChild(l1)

        l2 = doc.createElement("interp")
        l2.setAttribute("default", "LINEAR")
        l2.setAttribute("valid", "{FORBIDDEN, LINEAR}")
        l1.appendChild(l2)
        l3 = doc.createElement("extrap")
        l3.setAttribute("default", "LINEAR")
        l3.setAttribute("valid", "{FORBIDDEN, LINEAR, SPLINE}")
        l1.appendChild(l3)
        
        l4 = doc.createElement("axis1")
        l4.setAttribute("id", "Index")
        l4.setAttribute("id", "Field indicator (-)")
        l4.setAttribute("value", "{ 1,2,3 }")
        l1.appendChild(l4)
        
        # create field for corrected speed
        l5 = doc.createElement("axis2")
        l5.setAttribute("id", "index")
        l5.setAttribute("description", "Corrected speed (-)")
        l5.setAttribute("value", np.array2string(self.N, separator=",", max_line_width=99999, precision=5, floatmode="fixed").replace("[","{").replace("]","}"))
        l1.appendChild(l5)
        
        # create field for beta
        l6 = doc.createElement("axis3")
        l6.setAttribute("id", "Beta")
        l6.setAttribute("description", "Beta Parameter (-)")
        l6.setAttribute("value", np.array2string(self.beta, separator=",", max_line_width=99999, precision=5, floatmode="fixed").replace("[","{").replace("]","}"))
        l1.appendChild(l6)
        
        # create field for massflow, efficiency and pressure ratio
        l7 = doc.createElement("return")
        l7.setAttribute("id", "Value")
        l7.setAttribute("description", "Field value: [Index=1] Corrected mass flow rate (kg/s); [Index=2] Adiabatic efficiency (-); [Index=3] Pressure ratio t-t (-)")
        # compose string massflowrate
        sM = np.array2string(self.massFlow, separator=",", max_line_width=99999, precision=5, floatmode="fixed").replace("[","{").replace("]","}")
        # compose string efficieny
        sE = np.array2string(self.efficiency, separator=",", max_line_width=99999, precision=5, floatmode="fixed").replace("[","{").replace("]","}")
        # compose string pressure ratio
        sC = np.array2string(self.compression, separator=",", max_line_width=99999, precision=5, floatmode="fixed").replace("[","{").replace("]","}")

        l7.setAttribute("value", "{\n" + sM + ",\n" + sE + ",\n" + sC + "\n}")

        l1.appendChild(l7)

        #create folder to save converted map if the folder dos not exist
        if not os.path.exists("convertedMaps"):
            os.mkdir("convertedMaps")
            
        # save files 
        xml_str = doc.toprettyxml(indent="\t")
        with open("convertedMaps/" + fileLoc, "w") as f:
            f.write(xml_str)
    
    
    def printMap(self):
        """
        Display the scaled map. formatting is autamically performed.
        After calling this function more plotting actions can be performed on the same figure.
        
        Parameters
        ----------
        None
        """
        #always start a new figure when plotting a map
        plt.figure()

        #for compressors, plot the efficiency in terms of massflow and pressure ratio
        if self.isCompressor:
            plt.contourf(self.massFlow*self.Km, (self.compression-1)*self.Kp+1, self.efficiency*self.Keta, 20, cmap="gray")
            # plot the isospeed lines
            for i in range(np.size(self.N)):
                plt.plot(self.massFlow[i,:]*self.Km, (self.compression[i,:]-1)*self.Kp+1, "#1c1c1c", linewidth=1)
        # for turbine, plot efficiency in terms of massflow times rotational speed and pressure ratio
        if self.isTurbine:
            plt.contourf((self.massFlow.T*(self.N)).T*self.Km, (self.compression-1)*self.Kp+1, self.efficiency*self.Keta, 20, cmap="gray")
            # plot the isospeed lines
            for i in range(np.size(self.N)):
                plt.plot(self.massFlow[i,:]*self.N[i]*self.Km, (self.compression[i,:]-1)*self.Kp+1, "#1c1c1c", linewidth=1)

        # formatting of the plot()
        plt.grid(True)
        plt.title(self.name)
        plt.ylabel(r"$\Pi$ [-]")
        if self.isCompressor:
            plt.xlabel(r"$\tilde{\dot{m}}$ [kg/s]")
        else:
            plt.xlabel(r"$\tilde{\dot{m}} \cdot{} \tilde{N}$ [kg/s]")
        plt.colorbar(label=r"efficiency [-]")
    
    
    def plotDesingPoint(self, beta, N):
        """
        Plots the design point (in terms of N and beta) on the map.
        
        Parameters
        ----------
        beta:float
            Design beta value of the machine.
        
        N:float
            Design corrected rotaional speed of the machine.
            
        """
        
        # covert design point from beta and N to massflow and compression ratio
        mf = interp2d(self.N, self.beta, self.massFlow.T)
        pr = interp2d(self.N, self.beta, self.compression.T)
        # plot the actual point
        if self.isCompressor:
            plt.plot(mf(N,beta)*self.Km, pr(N,beta)*self.Kp, 'ro', markersize=10)
        if self.isTurbine:
            plt.plot(mf(N,beta)*self.Km*N, pr(N,beta)*self.Kp, 'ro', markersize=10)
    
    
    def plotExperiment(self, test, machine, Color="rainbow", res=50):
        """
        This function plots the movement of the operating point on the map.
        The color of the line can vary with time if specified color is rainbow.
        
        Parameters
        ----------
        test:DataFrame
            Pandas dataframe containing the results of a simulation performed in EcosimPro.
        
        machine:str
            name of the machine, as reported in the EcosimPro report file.
            
        Color:str
            Color of the operating point locus. 
            If not specified, use rainbow, where the line varies 
            from red ot green in finction of time.
            
        res:int
            Color resolution of the line in rainbow mode. Default is 50. 
            Choosing unreasonably large values might resutl in slow plotting
        """
        # the first two are boring monochrome plots
        if Color!="rainbow" and self.isCompressor:
            plt.plot(test[machine+".Wcorr(Kg/s)"], test[machine+".pi(-)"], color=Color)
        elif Color!="rainbow" and self.isTurbine:
            plt.plot(test[machine+".Wcorr(Kg/s)"]*test[machine+".Ncorr(-)"], test[machine+".pi(-)"], color=Color)
        # these two are gradient plots, varying the linecolor from red to green with time
        else:
            # min and max time to calculate color gradient
            Tmin = test["TIME"].iloc[0]
            Tmax = test["TIME"].iloc[-1]
            L = test.shape[0]
            #copy required data from pandas to numpy for efficiency
            time = test["TIME"].to_numpy()
            wcorr = test[machine+".Wcorr(Kg/s)"].to_numpy()
            pi = test[machine+".pi(-)"].to_numpy()
            # if turbine, scale massflow woth rotational speed
            if self.isTurbine:
                wcorr = wcorr*test[machine+".Ncorr(-)"].to_numpy()
            for i in range(res):
                # calculate the percentage along the experiment timeline
                p = (time[int(i/res*L)]-Tmin)/(Tmax-Tmin)
                # set the color in rgb values
                Color = [1-p, p, 0.0]
                # plot a line segment
                plt.plot(wcorr[int((i/res)*L):int((i+1)/res*L)], pi[int((i)/res*L):int((i+1)/res*L)], color=Color)
        
    
    
    def __readSection(self, width, height, fid):
        """
        Function to read a block of the Gasturb text file.
        Helper function of readMapGT only, as such it can be left private
        
        Parameters
        ----------
        width:int
            width of the map
            
        height:int
            height of the map
            
        fid:file
            file handler of the Gasturb map file, to be read line by line
            
        Returns
        -------
        beta:array
            1D array of beta values describing the map.
        
        N:array
            1D array of corrected rotational speed values describing the map.
        
        section:2D array
            2D array of the field which has been read (pressure, massflow or efficiency)
        
        """
        
        totalLines = math.ceil((width+1)/5)*(height+1)
        section = np.zeros([0])
        for i in range(totalLines):
            tline=fid.readline() # capture raw line
            tline = tline.replace(",", ".") # replace comma by dot as correct decimal
            tline = np.array(tline.split()) # divide line at whitespaces
            tline[tline==""]=np.array([]) # replace empty cells
            tline = tline.astype(np.float) # convert text to float
            section = np.concatenate((section, tline))
            
        # beta values are the first "column"
        beta = section[1:width+1]
        # delete the column with betas
        section=section[width+1:]
        
        # the first element of each "row" is N
        N = section[np.arange(0,width*height,width+1)]
        # delete the N elemets from the rows
        section = np.delete(section, np.arange(0,width*height,width+1))
        
        # reshape the vector to a matrix
        section = np.reshape(section, [int(np.size(section)/width), width])
        
        # returnb beta, N and the map
        return beta, N, section
    
    
    def __readTurbineSection(self, fid, width, height):
        """
        Reads minimum and maximum pressure ratio for turbine maps. 
        Helper function of readMapGT. Private as it is only needed within the class.
        
        Parameters
        ----------
        width:int
            width of the map
            
        height:int
            height of the map
            
        fid:file
            file handler of the Gasturb map file, to be read line by line
        
        
        Returns
        -------
        section:2D array
            return the section of the gasturb map, saved as 2D array
        """
        # the methods are in general completely analogous to the readsection class, except with different bounds
        totalLines = math.ceil((height+1)/5)*2
        section = np.zeros([0])
        for i in range(totalLines):
            tline = fid.readline()
            tline = tline.replace(",",".")
            tline = np.array(tline.split())
            tline[tline==""] = []
            tline = tline.astype(np.float)
            section = np.concatenate([section, tline])
        # drop the first column, as it has useless data
        section=section[height+2:]
        return section