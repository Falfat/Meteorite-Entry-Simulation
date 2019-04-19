# Armageddon - Modelling atmospheric entry and disruption of asteroids

## Purpose

The dynamics of an asteroid in Earth’s atmosphere prior to break-up is governed by a coupled set of ordinary differential equations.
We have implemented a model that solves these equations for a set of given input parameters.

The key results that can be drawn from out model are:

-  Whether the scenario is an airburst or a cratering event
-  The peak kinetic energy loss per km (airburst)
-  The burst altitude (airburst)
-  The total kinetic energy loss at burst (airburst)
-  The mass and speed of the remnants of the asteroid that strike(s) the ground (airburst and cratering)

#### Functionality:

There is an extended range of functionality which enable the user to deploy our software for a wide range of cases.

-  Ability to apply the model to an arbitrary planetary atmosphere and to use a tabulated atmospheric density profile

-  Alternative or more complex treatments of fragmentation akin to (Passey and Melosh 1980; Chyba et al. 1993; Avramenko et al. 2014; Wheeler et al. 2017).

## Software Usage:
-  As an **installable python package**: <br>

    -  To install the package to your python files run the following command in a terminal in the directory containing "setup.py": <br> **`$ python setup.py install`**  
    -  This will allow you to utilise the python modules contained within /MeteorPackage without being located in our file directory.  <br> <br>
    
-  As a **command line programme**:

`$ python start.py v m theta r [-p DENSITY] [-z ALTITUDE] [-g GRAVITY] [-o RHO] [-s SIGMA0] [-t PA_TAB] [-l PLANET] 
[-x PLOT] [-y RESULTS]` <br> <br>
    -  Specifying the -x flag gives the user plots for the key model results <br>
    -  Specyfying the -y flag saves the model data to a file "Data.csv" with column headers: time, v, m, theta, r and a text readout specifying exact numbers for key events and values.
   
   # The  meteor class: 
The meteor class, within "meteor.py", contains the core functionalities of our software <br> 
input arguments as follows:
<br>
-  _Required parameters: ( **v** = velocity, **m** = mass, **theta** = meteorite entry angle, **r** = meteorite radius)_ <br>
-  Optional parameters: (z0 = initial altitude, g = gravitional acceleration, Rp = planetry radius, sigma0 =  material strength of asteroid, pa_tab=False, planet='Earth' )

#### pa_tab parameter:
-  Allows for user input of tabulated atmospheric pressure data files **or** an atmospheric input function <br>
-  May be passed a function or a file location of .csv tabulated data as a string in US Standard Atmosphere table format
-  If not specified an exponential function is assumed, as default representing Earth's atmosphere
-  If specified a cubic spline is fitted to the tabulated data and used to interpolate a density function, when passed a z value.

#### Planet parameter:
-  Provides a quick way to specify gravitational accelaration, planetry radius and atmospheric pressure function for Earth or Mars. 
-  If not specified Earth is set as defalt.
-  If specific g, rp or pa_tab(atmospheric tabular data) are specified, these override the specific Earth or Mars cases.

## Class methods:
**meteor.fire():**
-  This method creates the data from running the solver. 
-  The method returns: time, yfinal, ablated, t_events
    -  time = an array of the time values produced by the RK45 solver. Note the values are not of even timestep
    -  yfinal = a list of arrays with v, m, theta, z, x, r data respectively
    -  ablated = Bool denoting whether fracturing occurs
    -  t_events = list of arrays with values denoting time when the fracturing condition is met and impact respectively

**meteor.plot():**
-  This method created two seperate pop up plots first containing plots of: altitude vs time; speed vs altitude; mass vs. altitude; diameter vs. altitude. The second kinetic energy loss per km vs. altitude.
-  If called with .plot(save=True) the method will save the plots as .png files to your directory <br>

**meteor.results():**
-  Provides readout of key results of the model: 
    -  Whether the scenario is an airburst or a cratering event
    -  The peak kinetic energy loss per km (airburst)
    -  The burst altitude (airburst)
    -  The total kinetic energy loss at burst (airburst)
    -  The mass and speed of the remnants of the asteroid that strike(s) the ground (airburst and cratering) <br>

**meteor.analytical():**
-  Provides functionality to return an symplistic analytical solution
- returns a list of arrays containing analytical solutions for variables; velocity, trajectory angle, altitude and downrage distance of the asteroid.
    
# The meteorFCM Class:
The meteorFCM class, within "meteor_FCM.py" is an alternate class which uses the fragment cloud model to control fragmentation and radius expansion. The class functions similarly to the basic meteor class. The class takes <br>
-  _Required parameters: ( **v** = velocity, **m** = mass, **theta** = meteorite entry angle, **r** = meteorite radius)_ <br>
-  Optional parameters: (z0 = initial altitude, g = gravitional acceleration, Rp = planetry radius, sigma0 =  material strength of asteroid, pa_tab=False, planet='Earth' ) 

# The ideal_meteor Class:
The ideal_meteor class, within "ideal_meteor.py" is child class of meteor class, which used to apply the model in an ideal condition and conduct the comparison between the numerical solution and analytical solution. The class inherit all the properties of meteor class. So the inputs are the same with meteor class.

# The statistic module:
The statistic module, within "sta.py" is an module which to perform an ensemble of simulations with probability distributions as inputs and calculate the uncertainty in altitude and magnitude of peak energy release. The function in the module takes <br>
-  _Required parameters: ( **v** = mean value of velocity distribution, **theta** = mean value of theta distribution, **r** = mean value of radius distribution, **sigma_v** = sigma of velocity distribution, **sigma_theta** = sigma of theta distribution, **sigma_r** = sigma of radius distribution)_ <br>
-  Optional parameters: (PLOT=True ) 
####

