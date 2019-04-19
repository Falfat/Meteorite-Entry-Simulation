import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt
import scipy.interpolate as si
import Read_AltitudeDensityTable as rd


class meteorFCM:
    
    """
    A class which represents state of debris of an asteroid

    ...

    Attributes
    ----------

    pa_tab: bool
        "True" if pressure table is to be used, "False" if otherwise.
    planet: string
        "Earth" if we model entry on Earth, "Mars" if we model entry on Mars
    v: float
        asteroid velocity
    m: float
        asteroid mass
    theta: float
        meteoroid trajectory angle to the surface
    r: float
        asteroid radius
    z0: int
        altitude at which the asteroid enters the earth
    g0: float
        the gravitational acceleration at sea level.
    Rp: int
        planetary radius
    Cd: int
        drag coefficient
    sigma: float
        ablation parameter
    Cdisp: float
        dimensionless dispersion  coefficient
    P: float
        stagnation pressure
    pm: int
        asteroid density
    A: float
        cross sectional area of asteroid
    pressure: float
        earth's air ram pressure

    ALL TERMS ARE IN MKS UNITS

    Methods
    -------

    pressurecalc(z):
        Returns ram pressure values depending on if z is a single value
        or list.
    gcalc(z):
        Returns g-values depending on if z is a single value
        or list.
    func(t, y):
        A differentiable function pre ablation for scipy solver.
    funcablated(t, y):
        A differentiable fucntion post ablation for scipy solver.
        Here, the radius of asteroid is increasing due to ablation.
    fire():
        Generates the solutions to the ODEs in an array.
        solves for both airbust and crater events.
    ablation(t, y):
        Determines when ablation occurs. solve_ivp will minimise
        this to return t at ablation.
    carter(t, y):
        A function to determine if crater event occured or not.
    plot():
        Plots graphs of asteroid's velocity, mass, trajectory, altitude
        and downward range versus time.
    energypatch(plot=False):
        Prints a list of arrays of energy loss(kilotonnes) per km, asteroid
        altitude(km), time(seconds). A plot of energy loss per km
        is returned if plot value is "True".
    energyrelease():
        Determines energies and losses in case of an airbust event.
        Prints peak energy loss, maximum energyloss per km and a list of
        energy losses per km.
    results(readout=True):
        Provides readout of summary for an asteroid's fate upon entering
        the earth. Prints out of the time of aurburst/crater event, debris'/
        cater's touchdown time, weight and speed. Also, burst
        altitude.

        Example
        -------
        Creating an instance of the class and printing key results readout
        ::

        mathilde = meteorFCM(v=19.5e3, m=mass, theta=18.3/180 * np.pi, r=8.5, z0=100E3, g=9.81, Rp=6400E3, sigma0=4e6,
                          pa_tab=False, planet="Earth")

        mathilde.results()
    """

    def __init__(self, v, m, theta, r, z0=100E3, g0=9.81, Rp=6400E3, P=64*1E6, pa_tab = False, planet = 'Earth'):
        
        """
    
        Parameters
        ----------
        pa_tab: bool
           "True" if pressure table is to be used, "False" if otherwise.
        planet: string
           "Earth" if we model entry on Earth, "Mars" if we model entry on Mars
        v: float
            asteroid velocity
        m: float
            asteroid mass
        theta: float
            meteoroid trajectory angle to the surface
        r: float
            asteroid radius
        z0: int
            altitude at which the asteroid enters the earth
        g0: float
            the gravitational acceleration at sea level.
        Rp: int
            planetary radius
        Cd: int
            drag coefficient
        sigma: float
            ablation parameter
        Cdisp: float
            dimensionless dispersion  coefficient
        P: float
            stagnation pressure
        pm: int
            asteroid density
        A: float
            cross sectional area of asteroid
        pressure: float
            earth's air ram pressure

        ALL TERMS ARE IN MKS UNITS

        """
        
        if planet == 'Mars':
            g0=3.71
            Rp=3390E3
            
        self.pa_tab = pa_tab
        self.planet = planet
        self.v = v
        self.m = m
        self.theta = theta
        self.r = r
        self.z0 = z0
        self.g0 = g0
        self.Rp = Rp
        self.Cd = 1
        self.sigma = 10**(-8)
        self.Cdisp= 3.5
        self.P= P
        self.pm = 3300
        self.A = np.pi * self.r**2
        
        if self.pa_tab != False:
            if isinstance(self.pa_tab, str):
                z_list, pa_list, H_list = rd.read(pa_tab)
                zarray=np.array(z_list)
                parray=np.array(pa_list)
                Harray=np.array(H_list)
                self.pressure = si.CubicSpline(zarray, parray)
            elif hasattr(self.pa_tab, '__call__') is True: #checks if callable
                self.pressure = self.pa_tab

    def pressurecalc(self, z):
        """
        Function for pressure values using the FC Model.
        
        Parameters
        ----------
        z(int): a single altitude value if pa_tab is "False"
                or an array_like if pa_tab is "True"
        
        Returns
        -------
        a single ram pressure value if pa_tab is "False"
            or a list of ram pressures id pa_tab is "True"
            using the FC Model.
        """
        if self.pa_tab is False:
            #calculate the atmosphere density on Mars
            if self.planet is 'Mars':
                if z>7000:
                    T = -23.4 - 0.00222 * z
                else:
                    T = -31 - 0.000998 * z
                p = 0.699 * np.exp(-0.00009 * z)
                return p/(0.1921*(T+273.1))
            #calculate the atmosphere density on Earth(default)
            else:
                return -140.2*np.exp(-0.000187*z)+141.4*np.exp(-0.000186*z)
        else:
            return self.pressure(z)
        
    def gcalc(self, z):
        """
        Function to calculate accelaration due to gravity values
            at different altitudes for Fragment-Cloud model.
            
        Paramters
        ---------
        z(int)： values for altitudes
        
        Returns
        -------
        g-values at different altitudes for the Fragment-Cloud model.
        """
        return self.g0*((self.Rp/(self.Rp+z))**2)
   

    
    def func(self, t, y):
        """
        A differentiable function pre ablation for scipy solver.

        Parameters
        ----------
        t(float): an array_like of time-steps

        y(float): a list of initial conditions for differentiable variables

        Returns
        -------
        An array containing the values of differentiable variables for each
            desired time in t, with the initial values in the first row.

        """
        v = y[0]
        m = y[1]
        theta = y[2]
        z = y[3]
        pa = self.pressurecalc(z)
        g = self.gcalc(z)
        dv = ((-self.Cd * pa * self.A * v**2) / (2*m)) + (g * np.sin(theta))
        dm = (-self.sigma * pa * self.A * v**3)/ 2 
        dtheta = (v/(self.Rp+z)+g/v)*np.cos(theta)
        dz = -v * np.sin(theta)
        return [dv, dm, dtheta, dz]  

    def funcablated(self, t, y):
        """
        differentiable fucntion post ablation for scipy solver.
        
        Parameters
        ----------
        t(float): an array_like of time-steps

        y(float): a list of conditions for differentiable variables

        Returns
        -------
        An array containing the values of differentiable variables for each
            desired time in t, with the initial values in the first row.

        """
        v = y[0]
        m = y[1]
        theta = y[2]
        z = y[3]
        r = y[4]
        pa = self.pressurecalc(z)
        g = self.gcalc(z)

        dv = ((-self.Cd * pa * self.A * v**2) / (2*m)) + (g * np.sin(theta))
        dm = (-self.sigma * pa * self.A * v**3)/ 2 
        dtheta = (v/(self.Rp+z)+g/v)*np.cos(theta)
        dz = -v * np.sin(theta)
        # change z0 to height at y1
        if (self.pressurecalc(y[3]) * y[0]**2) > self.P:
            dr = (np.sqrt((self.Cdisp * (pa/self.pm))))*v
        else:
            dr=0
        return [dv, dm, dtheta, dz, dr]  
    
    def fire(self):  
        """
        Generates the solutions to the ODEs in an array.
            solves for both airbust and crater events.

        Parameters
        ----------
        None

        Returns
        -------
        The solution time grid, numerical solutions for all differentiable
            variables, a list of arrays at which event/s (ablation/crater)
            was detected and a confirmation if ablation occurs or not.

        """
        y0 = [self.v, self.m, self.theta, self.z0]
        soln = sp.solve_ivp(self.func, [0, 20000], y0, dense_output=True, max_step=0.1, events=[self.ablation, self.crater])
        if len(soln.t_events[0]) != 0:  # ablation has occured
            y1 = np.zeros(len(y0) + 1)
            y1[:-1] = soln.y[:, -1]  # new initial conditions are output from old
            y1[-1] = self.r
            ablated = True
            soln1 = sp.solve_ivp(self.funcablated, [soln.t[-1], 20000], y1, dense_output=True, max_step=0.1, events=self.crater)
            # give from ablation to impact
            # need to stitch together soln and soln1
            yfinal = [np.append(soln.y[i], soln1.y[i]) for i in range(len(soln.y))]
            yfinal.append(np.append(np.full_like(soln.y[0], self.r), soln1.y[-1]))
            time = np.append(soln.t, soln1.t)
            t_events = np.append(soln.t_events[0],soln1.t_events[0]) # ablation time and crater time
        else:
            ablated = False
            yfinal = [soln.y[i] for i in range(len(soln.y))]
            yfinal.append(np.full_like(soln.y[0], self.r))
            time = soln.t
            t_events = soln.t_events[1] # crater event
        return time, yfinal, ablated, t_events

    

    def ablation(self, t, y):
        """
        A function to determine when ablation occurs.
            solve_ivp will minimise this func to return t at ablation.

        Parameters
        ----------
        t(float): an array_like of time-steps

        y(float): a list of initial conditions for differentiable variables

        Returns
        -------
        zero when ablation occurs and absolute value of the difference
            between ram pressure and asteroid strength prior to ablation.

        """
        if (self.pressurecalc(y[3]) * y[0]**2) > self.P:
            return 0
        else:
            return abs((self.pressurecalc(y[3]) * y[0]**2) - self.P)
   
    def crater(self, t, y):
        """
        A function to determine if crater event occured or not.

        Parameters
        ----------
        t(float): an array_like of time-steps

        y(float): a list of initial conditions for differentiable variables

        Returns
        -------
        zero once crater event occurs, else returns asteroid's altitiude.

        """
        if y[3] < 0:
            return 0
        else:
            return y[3]  # will minimise when z goes to 0
    
    def plot(self, save="False"):
        """
        Function to plot graphs of desired variables versus time.

        Parameters
        ----------
        save(bool) : saves plots to .png (default=False)

        Returns
        -------
        Plots of asteroid's velocity, mass, trajectory, altitude,
            downward range against time.
        """        
        # create plots for z/t, v/z, m/z, r/z and Energy Release/z
        data = self.fire()
        energydata = self.energypatch()

        plt.subplot(2, 1, 1)
        plt.plot(data[0], data[1][3], 'r-')
        plt.title('Meterorite Model')
        plt.ylabel('Altitude (m)')
        plt.xlabel('Time (s)')
        plt.xticks(size='small')
        plt.yticks(size='small')

        plt.subplot(2, 3, 4)
        plt.plot(data[1][3], data[1][0], 'g-')
        plt.xlabel('Altitude')
        plt.ylabel('Velocity (m/s)')
        plt.xticks(size='small')
        plt.yticks(size='small')
        plt.gca().invert_xaxis()

        plt.subplot(2, 3, 5)
        plt.plot(data[1][3], data[1][1], 'g-')
        plt.xlabel('Altitude')
        plt.ylabel('Mass (kg)')
        plt.xticks(size='small')
        plt.yticks(size='small')
        plt.gca().invert_xaxis()

        plt.subplot(2, 3, 6)
        plt.plot(data[1][3], data[1][4], 'g-')
        plt.xlabel('Altitude')
        plt.ylabel('Radius (metres)')
        plt.xticks(size='small')
        plt.yticks(size='small')
        plt.gca().invert_xaxis()
        plt.tight_layout()

        if save == True:
            plt.savefig('KeyResultsPlot.png')

        plt.show()
        plt.title('Meterorite Energy Release')
        plt.subplot(1, 1, 1)
        plt.plot(energydata[1], energydata[2])
        plt.xlabel('Energy Release Kt/Km')
        plt.ylabel('Altitude')
        plt.tight_layout()
        plt.show()

        if save == True:
            plt.savefig('EnergyRelease.png')

    def energypatch(self, plot=False):
        """
        Function to plot graphs of desired variables versus time.

        Parameters
        ----------
        None

        Returns
        -------
        Plots of asteroid's velocity, mass, trajectory, altitude,
            downward range against time.
        """     
        time, yfinal, ablated, t_events = self.fire()
        v = yfinal[0]
        m = yfinal[1]
        z = yfinal[3]
        energy = 0.5 * m * v**2
        energydiff = np.diff(energy)
        zdiff = np.diff(z)
        singular = np.where(zdiff==0)[0][0]
        zdiff[singular] = (zdiff[singular-1]+zdiff[singular+1])
        energydiff = energydiff / zdiff
        energydiff[singular] = (energydiff[singular+1] + energydiff[singular-1])/2
        energydiff = np.append(energydiff, 0)
        energydiff = energydiff / (4.184*1E9)  # convert kt per kilometer
        z = z/1E3
        if plot is True:
            plt.figure()
            plt.plot(energydiff, z)
        return time, energydiff, z

    def energyrelease(self):
        """
        Determines energies and losses in case of an airbust event.

        Parameters
        ----------
        None

        Returns
        -------
        The peak energy loss, maximum energyloss per km and a list of
            energy losses per km.
        """
        time, yfinal, ablated, t_events = self.fire()
        v = yfinal[0]
        m = yfinal[1]
        z = yfinal[3]
        energy = 0.5 * m * v**2
        energydiff = np.diff(energy)
        zdiff = np.diff(z)
        singular = np.where(zdiff == 0)[0][0]
        zdiff[singular] = (zdiff[singular-1]+zdiff[singular+1])
        energydiff = energydiff / zdiff
        energydiff[singular] = (energydiff[singular+1] + energydiff[singular-1])/2
        
        energydiff = np.append(energydiff, 0)
        energydiff = energydiff * 1E3 / (4.184*1E12)  # convert kt per kilometer
        z = z/1E3
        maxenergydiff = energydiff.max()  # RETURNED IN MTONNES AND KM
        peakenergyloss = z[energydiff.argmax()]

        return energydiff, maxenergydiff, peakenergyloss


    def results(self):
        """
        Provides readout of summary for an asteroid's fate upon entering
            the earth.

        Parameters
        ----------
        None

        Returns
        -------
        Print out of the time of aurburst/crater event, debris'/
            cater's touchdown time, weight and speed. Also, burst
            altitude.

        """
        time, yfinal, ablated, t_events = self.fire()

        if ablated is True:
            print("Airbust event occured at %6.2f seconds" % (t_events[0]))
            print("Debris field touchdown at %6.2f seconds" % (t_events[1]))
            print("Debris field was %6.2f tons and was travelling at %6.2f" % (yfinal[1][-1]/1000, yfinal[0][-1]))

        else:
            print("Cratering occured without airburst, debris touchdown at %7.2f seconds" % (t_events[0]))
            print("Meteoroid was of mass %6.2f tons and travelling at %6.2f m/s at touchdown" % (yfinal[1][-1]/1000, yfinal[0][-1]))

        ediff, maxenergyloss, elossaltitude = self.energyrelease()

        if ablated is True:
            print("maximum energy loss was %6.2f kt/km at %6.2f KM" % (maxenergyloss, elossaltitude))


###############################
meteorFCM.crater.terminal = True  # MONKEY PATCH
meteorFCM.ablation.terminal = True
