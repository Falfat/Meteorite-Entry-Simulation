# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 13:12:45 2018
@author: Gareth
"""
import numpy as np
import scipy.integrate as sp
import matplotlib.pyplot as plt
import scipy.interpolate as si
import Read_AltitudeDensityTable as rd
import scipy.optimize as so


class meteor:
    """
    A class which represents state of debris of an asteroid

    ...

    Attributes
    ----------

    pa_tab: bool
        "True" if pressure table is to be used, "False" if otherwise.
    v: float
        asteroid velocity
    m: float
        asteroid mass
    theta: float
        meteoroid trajectory angle to the surface
    x: int
        downrange distance of the meteoroid from its entry position
    r: float
        asteroid radius
    z0: int
        altitude at which the asteroid enters the earth
    g: float
        acceleration due to gravity
    Rp: int
        planetary radius
    Cd: int
        drag coefficient
    Ch: float
        ablation efficiency coefficient
    Q: int
        specific heat of ablation
    Cl: float
        lift coefficient
    a: float
        dispersion coefficient
    H: int
       atmospheric scale height
    sigma0: int
        strength of the asteroid
    pm: int
        asteroid density
    A: float
        cross sectional area of asteroid
    p0: float
        atmospheric pressure at z0
    pressure: float
        earth's air ram pressure

    ALL TERMS ARE IN MKS UNITS

    Methods
    -------

    pressurecalc(z):
        Returns ram pressure value/s depending on if z is a single value
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
    analytical():
        Provides analytucal solutions to the ODEs, considering simplifying
        assumptions; pa = p0 * exp(-z(t)/H), g = 0, Rp = ∞, Cl = 0,
        dm/dt = 0, sigma0 = ∞, dr/dt = 0
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

        mathilde = meteor(v=19.5e3, m=mass, theta=18.3/180 * np.pi, r=8.5, z0=100E3, g=9.81, Rp=6400E3, sigma0=4e6,
                          pa_tab=False, planet="Earth")

        mathilde.results()




    """

    def __init__(self, v, m, theta, r, pm=3300, z0=100E3, g=9.81, Rp=6400E3, sigma0=64*1E6, pa_tab=False, planet='Earth'):
        """
        Parameters
        ----------
        pa_tab: bool
        "True" if pressure table is to be used, "False" if otherwise.
        v: float
            asteroid velocity
        m: float
            asteroid mass
        theta: float
            meteoroid trajectory angle to the surface
        r: float
            asteroid radius
        z0: int
            altitude at which asteroid enters the earth (default value: 100E3)
        g: float
            acceleration due to gravity (default value: 9.81)
        Rp: int
            planetary radius (default value: 6400E3)
        Cd: int
            drag coefficient (constant value: 1)
        Ch: float
            ablation efficiency coefficient (constant value: 0.1)
        Q: int
            specific heat of ablation (constant value: 10E7)
        Cl: float
            lift coefficient (constant value: 10E-3)
        a: float
            dispersion coefficient (constant value: 0.3)
        H: int
            atmospheric scale height (constant value: 8000)
        sigma0: int
            strength of the asteroid
        pm: int
            asteroid density (constant value: 3000)
        A: float
            cross sectional area of asteroid
        p0: float
            atmospheric pressure at z0 (constant value: 1.2)
        pressure: float
            earth's air ram pressure

        ALL TERMS ARE IN MKS UNITS

        """

        self.planet = planet

        if planet == 'Mars':
            g = 3.71
            Rp = 3390E3

        self.pa_tab = pa_tab
        self.v = v
        self.m = m
        self.theta = theta
        self.x = 0
        self.r = r
        self.z0 = z0
        self.g = g
        self.Rp = Rp
        self.Cd = 1
        self.Ch = 0.1
        self.Q = 1e7
        self.Cl = 10**-3
        self.a = 0.3
        self.H = 8000
        self.sigma0 = sigma0
        self.pm = pm
        self.A = np.pi * self.r**2
        self.p0 = 1.2
        # initialise z and p arrays for tabulated atmosphere data input
        if self.pa_tab != False:
            if isinstance(self.pa_tab, str):  # checks if string is supplied
                z_list, pa_list, H_list = rd.read(pa_tab)
                zarray = np.array(z_list)
                parray = np.array(pa_list)
                Harray = np.array(H_list)
                self.pressure = si.CubicSpline(zarray, parray)
            elif hasattr(self.pa_tab, '__call__') is True:  # checks if callable
                self.pressure = self.pa_tab

    def pressurecalc(self, z):
        """
        Returns ram pressure value/s depending on
            if z is a single value or list.

        Parameters
        ----------
        z(int): a single altitude value if pa_tab is "False"
            or an array_like if pa_tab is "True".

        Returns
        -------
        a single ram pressure value if pa_tab is "False"
            or a list of ram pressure values if pa_tab is "True".

        """
        if self.pa_tab is False:
            if self.planet == 'Mars':  # calculate pressures for mars conditions
                if z > 7000:
                    T = -23.4 - 0.00222 * z
                else:  # calculate the atmosphere density on Earth(default)
                    T = -31 - 0.000998 * z
                p = 0.699 * np.exp(-0.00009 * z)
                return p / (0.1921 * (T + 273.1))
            else:
                return self.p0 * np.exp(-z / self.H)

        else:
            return self.pressure(z)  # atmospheric tabulated data function

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
        v = y[0]  # initialises values for solver
        m = y[1]
        theta = y[2]
        z = y[3]
        x = y[4]
        r = y[5]
        A = np.pi * r**2
        pa = self.pressurecalc(z)  # calls externally specified pressure func
        dr = 0
        if pa * v**2 > self.sigma0:  # sets fracturing condition
            dr = ((3.5 * self.a * (pa/self.pm))**0.5)*v
        dv = ((-self.Cd * pa * A * v**2) / (2*m)) + (self.g * np.sin(theta))
        dm = (-self.Ch * pa * A * v**3) / (2 * self.Q)
        dtheta = ((self.g * np.cos(theta))/v) - ((self.Cl * pa * A * v)/(2 * m)) - (v * np.cos(theta)/(self.Rp + z))
        dz = -v * np.sin(theta)
        dx = (v * np.cos(theta))/(1+z/self.Rp)
        return [dv, dm, dtheta, dz, dx, dr]

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
        y0 = [self.v, self.m, self.theta, self.z0, self.x, self.r]  # initialises first solver

        soln = sp.solve_ivp(self.func, [0, 20000], y0, dense_output=True,  rtol=1e-8, events=[self.ablation, self.crater])
        # will run until the solver hits the floor.
        yfinal = soln.y
        time = soln.t
        t_events = soln.t_events
        if len(soln.t_events[0]) != 0:  # checks if event function has converged
            ablated = True
        else:
            ablated = False
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
        if (self.pressurecalc(y[3]) * y[0]**2) > self.sigma0:
            return 0  # 0 condition given for any convergence over to make location of event parameters easier to hit
        else:
            return abs((self.pressurecalc(y[3]) * y[0]**2) - self.sigma0)

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

    def analytical(self):
        """
        Provides analytical solutions to the ODEs, considering simplifying
            assumptions; pa = p0 * exp(-z(t)/H), g = 0, Rp = ∞, Cl = 0,
            dm/dt = 0, sigma0 = ∞, dr/dt = 0

        Parameters
        ---------
        None

        Returns
        -------
        A list of arrays containing analytical solutions for variables;
            velocity, trajectory angle, altitude and downrage distance
            of the asteroid.

        """

        p0 = self.p0
        v0 = self.v
        theta0 = self.theta
        z0 = self.z0
        Cd = self.Cd
        A = self.A
        Cl = self.Cl
        m = self.m
        H = self.H

        K = - (Cd * A * p0 * H)/(2 * m * np.sin(theta0))
        C = np.log(v0) - K * np.exp(-z0/H)
        print('work')
        datanum = self.fire()
        z1 = datanum[1][3]
        v = np.exp(K * np.exp(-z1/H) + C)
        v1 = datanum[1][0]
        z = -H*np.log((np.log(v1)-C)/K)
        return v, z, datanum

    def plot(self, save=False):
        """
        Function to plot graphs of desired variables versus time.

        Parameters
        ----------
        save(bool) :  saves plots to .png files (defualt=false)

        Returns
        -------
        Plots of asteroid's velocity, mass, trajectory, altitude,
            downward range against time.
        """

        data = self.fire()
        energydata = self.energypatch()

        plt.subplot(2, 1, 1)
        plt.plot(data[0], data[1][3], 'r-')
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
        plt.plot(data[1][3], data[1][5], 'g-')
        plt.xlabel('Altitude')
        plt.ylabel('Radius (metres)')
        plt.xticks(size='small')
        plt.yticks(size='small')
        plt.gca().invert_xaxis()

        plt.tight_layout()

        if save:
            plt.savefig('KeyResultsPlot.png')

        plt.show()

        plt.title('Meterorite Energy Release')
        plt.subplot(1, 1, 1)
        plt.plot(energydata[1], energydata[2])
        plt.xlabel('Energy Release Kt/Km')
        plt.ylabel('Altitude')
        plt.xlim(left = 0)

        plt.tight_layout()
        plt.show()

        if save:
            plt.savefig('EnergyRelease.png')

    def energypatch(self, plot=False):
        """
        Function for the asteroid energy loss to the air.

        Parameters
        ----------
        plot(bool): permits the return of plot for kinetic energy loss per km
                    if "True" (Default value: "False")

        Returns
        -------
        A list of arrays of energy loss(kilotonnes) per km, asteroid
            altitude(km), time(seconds). A plot of energy loss per km
            is returned if plot value is "True"

        """
        time, yfinal, ablated, t_events = self.fire()
        v = yfinal[0]
        m = yfinal[1]
        z = yfinal[3]
        energy = 0.5 * m * v**2

        energydiff = np.diff(energy)
        energydiff = np.append(energydiff, 0)  # pads out diff(array) to previous size
        energydiff = energydiff * -1  # corrects backwards difference give negative sign
        energydiff = energydiff / (4.184*1E12)  # convert kt per kilometer

        z = z/1E3  # convert to kilometer
        energydiff = -1 * energydiff / np.append(np.diff(z), 1)
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
        energydiff = np.append(energydiff, 0)
        energydiff = energydiff * -1
        energydiff = energydiff / (4.184*1E12)  # convert kt per kilometer

        z = z/1E3
        energydiff = -1 * energydiff / np.append(np.diff(z), 1)
        maxenergydiff = energydiff.max()  # returned in kt TNT / km
        peakenergyloss = z[energydiff.argmax()]

        return energydiff, maxenergydiff, peakenergyloss

    def results(self, readout=True):
        """
        Provides readout of summary for an asteroid's fate upon entering
            the earth.

        Parameters
        ----------
        readout(bool): ensures the print out of result summary if "True"
                       and otherwise if false. (Default value: True)

        Returns
        -------
        Print out of the time of aurburst/crater event, debris'/
            cater's touchdown time, weight and speed. Also, burst
            altitude.

        """
        time, yfinal, ablated, t_events = self.fire()
        ediff, maxenergyloss, elossaltitude = self.energyrelease()
        lossindex = np.where(yfinal[3] == elossaltitude*1e3)  # finding index of burst point
        keloss = (0.5 * yfinal[1][0] * yfinal[0][0]**2) - (0.5 * yfinal[1][lossindex] * yfinal[0][lossindex]**2)

        if ablated is True:
            print("Airbust event occured at %6.2f seconds" % (t_events[0][0]))
            print("Debris field touchdown at %6.2f seconds" % (t_events[1][0]))
            print("Debris field was %6.2f tons and was travelling at %6.2f meters per second" % (yfinal[1][-1]/1000, yfinal[0][-1]))
            print("Total kinetic energy loss at burst %6.2f joules" % (keloss))
            if yfinal[3][np.where(time == t_events[0][0])[0][0]] < 5: # checks z value of first fracture
                print("Fracturing occured below 5km producing cratering")

        else:
            print("Cratering occured without airburst, debris touchdown at %7.2f seconds" % (t_events[1][0]))
            print("Meteoroid was of mass %6.2f tons and travelling at %6.2f m/s at touchdown" % (yfinal[1][-1]/1000, yfinal[0][-1]))

        if ablated is True:
            print("maximum energy loss was %6.2f kt/km at %6.2f KM" % (maxenergyloss, elossaltitude))

    def fireiterate(self, r, sigma0):  # m and r should be coupled
        self.sigma0 = sigma0
        v = 19.2E3
        theta = (18.3/180) * np.pi
        z0=97E3
        x=0
        m = (4/3) * np.pi * r**3 * self.pm  # ansatz constraining m with initial spherical geometry
        y0 = [v, m, theta, z0, x, r]
        soln = sp.solve_ivp(self.func, [0, 20000], y0, dense_output=True, rtol = 1e-8, events=[self.ablation, self.crater])
        yfinal = soln.y
        time = soln.t
        t_events = soln.t_events  # times of fracturing and cratering. Note may be more than one fracturing event
        if len(soln.t_events[0]) != 0:  # if an ablation time exists
            ablated = True
        else:
            ablated = False

        v = yfinal[0]
        m = yfinal[1]
        z = yfinal[3]
        r = yfinal[-1]
        energy = 0.5 * m * v**2
        energydiff = np.diff(energy)
        energydiff = np.append(energydiff, 0)
        energydiff = energydiff * -1
        energydiff = energydiff / (4.184*1E12)  # convert kt per kilometer
        z = z/1E3
        energydiff = -1 * energydiff / np.append(np.diff(z), 1)
        # returns energy diff for use in fitting function.
        maxenergydiff = energydiff.max()
        maxenergylossalt = z[energydiff.argmax()]
        return z, energydiff, v, m, r  # decide whether this is useful or not

    def errorplot(self):
        tolerances = np.logspace(-2, -9, num=49)
        ablate = []
        zmetric = []
        y0 = [self.v, self.m, self.theta, self.z0, self.x, self.r]
        for rt in tolerances:
            soln = sp.solve_ivp(self.func, [0, 20000], y0, dense_output=True, rtol=rt, events=[self.ablation, self.crater])
            ablate.append(soln.t_events)
            zmetric.append(soln.y[3])
        ideal = sp.solve_ivp(self.func, [0, 20000], y0, dense_output=True, rtol=10e-11, events=[self.ablation, self.crater])

        idealt = ideal.t_events[1][0]  # time of collision.

        collision = [ablation[1][0] for ablation in ablate]
        corrected = [abs(i - idealt)/i for i in collision]
        corrected = [np.log10(i) for i in corrected]
        tolerances = [np.log10(i) for i in tolerances]

        popt, pcov = so.curve_fit(self.lobf, tolerances, corrected)

        plt.figure()
        plt.xlabel(r'$log_{10}(Solver \ Tolerance)$')
        plt.ylabel(r'$log_{10}(Impact \ time \ error)$')
        plt.xticks(size='small')
        plt.yticks(size='small')
        plt.title(r'Convergence Error vs Solver Tolerance $\sigma_{0} = \inf$')
        plt.plot(tolerances, corrected, 'o')
        plt.plot(tolerances, corrected, label='Actual error')
        plt.plot(tolerances, [self.lobf(i, popt[0], popt[1]) for i in tolerances], label='Best fit through error plot')
        plt.legend()
        plt.savefig('nonfractureerror.png')
        plt.show()
        return tolerances, ablate, idealt, popt

    def lobf(self, x, mm, c):
        return mm * x + c

######################


meteor.crater.terminal = True  # monkey patches the terminal condition into the event function for the main solver



