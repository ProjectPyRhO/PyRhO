r"""
Rhodopsin model class and functions
    * Three-state model :math:`\{C, O, D\}`
    * Four-state model :math:`\{C_1, O_1, O_2, C_2\}`
    * Six-state model :math:`\{C_1, I_1, O_1, O_2, I_2, C_2\}`
"""

import warnings
import logging
import abc
import itertools

import matplotlib.pyplot as plt
#import matplotlib as mpl # For plotStates
import numpy as np
from scipy.integrate import odeint

from pyrho.utilities import calcV1
from pyrho.parameters import PyRhOobject, modelParams, stateLabs, rhoType
from pyrho import config

__all__ = ['models', 'selectModel']

logger = logging.getLogger(__name__)

# Model class definitions #


class RhodopsinModel(PyRhOobject):
    """Common base class for all models."""
    # This an abstract base class since it is never directly instantiated
    __metaclass__ = abc.ABCMeta
    # TODO: Revise to be stateless and store data in PhotoCurrent objects
    phi = 0.0  # Instantaneous Light flux [photons * mm^-2 * s^-1]

    def __init__(self, params=None, rhoType=rhoType):

        if params is None:
            try:
                params = modelParams[str(self.nStates)]
            except:
                pass
        self.rhoType = rhoType  # E.g. 'ChR2' or 'ArchT'

        self.setParams(params)

        # Ensure v1 is scaled correctly so that f(V=-70) = 1
        v1 = calcV1(self.E, self.v0)
        if not np.isclose(self.v1, v1, rtol=1e-3, atol=1e-5):
            warnings.warn(f"Correcting v1 scaling: {self.v1} <-- {v1}")
            self.v1 = v1

        self.initStates(phi=self.phi_0, s0=self.s_0)
        # self.transRates = {r: getattr(self, r) for r in itertools.chain(self.photoRates, self.constRates)}

        if config.verbose > 1:
            print(f"PyRhO {self.nStates}-state {self.rhoType} model initialised!")

        if config.verbose > 2:
            self.printParams()

    def __str__(self):
        return f"{stateLabs[self.nStates]}-state {self.rhoType}"  # self.__name__+
        # return self.brian_phi_t

    def __repr__(self):
        return f"<PyRhO {stateLabs[self.nStates]}-state {self.rhoType} Model object>"

    def __call__(self):
        """When a rhodopsin is called, return its internal state at that instant."""
        return self.calcI(self.V, self.states[-1, :])

    def storeStates(self, soln, t):
        # TODO: Make array dimensions consistent e.g. both row vectors
        self.states = np.vstack((self.states, soln))  # np.append(self.states, soln, axis=0)
        self.t = np.hstack((self.t, t))  # np.append(self.t, t, axis=1)
        #self.pulseInd = np.append(self.pulseInd, _idx_pulses_, axis=0)

    def getStates(self):
        """Returns (states, t)."""
        return self.states, self.t

    def getRates(self):
        """Returns an ordered dictionary of all transition rates."""
        return {r: getattr(self, r) 
                for r in itertools.chain(self.photoRates, self.constRates)}

    def reportState(self):
        self.dispRates()

    def initStates(self, phi, s0=None):
        """Clear state arrays and set transition rates."""
        if s0 is None:
            s0 = self.s_0
        #assert(len(s0) == self.nStates)
        #if self.nStates=='6K':
        #    self.states = np.vstack((np.empty([0, 6]), s0))
        # 
        #else:
        #    self.states = np.vstack((np.empty([0, self.nStates]), s0))
        assert len(s0) == self.nStates
        self.states = np.vstack((np.empty([0, self.nStates]), s0))
        self.t = [0]
        self.pulseInd = np.empty([0, 2], dtype=int)  # Light on and off indexes for each pulse
        self.ssInf = []
        self.setLight(phi)


    # Implement this universal function
    def calcI(self, V, states=None):
        r"""Takes Voltage [mV] and open state variables {O} or {O1, O2} to
        calculate current [nA].

        By convention:
        Negative current: positive ions entering the cell (e.g. Na^+ influx).
        Positive current: positive ions exiting (or negative ions entering) the cell (e.g. Cl^- in or K^+ out).
        """

        if states is None:
            states = self.states

        g_RhO = self.g0 * self.calcfphi(states) * self.calcfV(V)
        I_RhO = g_RhO * (V - self.E)  # Photocurrent: (pS * mV)
        return I_RhO * (1e-6)  # 10^-12 * 10^-3 * 10^-6 (nA)

    def calcfV(self, V):
        """Calculate the voltage-dependent conductance scaling factor, f(v)."""
        # TODO: Try this solution
        # http://stackoverflow.com/questions/26248654/numpy-return-0-with-divide-by-zero
        #with np.errstate(divide='ignore', invalid='ignore'):
        #    fV = (self.v1/(V-self.E))*(1-np.exp(-(V-self.E)/self.v0))
        #    fV[~ np.isfinite(fV)] = self.v1/self.v0  # -inf inf NaN

        #with np.errstate(divide='ignore', invalid='ignore'):
        #    c = np.true_divide(a,b)
        #    c[c == np.inf] = 0
        #    c = np.nan_to_num(c)

        if self.v0 == 0:
            raise ZeroDivisionError("f(V) undefined for v0 = 0")
        try:
            fV = (self.v1/(V-self.E))*(1-np.exp(-(V-self.E)/self.v0))  # Dimensionless
        except ZeroDivisionError:
            if np.isscalar(V):
                if np.isclose(V, self.E):  # (V == self.E):
                    fV = self.v1/self.v0
            else:  # type(fV) in (tuple, list, array)
                # fV[np.isnan(fV)] = self.v1/self.v0 # Fix the error when dividing by zero
                fV[np.isclose(V-self.E, np.zeros_like(V))] = self.v1/self.v0
        #else:
        #    fV=1 ### Extend to vector
        return fV

    #@property
    def calcIss(self, V):
        """Calculate the steady-state current for a given voltage (and model parameters).
        """
        return self.calcI(V, states=self.calcSteadyState())

    '''
    @property
    def T(self):
        return self._T
    @T.setter
    def T(self, T):
        self._T = T
    '''

    def calcfT(self, T):
        raise NotImplementedError

    def calcfpH(self, pH):
        raise NotImplementedError

    def calcSoln(self, t, s0=None):
        if s0 is None:
            s0 = self.s_0
        return odeint(self.solveStates, s0, t, args=(None,), Dfun=self.jacobian)

    def plotActivation(self, actFunc, label=None, phis=np.logspace(12, 21, 1001), ax=None):
        if ax is None:
            ax = plt.gca()
        else:
            plt.sca(ax)

        if label is not None:
            ax.plot(phis, actFunc(phis), label=label)
            ax.legend()
        else:
            ax.plot(phis, actFunc(phis))
        ax.set_xscale('log')
        ax.set_xlabel(r'$\phi \ \mathrm{[photons \cdot mm^{-2} \cdot s^{-1}]}$')
        ax.set_ylabel(r'$\mathrm{Transition\ rate \ [ms^{-1}]}$')
        return

    def plotRates(self, phis=np.logspace(12, 21, 1001), logscale='both'):
        """Plot all activation rates."""
        rFig, axR = plt.subplots()

        for r, l in zip(self.photoFuncs, self.photoLabels):
            #self.plotActivation(getattr(self, r), label=l, phis=phis, ax=axR)
            axR.plot(phis, getattr(self, r)(phis), label=l)
            #for phi in phis:
            #    self.setLight(phi)
            #    axR.plot()

        for r, l in zip(self.constRates, self.constLabels):
            #axR.axhline(y=r, label=l)
            axR.plot(phis, np.ones(len(phis)) * self.__dict__[r], label=l, ls='--')

        if logscale == 'both' or logscale == 'x':
            axR.set_xscale('log')
        if logscale == 'both' or logscale == 'y':
            axR.set_yscale('log')

        axR.set_xlabel(r'$\phi \ \mathrm{[photons \cdot mm^{-2} \cdot s^{-1}]}$')
        axR.set_ylabel(r'$\mathrm{Transition\ rate \ [ms^{-1}]}$')
        axR.legend(loc='best')

        # Tick locators take too long with log scales
        #axR.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        #axR.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
        axR.grid(visible=True, which='minor', axis='both', linewidth=.2)
        axR.grid(visible=True, which='major', axis='both', linewidth=1)

        plt.tight_layout()
        plt.show()

        return


class RhO_3states(RhodopsinModel):
    """Class definition for the 3-state model."""

    # Class attributes
    nStates = 3
    useAnalyticSoln = True

    s_0 = np.array([1, 0, 0])               # Default: Initialise in dark
    phi_0 = 0.0                             # Default flux level in dark-adapted state
    stateVars = ['C', 'O', 'D']             # stateVars[0] is the 'ground' state
    stateLabels = ['$C$', '$O$', '$D$']
    photoFuncs = ['_calcGa', '_calcGr']     # {'Ga':'_calcGa', 'Gr':'_calcGr'} --> photoRates['Ga'](phi)
    photoRates = ['Ga', 'Gr']
    photoLabels = ['$G_a$', '$G_r$']        # {'Ga':'$G_a$', 'Gr':'$G_r$'}
    constRates = ['Gd']
    constLabels = ['$G_d$']

    paramsList = ['g0', 'phi_m', 'k_a', 'p', 'Gd', 'Gr0', 'k_r', 'q',
                  'E', 'v0', 'v1']  # List of model constants

    # TODO: Create lists of parameters, variables and constants
    connect = [[0, 1, 0],
               [0, 0, 1],
               [1, 0, 0]]

    conDir  = [[ 0, -1,  1],
               [ 1,  0, -1],
               [-1,  1,  0]]

    # TODO: Try using \begin{align*} now this is a raw string, or leave a comment if it does not work
    equations = r"""
                $$ \dot{C} = G_{r}(\phi)D - G_{a}(\phi)C $$
                $$ \dot{O} = G_{a}(\phi)C - G_{d}O $$
                $$ \dot{D} = G_{d}O - G_{r}(\phi)D $$
                $$ C + O + D = 1 $$
                $$ $$
                $$ G_a(\phi) = k_a\frac{\phi^p}{\phi^p + \phi_m^p} $$
                $$ G_r(\phi) = k_r\frac{\phi^q}{\phi^q + \phi_m^q} + G_{r0} $$
                $$ $$
                $$ f_{\phi}(\phi) = O \qquad \qquad $$
                $$ f_v(v) = v_1\frac{1-e^{-(v-E)/v_0}}{(v-E)} $$
                $$ I_{\phi} = g_0 \cdot f_{\phi}(\phi) \cdot f_v(v) \cdot (v-E) $$
                """

    _latex =   [r"$\dot{C} = G_{r}(\phi)D - G_{a}(\phi)C$",
                r"$\dot{O} = G_{a}(\phi)C - G_{d}O$",
                r"$\dot{D} = G_{d}O - G_{r}(\phi)D$",
                r"$C + O + D = 1$",
                r"$G_a(\phi) = k\frac{\phi^p}{\phi^p + \phi_m^p}$",
                r"$G_r(\phi) = \mathcal{H}(\phi) \cdot G_{r1} + G_{r0}$",
                r"$f_{\phi}(\phi) = O$"
                r"$f_v(v) = \frac{1-\exp({-(v-E)/v_0})}{(v-E)/v_1}$",
                r"$I_{\phi} = g_0 \cdot f_{\phi}(\phi) \cdot f_v(v) \cdot (v-E)$"]

    eqIss = r"""$I_{SS} = \bar{g_0} \cdot \frac{G_a \cdot G_r}{G_d \cdot (G_r + G_a) + G_a \cdot G_r} \cdot (v-E)
    = \bar{g_0} \cdot \frac{\tau_d}{\tau_d + \tau_r + \tau_\phi} \cdot (v-E)$"""

    brianStateVars = ['C', 'O', 'D']

    brian = '''
            dC/dt = Gr*D - Ga*C                     : 1
            dO/dt = Ga*C - Gd*O                     : 1
            dD/dt = Gd*O - Gr*D                     : 1
            Ga = Theta*k_a*((phi**p)/(phi**p + phi_m**p))             : second**-1
            Gr = Theta*k_r*((phi(t)**q)/(phi(t)**q + phi_m**q)) + Gr0 : second**-1
            f_phi = O                               : 1
            f_v = (1-exp(-(v-E)/v0))/((v-E)/v1)     : 1
            I = g0*f_phi*f_v*(v-E)                  : amp
            phi                                     : metre**-2*second**-1 (shared)
            Theta = int(phi > 0*phi)                : 1 (shared)
            '''
    #S_D = 1 - S_C - S_O : 1 #
    #int(phi>0)
    #Ga = stimulus*k*((phi**p)/(phi**p+phi_m**p)) : second**-1
    #H = stimulus*((phi**p)/(phi**p+phi_m**p)) : 1
    #Ga = k * H : second**-1
    #stimulus = int(ceil(clip(phi(t), 0, 1)))

    brian_phi_t = '''
            dC/dt = Gr*D - Ga*C                                       : 1
            dO/dt = Ga*C - Gd*O                                       : 1
            dD/dt = Gd*O - Gr*D                                       : 1
            Ga = Theta*k_a*((phi(t)**p)/(phi(t)**p + phi_m**p))       : second**-1
            Gr = Theta*k_r*((phi(t)**q)/(phi(t)**q + phi_m**q)) + Gr0 : second**-1
            f_phi = O                                                 : 1
            f_v = (1-exp(-(v-E)/v0))/((v-E)/v1)                       : 1
            I = g0*f_phi*f_v*(v-E)                                    : amp
            Theta = int(phi(t) > 0*phi(t))                            : 1 (shared)
            '''

    """
    @property
    def Ga(self):
        return self._calcGa(self.phi)

    @Ga.setter
    def Ga(self, value):
        return

    @property
    def Gr(self):
        return self._calcGr(self.phi)
    """

    def _calcGa(self, phi):
        return self.k_a * phi**self.p/(phi**self.p + self.phi_m**self.p)

    def _calcGr(self, phi):
        #return self.Gr0 + (phi>0)*self.Gr1
        return self.Gr0 + self.k_r * phi**self.q/(phi**self.q + self.phi_m**self.q)
        #return self.Gr0 + self.Gr1 * np.log(1 + phi/self.phi0) # self.Gr0 + self.Gr1
        #return self.Gr_dark + self.Gr_light * np.log(1 + phi/self.phi0) # self.Gr0 + self.Gr1
        # return 1/(taur_dark*exp(-log(1+phi/phi0))+taur_min) # Fig 6 Nikolic et al. 2009
        ### return Gr_dark + kr*(1-exp(-phi/phi0)) # a = Gr_max - Gr_dark

    #def set_Gd(self, phi):
    #    return 1/(taud_dark - a*log(1+phi/phi0)) # Fig 6 Nikolic et al. 2009

    def setLight(self, phi):
        """Set transition rates according to the instantaneous photon flux density"""
        if phi < 0:
            phi = 0
        self.phi = phi
        self.Ga = self._calcGa(phi)
        self.Gr = self._calcGr(phi)
        if config.verbose > 1:
            self.dispRates()

    def dispRates(self):
        print(f"Transition rates (phi={self.phi:.3g}):\n"
            f"C --[Ga={self.Ga:.3g}]--> O --[Gd={self.Gd:.3g}]--> D --[Gr={self.Gr:.3g}]--> C")

    def solveStates(self, s_0, t, phi_t=None):
        """Differential equations of the 3-state model to be solved by odeint"""
        # Add interpolation of values for phi(t) to initialisation f_phi = interp1d(t,sin(w*t),kind='cubic')
        # Then pass as an argument to integrator: odeint(func, y0, t, args=())

        if phi_t is not None:
            self.setLight(float(phi_t(t)))
        C, O, D = s_0  # Split state vector into individual variables
        dCdt = -self.Ga*C +             self.Gr*D   # C'
        dOdt =  self.Ga*C - self.Gd*O               # O'
        dDdt =              self.Gd*O - self.Gr*D   # D'
        #f2 = -(f0+f1)
        return np.array([dCdt, dOdt, dDdt])

    def jacobian(self, s_0, t, phi_t=None):  # jacobianPhi_t
        """
        Jacobian matrix used to improve precision / speed up ODE solver
        jac[i,j] = df[i]/dy[j]; where y'(t) = f(t,y)
        """
        #self.setLight(phi_t(t))
        return np.array([[-self.Ga, 0, self.Gr],  # [dCdt/dC, dCdt/dO, dCdt/dD]
                         [self.Ga, -self.Gd, 0],  # [dOdt/dC, dOdt/dO, dOdt/dD]
                         [0, self.Gd, -self.Gr]]) # [dDdt/dC, dDdt/dO, dDdt/dD]

    def hessian(self, s_0, t):
        """
        Hessian matrix for scipy.optimize.minimize.
        (Only for Newton-CG, dogleg, trust-ncg.)
        H(f)_ij(X) = D_iD_jf(X)
        """
        # return np.array([[0, 0, 0],
        #                  [0, 0, 0],
        #                  [0, 0, 0]])
        return np.zeros((3, 3))

    def calcfphi(self, states=None):
        if states is None:
            states = self.states
        C, O, D = states.T
        return O

    # def calcOn(self,t):
        # """Calculate the on phase current for square light pulses from the analytic solution"""
        # r = np.array([lam1, lam2])
        # k = np.array([a1, a2])
        # I = k * np.exp(-r*t)
        # -(a0 + a1*(1-np.exp(-t/tau_act)) + a2*np.exp(-t/tau_deact))
        # pass

    # def calcOff():
        # """Calculate the off phase current for square light pulses from the analytic solution"""
        # -(A*np.exp(-Gd*t))
        # pass

    def calcSteadyState(self, phi):
        self.setLight(phi)
        denom = self.Gd * (self.Gr + self.Ga) + self.Ga * self.Gr
        Css = self.Gd * self.Gr
        Oss = self.Ga * self.Gr
        Dss = self.Ga * self.Gd
        self.steadyStates = np.array([Css, Oss, Dss]) / denom
        return self.steadyStates

    def calcSoln(self, t, s0=None):
        if s0 is None:
            s0 = self.s_0
        [C_0, O_0, D_0] = s0
        Ga = self.Ga
        Gd = self.Gd
        Gr = self.Gr

        t = t - t[0]  # Shift time array forwards or backwards to start at 0

        SP = Ga*Gd + Ga*Gr + Gd*Gr
        SQ = Ga**2 + Gd**2 + Gr**2
        if 2 * SP > SQ:
            if config.verbose > 1:
                print(f'Imaginary solution! {SP = }; {SQ = }'
                      ' --> (SQ-2*SP)**(1/2) = NaN')
            return odeint(self.solveStates, s0, t, Dfun=self.jacobian)
            #raise ValueError()  # Uncomment this when error catching is implemented

        RSD = (SQ - 2 * SP)**(1/2)  # xi
        lambda_1 = (Ga + Gd + Gr + RSD)/2
        lambda_2 = (Ga + Gd + Gr - RSD)/2
        Z_1 = C_0 * Gd * Ga + O_0 * (Gd * (Ga - lambda_1)) + D_0 * Gr * (Gr-lambda_2)
        Z_2 = C_0 * Gd * Ga + O_0 * (Gd * (Ga - lambda_2)) + D_0 * Gr * (Gr-lambda_1)
        Exp_1 = np.exp(-t * lambda_1)
        Exp_2 = np.exp(-t * lambda_2)

        C = (Z_1 * lambda_2 * (lambda_1-Gd-Gr) * Exp_1 - Z_2 * lambda_1 * (lambda_2-Gd-Gr) * Exp_2
            + (RSD * Gd**2 * Gr * (C_0 + D_0 + O_0))) / (Gd * SP * RSD)
        O = (-Z_1 * lambda_2 * (lambda_1-Gr) * Exp_1 + Z_2 * lambda_1 * (lambda_2-Gr) * Exp_2
            + (RSD * Ga * Gd * Gr * (C_0+D_0+O_0))) / (Gd * SP * RSD)
        D = (Z_1 * lambda_2 * Exp_1 - Z_2 * lambda_1 * Exp_2
            + (RSD * Gd * Ga * (C_0+D_0+O_0))) / (SP * RSD)

        return np.column_stack((C, O, D))


class RhO_4states(RhodopsinModel):
    """Class definition for the 4-state model."""

    # Class attributes
    nStates = 4
    useAnalyticSoln = False

    phi_0 = 0.0                             # Instantaneous Light flux
    s_0 = np.array([1, 0, 0, 0])            # Default: Initialise in the dark
    stateVars = ['C1', 'O1', 'O2', 'C2']    # stateVars[0] is the 'ground' state
    stateLabels = ['$C_1$', '$O_1$', '$O_2$', '$C_2$']  # [texIt(s) for s in stateVars]

    photoFuncs = ['_calcGa1', '_calcGa2', '_calcGf', '_calcGb']
    photoRates = ['Ga1', 'Ga2', 'Gf', 'Gb']
    photoLabels = ['$G_{a1}$', '$G_{a2}$', '$G_{f}$', '$G_{b}$', '$G_{d1}$', '$G_{d2}$']
    constRates = ['Gd1', 'Gd2', 'Gr0']
    constLabels = ['$G_{d1}$', '$G_{d2}$', '$G_{r0}$']

    paramsList = ['g0', 'gam', 'phi_m', 'k1', 'k2', 'p',
                  'Gf0', 'k_f', 'Gb0', 'k_b', 'q', 'Gd1', 'Gd2', 'Gr0',
                  'E', 'v0', 'v1']  # List of model constants

    connect = [[0, 1, 0, 0],
               [1, 0, 1, 0],
               [0, 1, 0, 1],
               [1, 0, 1, 0]]

    equations = r"""
                $$ \dot{C_1} = G_{d1}O_1 + G_{r0}C_2 - G_{a1}(\phi)C_1 $$
                $$ \dot{O_1} = G_{a1}(\phi)C_1 + G_{b}(\phi)O_2 - (G_{d1} + G_{f}(\phi))O_1 $$
                $$ \dot{O_2} = G_{a2}(\phi)C_2 + G_{f}(\phi)O_1 - (G_{d2} + G_{b}(\phi))O_2 $$
                $$ \dot{C_2} = G_{d2}O_2 - (G_{r0} + G_{a2}(\phi))C_2 $$
                $$ C_1 + O_1 + O_2 + C_2 = 1 $$
                $$$$
                $$ G_{a1}(\phi) = k_1 \frac{\phi^p}{\phi^p + \phi_m^p} $$
                $$ G_{f}(\phi)  = k_{f} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{f0} $$
                $$ G_{b}(\phi)  = k_{b} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{b0} $$
                $$ G_{a2}(\phi) = k_2 \frac{\phi^p}{\phi^p + \phi_m^p} $$
                $$$$
                $$ f_{\phi}(\phi) = O_1+\gamma O_2 $$
                $$ f_v(v) = v_1\frac{1-e^{-(v-E)/v_0}}{(v-E)} $$
                $$ I_{\phi} = g_0 \cdot f_{\phi}(\phi) \cdot f_v(v) \cdot (v-E) $$
                """

    brianStateVars = ['C_1', 'O_1', 'O_2', 'C_2']

    brian = '''
            dC_1/dt = Gd1*O_1 + Gr0*C_2 - Ga1*C_1     : 1
            dO_1/dt = Ga1*C_1 + Gb*O_2 - (Gd1+Gf)*O_1 : 1
            dO_2/dt = Ga2*C_2 + Gf*O_1 - (Gd2+Gb)*O_2 : 1
            C_2 = 1 - C_1 - O_1 - O_2                 : 1
            H_p = Theta*((phi**p)/(phi**p+phi_m**p))  : 1
            H_q = Theta*((phi**q)/(phi**q+phi_m**q))  : 1
            Ga1 = k1*H_p                              : second**-1
            Ga2 = k2*H_p                              : second**-1
            Gf  = k_f*H_q + Gf0                       : second**-1
            Gb  = k_b*H_q + Gb0                       : second**-1
            f_v   = (1-exp(-(v-E)/v0))/((v-E)/v1)     : 1
            f_phi = O_1+gam*O_2                       : 1
            I     = g0*f_phi*f_v*(v-E)                : amp
            phi                                       : metre**-2*second**-1 (shared)
            Theta = int(phi > 0*phi)                  : 1 (shared)
            '''

    brian_phi_t = '''
            dC_1/dt = Gd1*O_1 + Gr0*C_2 - Ga1*C_1                  : 1
            dO_1/dt = Ga1*C_1 + Gb*O_2 - (Gd1+Gf)*O_1              : 1
            dO_2/dt = Ga2*C_2 + Gf*O_1 - (Gd2+Gb)*O_2              : 1
            C_2 = 1 - C_1 - O_1 - O_2                              : 1
            Theta = int(phi(t) > 0*phi(t))                         : 1 (shared)
            Ga1   = Theta*k1*phi(t)**p/(phi(t)**p+phi_m**p)        : second**-1
            Ga2   = Theta*k2*phi(t)**p/(phi(t)**p+phi_m**p)        : second**-1
            Gf    = Theta*k_f*phi(t)**q/(phi(t)**q+phi_m**q) + Gf0 : second**-1
            Gb    = Theta*k_b*phi(t)**q/(phi(t)**q+phi_m**q) + Gb0 : second**-1
            f_v   = (1-exp(-(v-E)/v0))/((v-E)/v1)                  : 1
            f_phi = O_1+gam*O_2                                    : 1
            I     = g0*f_phi*f_v*(v-E)                             : amp
            '''

    def _calcGa1(self, phi):
        # N.B. making Ga a function of time (as in Appendix 1) results in the Six-state model
        # Gai = ei * F * f(t,tChR) See App 1
        # Ga = 1/tauChR
        #e = 0.5
        #sigma_ret = 1.2e-8 * 1e-6 # Convert from m^2 to mm^2
        #w_loss = 1.1
        #return self.k1 * phi/self.phi0 #e*phi*sigma_ret / w_loss
        #return self.k1 * (1-np.exp(-phi/self.phi0)) #e*phi*sigma_ret / w_loss
        return self.k1 * phi**self.p/(phi**self.p + self.phi_m**self.p)

    def _calcGa2(self, phi):
        #e = 0.15
        #sigma_ret = 1.2e-8  * 1e-6 # Convert from m^2 to mm^2
        #w_loss = 1.1
        #return self.k2 * phi/self.phi0 #e*phi*sigma_ret / w_loss
        #return self.k2 * (1-np.exp(-phi/self.phi0))
        return self.k2 * phi**self.p/(phi**self.p + self.phi_m**self.p)

    def _calcGf(self, phi):
        #return self.e12d + self.c1*np.log(1+(phi/self.phi0)) # e12(phi=0) = e12d
        return self.Gf0 + self.k_f * phi**self.q/(phi**self.q + self.phi_m**self.q)

    def _calcGb(self, phi):
        #return self.e21d + self.c2*np.log(1+(phi/self.phi0)) # e21(phi=0) = e21d
        return self.Gb0 + self.k_b * phi**self.q/(phi**self.q + self.phi_m**self.q)

    def setLight(self, phi):
        """
        Set transition rates according to the instantaneous photon flux density
        """
        if phi < 0:
            phi = 0
        self.phi = phi
        self.Ga1 = self._calcGa1(phi)
        self.Ga2 = self._calcGa2(phi)
        self.Gf = self._calcGf(phi)
        self.Gb = self._calcGb(phi)
        if config.verbose > 1:
            self.dispRates()

    def dispRates(self):
        print(f"Transition rates (phi={self.phi:.3g}):\n"
            f"C1 --[Ga1={self.Ga1:.3g}]--> O1 --[Gf={self.Gf:.3g}]--> O2\n"
            f"O1 <--[Gb={self.Gb:.3g}]-- O2 <--[Ga2={self.Ga2:.3g}]-- C2")

    def solveStates(self, s_0, t, phi_t=None):
        """Differential equations of the 4-state model to be solved by odeint"""
        if phi_t is not None:
            self.setLight(float(phi_t(t)))
        C1, O1, O2, C2 = s_0 # Split state vector into individual variables s1=s[0], s2=s[1], etc
        dC1dt = -self.Ga1*C1 +           self.Gd1*O1       +                  self.Gr0 *C2  # C1'
        dO1dt =  self.Ga1*C1 - (self.Gd1+self.Gf)*O1       +   self.Gb*O2                   # O1'
        dO2dt =                           self.Gf*O1 - (self.Gd2+self.Gb)*O2 + self.Ga2*C2  # O2'
        dC2dt =                                       self.Gd2*O2 - (self.Ga2+self.Gr0)*C2  # C2'
        #f3 = -(f0+f1+f2)
        return np.array([dC1dt, dO1dt, dO2dt, dC2dt])

    def jacobian(self, s_0, t, phi_t=None):
        """
        Jacobian matrix used to improve precision / speed up ODE solver.
        jac[i,j] = df[i]/dy[j]; where y'(t) = f(t,y)
        """
        return np.array([[-self.Ga1, self.Gd1, 0, self.Gr0],
                         [self.Ga1, -(self.Gd1+self.Gf), self.Gb, 0],
                         [0, self.Gf, -(self.Gd2+self.Gb), self.Ga2],
                         [0, 0, self.Gd2, -(self.Ga2+self.Gr0)]])

    def hessian(self, s_0, t):
        """
        Hessian matrix for scipy.optimize.minimize.
        (Only for Newton-CG, dogleg, trust-ncg.)
        H(f)_ij(X) = D_iD_jf(X)
        """
        return np.zeros((4, 4))

    def calcfphi(self, states=None):
        if states is None:
            states = self.states
        C1, O1, O2, C2 = states.T
        gam = self.gam
        return O1 + gam * O2

    def calcSteadyState(self, phi):
        self.setLight(phi)
        Ga1 = self.Ga1
        Ga2 = self.Ga2
        Gr0 = self.Gr0
        Gd1 = self.Gd1
        Gd2 = self.Gd2
        Gf = self.Gf
        Gb = self.Gb
        denom = Ga1 * (Gf * (Gr0 + Gd2 + Ga2) + Gb * (Gr0 + Ga2) + Gd2 * Gr0) \
              + Gd1 * (Gb * (Gr0 + Ga2) + Gd2 * Gr0) + Gf * Gd2 * Gr0
        C1ss = (Gd1 * (Gb * (Gr0 + Ga2) + Gd2 * Gr0) + Gf * Gd2 * Gr0)
        O1ss = (Ga1 * (Gb * (Gr0 + Ga2) + Gd2 * Gr0))
        O2ss = (Gf * Ga1 * (Gr0 + Ga2))
        C2ss = (Gf * Ga1 * Gd2)
        self.steadyStates = np.array([C1ss, O1ss, O2ss, C2ss]) / denom
        return self.steadyStates

    def calcSoln(self, t, s0=None):
        raise NotImplementedError(self.nStates)


class RhO_6states(RhodopsinModel):
    """Class definition for the 6-state model."""

    # Class attributes
    nStates = 6
    useAnalyticSoln = False
    s_0 = np.array([1, 0, 0, 0, 0, 0])  # [s1_0=1, s2_0=0, s3_0=0, s4_0=0, s5_0=0, s6_0=0] # array not necessary
    phi_0 = 0.0                         # Default initial flux
    stateVars = ['C1', 'I1', 'O1', 'O2', 'I2', 'C2']  # stateVars[0] is the 'ground' state
    stateLabels = ['$C_1$', '$I_1$', '$O_1$', '$O_2$', '$I_2$', '$C_2$']
    photoFuncs = ['_calcGa1', '_calcGa2', '_calcGf', '_calcGb']
    photoRates = ['Ga1', 'Ga2', 'Gf', 'Gb']
    photoLabels = ['$G_{a1}$', '$G_{a2}$', '$G_{f}$', '$G_{b}$', '$G_{d1}$', '$G_{d2}$']
    constRates = ['Go1', 'Go2', 'Gd1', 'Gd2', 'Gr0']
    constLabels = ['$G_{o1}$', '$G_{o2}$', '$G_{d1}$', '$G_{d2}$', '$G_{r0}$']

    paramsList = ['g0', 'gam', 'phi_m', 'k1', 'k2', 'p',
                  'Gf0', 'k_f', 'Gb0', 'k_b', 'q', 'Go1', 'Go2',
                  'Gd1', 'Gd2', 'Gr0', 'E', 'v0', 'v1']  # List of model constants

    connect = [[0, 1, 0, 0, 0, 0],  # s_1 --> s_i=1...6
               [0, 0, 1, 0, 0, 0],  # s_2 -->
               [1, 0, 0, 1, 0, 0],
               [0, 0, 1, 0, 0, 1],
               [0, 0, 0, 1, 0, 0],
               [1, 0, 0, 0, 1, 0]]

    equations = r"""
                $$ \dot{C_1} = G_{d1}O_1 + G_{r0}C_2 - G_{a1}(\phi)C_1 $$
                $$ \dot{I_1} = G_{a1}(\phi)C_1 - G_{o1}I_1 $$
                $$ \dot{O_1} = G_{o1}I_1 + G_{b}(\phi)O_2 - (G_{d1} + G_{f}(\phi))O_1 $$
                $$ \dot{O_2} = G_{o2}I_2 + G_{f}(\phi)O_1 - (G_{d2} + G_{b}(\phi))O_2 $$
                $$ \dot{I_2} = G_{a2}(\phi)C_2 - G_{o2}I_2 $$
                $$ \dot{C_2} = G_{d2}O_2 - (G_{r0} + G_{a2}(\phi))C_2 $$
                $$ C_1 + I_1 + O_1 + O_2 + I_2 + C_2 = 1 $$
                $$$$
                $$ G_{a1}(\phi) = k_{1} \frac{\phi^p}{\phi^p + \phi_m^p} $$
                $$ G_{f}(\phi)  = k_{f} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{f0} $$
                $$ G_{b}(\phi)  = k_{b} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{b0} $$
                $$ G_{a2}(\phi) = k_{2} \frac{\phi^p}{\phi^p + \phi_m^p} $$
                $$$$
                $$ f_{\phi}(\phi) = O_1+\gamma O_2 $$
                $$ f_v(v) = v_1\frac{1-e^{-(v-E)/v_0}}{(v-E)} $$
                $$ I_{\phi} = g_0 \cdot f_{\phi}(\phi) \cdot f_v(v) \cdot (v-E) $$
                """

    brianStateVars = ['C_1', 'I_1', 'O_1', 'O_2', 'I_2', 'C_2']

    brian = '''
            dC_1/dt = Gd1*O_1 + Gr0*C_2 - Ga1*C_1     : 1
            dI_1/dt = Ga1*C_1 - Go1*I_1               : 1
            dO_1/dt = Go1*I_1 + Gb*O_2 - (Gd1+Gf)*O_1 : 1
            dO_2/dt = Go2*I_2 + Gf*O_1 - (Gd2+Gb)*O_2 : 1
            dI_2/dt = Ga2*C_2 - Go2*I_2               : 1
            C_2 = 1 - C_1 - I_1 - O_1 - O_2 - I_2     : 1
            H_p = Theta*((phi**p)/(phi**p+phi_m**p))  : 1
            H_q = Theta*((phi**q)/(phi**q+phi_m**q))  : 1
            Ga1 = k1*H_p                              : second**-1
            Ga2 = k2*H_p                              : second**-1
            Gf = k_f*H_q + Gf0                        : second**-1
            Gb = k_b*H_q + Gb0                        : second**-1
            f_v = (1-exp(-(v-E)/v0))/((v-E)/v1)       : 1
            f_phi = O_1+gam*O_2                       : 1
            I = g0*f_phi*f_v*(v-E)                    : amp
            phi                                       : metre**-2*second**-1 (shared)
            Theta = int(phi > 0*phi)                  : 1 (shared)
            '''

    brian_phi_t = '''
            dC_1/dt = Gd1*O_1 + Gr0*C_2 - Ga1*C_1            : 1
            dI_1/dt = Ga1*C_1 - Go1*I_1                      : 1
            dO_1/dt = Go1*I_1 + Gb*O_2 - (Gd1+Gf)*O_1        : 1
            dO_2/dt = Go2*I_2 + Gf*O_1 - (Gd2+Gb)*O_2        : 1
            dI_2/dt = Ga2*C_2 - Go2*I_2                      : 1
            C_2 = 1 - C_1 - I_1 - O_1 - O_2 - I_2            : 1
            Theta = int(phi(t) > 0*phi(t))                   : 1 (shared)
            H_p   = Theta*((phi(t)**p)/(phi(t)**p+phi_m**p)) : 1
            H_q   = Theta*((phi(t)**q)/(phi(t)**q+phi_m**q)) : 1
            Ga1   = k1*H_p                                   : second**-1
            Ga2   = k2*H_p                                   : second**-1
            Gf    = k_f*H_q + Gf0                            : second**-1
            Gb    = k_b*H_q + Gb0                            : second**-1
            f_v   = (1-exp(-(v-E)/v0))/((v-E)/v1)            : 1
            f_phi = O_1+gam*O_2                              : 1
            I     = g0*f_phi*f_v*(v-E)                       : amp
            '''

    def _calcGa1(self, phi):
        #return self.a10*(phi/self.phi0)
        return self.k1 * phi**self.p/(phi**self.p + self.phi_m**self.p)

    def _calcGf(self, phi):
        #return self.a30 + self.a31*np.log(1+(phi/self.phi0))
        return self.Gf0 + self.k_f * phi**self.q/(phi**self.q + self.phi_m**self.q)

    def _calcGb(self, phi):
        #return self.b20 + self.b21*np.log(1+(phi/self.phi0))
        return self.Gb0 + self.k_b * phi**self.q/(phi**self.q + self.phi_m**self.q)

    def _calcGa2(self, phi):
        #return self.b40*(phi/self.phi0)
        return self.k2 * phi**self.p/(phi**self.p + self.phi_m**self.p)

    def setLight(self, phi):
        if phi < 0:
            phi = 0
        self.phi = phi
        self.Ga1 = self._calcGa1(phi)
        self.Gf = self._calcGf(phi)
        self.Gb = self._calcGb(phi)
        self.Ga2 = self._calcGa2(phi)
        if config.verbose > 1:
            self.dispRates()

    def dispRates(self):
        """Print photosensitive transition rates."""
        print(f"Transition rates (phi={self.phi:.3g}):\n"
            f"C1 --[Ga1={self.Ga1:.3g}]--> O1 --[Gf={self.Gf:.3g}]--> O2\n"
            f"O1 <--[Gb={self.Gb:.3g}]-- O2 <--[Ga2={self.Ga2:.3g}]-- C2")

    def solveStates(self, s_0, t, phi_t=None):
        """Differential equations of the 6-state model to be solved by odeint."""
        if phi_t is not None:
            self.setLight(float(phi_t(t)))
        C1, I1, O1, O2, I2, C2 = s_0  # Unpack state vector
        dC1dt = -self.Ga1*C1 + self.Gd1*O1 + self.Gr0*C2
        dI1dt =  self.Ga1*C1 - self.Go1*I1
        dO1dt =  self.Go1*I1 - (self.Gd1+self.Gf)*O1 + self.Gb*O2
        dO2dt =  self.Gf*O1 - (self.Gb+self.Gd2)*O2 + self.Go2*I2
        dI2dt = -self.Go2*I2 + self.Ga2*C2
        dC2dt =  self.Gd2*O2 - (self.Ga2+self.Gr0)*C2
        #d5 = - (f0+f1+f2+f3+f4)
        #dr0 = -2*pi*f*sin(2*pi*f*t)*C(1+cos(2*pi*f*t))      # d/dt (A(1+cos(2*pi*f*t)))
        return np.array([dC1dt, dI1dt, dO1dt, dO2dt, dI2dt, dC2dt])

    def jacobian(self, s_0, t, phi_t=None):
        return np.array([[-self.Ga1, 0, self.Gd1, 0, 0, self.Gr0],
                         [self.Ga1, -self.Go1, 0, 0, 0, 0],
                         [0, self.Go1, -(self.Gd1+self.Gf), self.Gb, 0, 0],
                         [0, 0, self.Gf, -(self.Gb+self.Gd2), self.Go2, 0],
                         [0, 0, 0, 0, -self.Go2, self.Ga2],
                         [0, 0, 0, self.Gd2, 0, -(self.Ga2+self.Gr0)]])

    def hessian(self, s_0, t):
        """
        Hessian matrix for scipy.optimize.minimize.
        (Only for Newton-CG, dogleg, trust-ncg.)
        H(f)_ij(X) = D_iD_jf(X)
        """
        return np.zeros((6, 6))

    def calcSteadyState(self, phi):
        self.setLight(phi)
        Ga1 = self.Ga1
        Go1 = self.Go1
        Gf = self.Gf
        Gd2 = self.Gd2
        Gr0 = self.Gr0
        Gd1 = self.Gd1
        Gb = self.Gb
        Go2 = self.Go2
        Ga2 = self.Ga2
        denom = (Ga1*Go1*(Gf*(Go2*(Ga2+Gd2)+Gd2*Ga2)+Gb*Go2*Ga2)+Gd1*(Go1*Gb*Go2*Ga2+Ga1*Gb*Go2*Ga2+Gr0*(Go1*(Gb*Go2+Gd2*Go2)+Ga1*(Gb*Go2+Gd2*Go2)))+Gr0*(Ga1*(Go1*(Gb*Go2+Gd2*Go2+Gf*Go2)+Gf*Gd2*Go2)+Go1*Gf*Gd2*Go2))
        C1ss = (Gd1*(Go1*Gb*Go2*Ga2 + Go1*Gr0*(Gb*Go2 + Gd2*Go2)) + Go1*Gf*Gd2*Gr0*Go2)
        I1ss = (Gd1*(Ga1*Gb*Go2*Ga2 + Ga1*Gr0*(Gb*Go2 + Gd2*Go2)) + Ga1*Gf*Gd2*Gr0*Go2)
        O1ss = (Ga1*Go1*Gb*Go2*Ga2 + Ga1*Go1*Gr0*(Gb*Go2 + Gd2*Go2))
        O2ss = (Ga1*Go1*Gf*Go2*Ga2 + Ga1*Go1*Gf*Gr0*Go2)
        I2ss = (Ga1*Go1*Gf*Gd2*Ga2)
        C2ss = (Ga1*Go1*Gf*Gd2*Go2)
        self.steadyStates = np.array([C1ss, I1ss, O1ss, O2ss, I2ss, C2ss]) / denom
        return self.steadyStates

    def calcfphi(self, states=None):
        """Calculate the conductance scalar from the photocycle."""
        if states is None:
            states = self.states
        C1, I1, O1, O2, I2, C2 = states.T
        gam = self.gam
        return O1 + gam * O2

    def calcSoln(self, t, s0=None):
        raise NotImplementedError(self.nStates)

class RhO_6Kstates(RhodopsinModel):
    """Class definition for the 6K-state model"""

    # Class attributes
    nStates = 6
    useAnalyticSoln = False
    s_0 = np.array([1, 0, 0, 0, 0, 0])  # [s1_0=1, s2_0=0, s3_0=0, s4_0=0, s5_0=0, s6_0=0] # array not necessary
    phi_0 = 0.0                         # Default initial flux
    stateVars = ['C1', 'I1', 'O1', 'O2', 'I2', 'C2']  # stateVars[0] is the 'ground' state
    stateLabels = ['$C_1$', '$I_1$', '$O_1$', '$O_2$', '$I_2$', '$C_2$']
    photoFuncs = ['_calcGa1', '_calcGa2', '_calcGf']
    photoRates = ['Ga1', 'Ga2', 'Gf']
    photoLabels = ['$G_{a1}$', '$G_{a2}$', '$G_{f}$']
    constRates = ['Go1', 'Go2', 'Ga3', 'Gd1', 'Gd2', 'Gb']
    constLabels = ['$G_{o1}$', '$G_{o2}$', '$G_{a3}$', '$G_{d1}$', '$G_{d2}$', '$G_{b}$']

    paramsList = ['g0', 'gam', 'phi_m', 'k1', 'k2', 'k3', 'p',
                  'Gf0', 'k_f', 'Gb0', 'k_b', 'q', 'Go1', 'Go2',
                  'Ga3', 'Gd1', 'Gd2', 'Gb', 'E', 'v0', 'v1']  # List of model constants

    connect = [[0, 0, 1, 1, 0, 1],  # s_1 --> s_i=1...6
               [1, 0, 1, 0, 0, 0],  # s_2 -->
               [0, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 1],
               [1, 0, 0, 1, 0, 0]]

    equations = r"""
                $$ \dot{C_1} = G_{d1}O_1 + G_{b}C_2 + G_{a3}O_2 - (G_{a1} + G_{f})(\phi)C_1 $$
                $$ \dot{I_1} = G_{a1}(\phi)C_1 - G_{o1}I_1 $$
                $$ \dot{O_1} = G_{o1}I_1 - G_{d1}O_1 $$
                $$ \dot{O_2} = G_{o2}I_2 - (G_{d2} + G_{a3})O_2 $$
                $$ \dot{I_2} = G_{a2}(\phi)C_2 - G_{o2}I_2 $$
                $$ \dot{C_2} = G_{d2}O_2 + G_{f}(\phi)C1 - (G_{b} + G_{a2}(\phi))C_2 $$
                $$ C_1 + I_1 + O_1 + O_2 + I_2 + C_2 = 1 $$
                $$$$
                $$ G_{a1}(\phi) = k_{1} \frac{\phi^p}{\phi^p + \phi_m^p} $$
                $$ G_{a2}(\phi) = k_{2} \frac{\phi^p}{\phi^p + \phi_m^p} $$
                $$ G_{a3}(\phi) = k_{3} \frac{\phi^p}{\phi^p + \phi_m^p} $$
                $$ G_{f}(\phi)  = k_{f} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{f0} $$
                $$ G_{b}(\phi)  = k_{b} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{b0} $$
                $$$$
                $$ f_{\phi}(\phi) = O_1+\gamma O_2 $$
                $$ f_v(v) = v_1\frac{1-e^{-(v-E)/v_0}}{(v-E)} $$
                $$ I_{\phi} = g_0 \cdot f_{\phi}(\phi) \cdot f_v(v) \cdot (v-E) $$
                """

    brianStateVars = ['C_1', 'I_1', 'O_1', 'O_2', 'I_2', 'C_2']

    brian = '''
            dC_1/dt = Gd1*O_1 + Gb*C_2 + Ga3*O_2 - (Ga1+Gf)*C_1 : 1
            dI_1/dt = Ga1*C_1 - Go1*I_1                         : 1
            dO_1/dt = Go1*I_1 - Gd1*O_1                         : 1
            dO_2/dt = Go2*I_2 - (Gd2+Ga3)*O_2                   : 1
            dI_2/dt = Ga2*C_2 - Go2*I_2                         : 1
            C_2 = 1 - C_1 - I_1 - O_1 - O_2 - I_2               : 1
            H_p = Theta*((phi**p)/(phi**p+phi_m**p))            : 1
            H_q = Theta*((phi**q)/(phi**q+phi_m**q))            : 1
            Ga1 = k1*H_p                                        : second**-1
            Ga2 = k2*H_p                                        : second**-1
            Ga3 = k3*H_p                                        : second**-1
            Gf = k_f*H_q + Gf0                                  : second**-1
            Gb = k_b*H_q + Gb0                                  : second**-1
            f_v = (1-exp(-(v-E)/v0))/((v-E)/v1)                 : 1
            f_phi = O_1+gam*O_2                                 : 1
            I = g0*f_phi*f_v*(v-E)                              : amp
            phi                                                 : metre**-2*second**-1 (shared)
            Theta = int(phi > 0*phi)                            : 1 (shared)
            '''

    brian_phi_t = '''
            dC_1/dt = Gd1*O_1 + Gb*C_2 + Ga3*O_2 - (Ga1+Gf)*C_1 : 1
            dI_1/dt = Ga1*C_1 - Go1*I_1                         : 1
            dO_1/dt = Go1*I_1 - Gd1*O_1                         : 1
            dO_2/dt = Go2*I_2 - (Gd2+Ga3)*O_2                   : 1
            dI_2/dt = Ga2*C_2 - Go2*I_2                         : 1        
            C_2 = 1 - C_1 - I_1 - O_1 - O_2 - I_2               : 1
            Theta = int(phi(t) > 0*phi(t))                      : 1 (shared)
            H_p   = Theta*((phi(t)**p)/(phi(t)**p+phi_m**p))    : 1
            H_q   = Theta*((phi(t)**q)/(phi(t)**q+phi_m**q))    : 1
            Ga1   = k1*H_p                                      : second**-1
            Ga2   = k2*H_p                                      : second**-1
            Ga3   = k3*H_p                                      : second**-1
            Gf    = k_f*H_q + Gf0                               : second**-1
            Gb    = k_b*H_q + Gb0                               : second**-1
            f_v   = (1-exp(-(v-E)/v0))/((v-E)/v1)               : 1
            f_phi = O_1+gam*O_2                                 : 1
            I     = g0*f_phi*f_v*(v-E)                          : amp
            '''

    def _calcGa1(self, phi):
        #return self.a10*(phi/self.phi0)
        return self.k1 * phi**self.p/(phi**self.p + self.phi_m**self.p)

    def _calcGa2(self, phi):
        #return self.b40*(phi/self.phi0)
        return self.k2 * phi**self.p/(phi**self.p + self.phi_m**self.p)

    def _calcGf(self, phi):
        #return self.a30 + self.a31*np.log(1+(phi/self.phi0))
        return self.Gf0 + self.k_f * phi**self.q/(phi**self.q + self.phi_m**self.q)

    def setLight(self, phi):
        if phi < 0:
            phi = 0
        self.phi = phi
        self.Ga1 = self._calcGa1(phi)
        self.Ga2 = self._calcGa2(phi)
        self.Gf = self._calcGf(phi)
        if config.verbose > 1:
            self.dispRates()

    def dispRates(self): # This needs rechecking
        """Print photosensitive transition rates"""
        print("Transition rates (phi={:.3g}): C1 --[Ga1={:.3g}]--> O1 --[Gf={:.3g}]--> O2".format(self.phi, self.Ga1, self.Gf))
        #print("Transition rates (phi={:.3g}): O1 <--[Gb={:.3g}]-- O2 <--[Ga2={:.3g}]-- C2".format(self.phi, self.Gb, self.Ga2))

    def solveStates(self, s_0, t, phi_t=None):
        """Differential equations of the 6-state model to be solved by odeint"""
        if phi_t is not None:
            self.setLight(float(phi_t(t)))
        C1, I1, O1, O2, I2, C2 = s_0  # Unpack state vector
        dC1dt = -(self.Ga1+self.Gf)*C1 + self.Gd1*O1 + self.Ga3*O2 + self.Gb*C2
        dI1dt =  self.Ga1*C1 - self.Go1*I1
        dO1dt =  self.Go1*I1 - self.Gd1*O1
        dO2dt = -(self.Ga3+self.Gd2)*O2 + self.Go2*I2
        dI2dt = -self.Go2*I2 + self.Ga2*C2
        dC2dt = self.Gf*C1 + self.Gd2*O2 - (self.Ga2+self.Gb)*C2
        return np.array([dC1dt, dI1dt, dO1dt, dO2dt, dI2dt, dC2dt])

    def jacobian(self, s_0, t, phi_t=None):
        return np.array([[-(self.Ga1+self.Gf),  0,          self.Gd1,   self.Ga3,               0,          self.Gb],
                         [self.Ga1,             -self.Go1,  0,          0,                      0,          0],
                         [0,                    self.Go1,   -self.Gd1,  0,                      0,          0],
                         [0,                    0,          0,          -(self.Ga3+self.Gd2),   self.Go2,   0],
                         [0,                    0,          0,          0,                      -self.Go2,  self.Ga2],
                         [self.Gf,              0,          0,          self.Gd2,               0,          -(self.Ga2+self.Gb)]])

    def hessian(self, s_0, t):
        """
        Hessian matrix for scipy.optimize.minimize.
        (Only for Newton-CG, dogleg, trust-ncg.)
        H(f)_ij(X) = D_iD_jf(X)
        """
        return np.zeros((6, 6))

    def calcSteadyState(self, phi):
        self.setLight(phi)
        Ga1 = self.Ga1
        Go1 = self.Go1
        Gf = self.Gf
        Gd2 = self.Gd2
        Gd1 = self.Gd1
        Gb = self.Gb
        Go2 = self.Go2
        Ga2 = self.Ga2
        Ga3 = self.Ga3
        
        denom = (Ga1*Ga2*Ga3*Gd1*Go2 + Ga1*Ga2*Ga3*Go1*Go2 + Ga1*Ga3*Gb*Gd1*Go2 + Ga1*Ga3*Gb*Go1*Go2 + Ga1*Gb*Gd1*Gd2*Go2 + Ga1*Gb*Gd2*Go1*Go2 + Ga2*Ga3*Gd1*Gf*Go1 + Ga2*Ga3*Gd1*Go1*Go2 + Ga2*Gd1*Gd2*Gf*Go1 + Ga2*Gd1*Gf*Go1*Go2 + Ga3*Gb*Gd1*Go1*Go2 + Ga3*Gd1*Gf*Go1*Go2 + Gb*Gd1*Gd2*Go1*Go2 + Gd1*Gd2*Gf*Go1*Go2)
        C1ss = (Ga2*Ga3*Gd1*Go1*Go2 + Ga3*Gb*Gd1*Go1*Go2 + Gb*Gd1*Gd2*Go1*Go2)
        I1ss = (Ga1*Ga2*Ga3*Gd1*Go2 + Ga1*Ga3*Gb*Gd1*Go2 + Ga1*Gb*Gd1*Gd2*Go2)
        O1ss = (Ga1*Ga2*Ga3*Go1*Go2 + Ga1*Ga3*Gb*Go1*Go2 + Ga1*Gb*Gd2*Go1*Go2)
        O2ss = (Ga2*Gd1*Gf*Go1*Go2)
        I2ss = (Ga2*Ga3*Gd1*Gf*Go1 + Ga2*Gd1*Gd2*Gf*Go1)
        C2ss = (Ga3*Gd1*Gf*Go1*Go2 + Gd1*Gd2*Gf*Go1*Go2)
        self.steadyStates = np.array([C1ss, I1ss, O1ss, O2ss, I2ss, C2ss]) / denom
        return self.steadyStates

    def calcfphi(self, states=None):
        """Calculate the conductance scalar from the photocycle"""
        if states is None:
            states = self.states
        C1, I1, O1, O2, I2, C2 = states.T
        gam = self.gam
        return O1 + gam * O2

    def calcSoln(self, t, s0=None):
        raise NotImplementedError(self.nStates)


models = {
    '3': RhO_3states, 3: RhO_3states,
    '4': RhO_4states, 4: RhO_4states,
    '6': RhO_6states, 6: RhO_6states,
    '6K': RhO_6Kstates
}


def selectModel(nStates):
    """Model selection function."""
    if int(nStates) == 3 or nStates == 'three':
        return RhO_3states()
    elif int(nStates) == 4 or nStates == 'four':
        return RhO_4states()
    elif int(nStates) == 6 or nStates == 'six':
        return RhO_6states()
    elif nStates == '6K':
        return RhO_6Kstates()
    else:
        print("Error in selecting model - please choose from 3, 4 or 6 states")
        raise NotImplementedError(nStates)
