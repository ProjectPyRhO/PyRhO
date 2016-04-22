.. _models:

Rhodopsin models
================

Photocurrent model
------------------

We assume that all light-sensitive ion channel currents (:math:`I`) can
be expressed in the classic form:

.. math:: I = g \cdot (v-E)~,

where :math:`g` is the channel conductance, :math:`v` the membrane
voltage and :math:`E` reversal potential for the specific opsin type.
Generally speaking the ionic conductance is a complex function of light
flux (:math:`\phi(t)`), wavelength (:math:`\lambda`) and the opsin's
photocycle, membrane voltage, temperature (:math:`T`), and intracellular
and extracellular pH. We use a simplified empirical form for the channel
conductance, introduced by Hodgkin and Huxley, expressing it as a
product of a constant (:math:`g_0`, in our case this is maximum
conductance at :math:`v=-70\,\mathrm{mV}`), and a numerical coefficient
(:math:`f>0`):

.. math:: g = g(\phi,\lambda,v,T,pH,t) = g_0 \cdot f(\phi,\lambda,v,T,pH,t) ~,

In this version of PyRhO we have implemented the photocycle and membrane
voltage dependencies and assumed that these two contributions can be
separated:

.. math:: g = g_0 \cdot f_{\rm \phi}(\phi,t)\cdot f_{\rm v}(v) ~.

These two dependences are considered to be most relevant for
physiological electrolyte conditions, when temperature and pH are
considered to be fixed. Other dependencies will be implemented in the
next version of PyRhO.

Photocycle models
-----------------

At the core of PyRhO are three functional Markov models of opsin
kinetics, namely the three, four and six-state models. These models vary
in complexity providing a range in the trade-off between biological
accuracy and computational efficiency to choose from. The key features
of these models, including an outline of their strengths and weaknesses,
are summarised in the table below.

Both four and six state models assume that there are two open states
(:math:`O_1` and :math:`O_2`), with channel conductances
:math:`g_{\rm O_1}` and :math:`g_{\rm O_2}`, respectively. The
Photocycle factor has the form:

.. math:: f_{\rm \phi}(\phi) = O_1+\gamma O_2 ~,

where :math:`O_1` and :math:`O_2` are the fractions of opsins in two
open states in the interval :math:`[0,1]`, and
:math:`\gamma = g_{\rm O_2}/g_{\rm O_1}`. In contrast, the three-state
model assumes only one open state (:math:`O`) making the photocycle
factor simply: :math:`f_{\rm \phi}(\phi) = O`.

Voltage dependence
------------------

Here we assume that the membrane voltage affects only the ion-channel
conductance but not the channel kinetics. By investigating experimental
results for Channelrhodopsin-2 steady state current versus clamped
voltage the I-V curve shows inwardly rectifying behaviour. We have
previously found that an exponential function gives a good fit for this
dependency, therefore for the voltage factor we adopt the form:

.. math:: f_{\rm v}(v) = \frac{v_1}{v-E}\cdot \left(1-e^{-\frac{v-E}{v_0}}\right) ~,

where :math:`v_0` and :math:`v_1` are fitting parameters, along with
:math:`E`, the channel's reversal potential. The exponential dependence
on :math:`v` transforms into a linear dependence for large values of
:math:`v_0` which cause the exponent to be small and the expression
reduces to :math:`f_{\rm v}(v) \approx v_1/v_0 = const`, i.e. no direct
dependence on membrane voltage, which may be a more appropriate form for
some opsins. The equation for :math:`f_v` therefore generalises to both
cases for appropriate choices of the parameters :math:`v_0` and
:math:`v_1`.

Furthermore, since the voltage dependence factor is defined to be equal
to 1 at :math:`-70\ \mathrm{mV}`
(:math:`f_{\rm v}(-70\ \mathrm{mV}) := 1`), the value of :math:`v_0` is
related to the other parameters through the following equation:

.. raw:: latex

   \begin{equation}
   v_1 = \frac{70+E}{e^{\frac{70+E}{v_0}}-1}
   \end{equation}

Substituting in the expression for :math:`v_1`, the full expression for
:math:`f_v(v\vert E, v_0)` is therefore:

.. math:: f_v(v) = \frac{\left(E + 70\right) \left(e^{\frac{E-v}{v_{0}} } - 1\right)}{\left(E - v\right) \left(e^{\frac{E+70}{v_{0}} } - 1\right)}

.. math:: f_v(v) = \frac{\left(-70 - E\right) \left(1 - e^{-\frac{v-E}{v_{0}} }\right)}{\left(v - E\right) \left(1 - e^{-\frac{-70-E}{v_{0}} }\right)}

Summary of opsin models
~~~~~~~~~~~~~~~~~~~~~~~

+----------+---------------+--------------+------------------------------------+--------------------------------------+
| States   | Transitions   | Parameters   | Pros                               | Cons                                 |
+==========+===============+==============+====================================+======================================+
| 3        | 3             | 11           | Efficient analytic solution        | Single exponential off-phase decay   |
+----------+---------------+--------------+------------------------------------+--------------------------------------+
| 4        | 7             | 17           | Balance of detail and efficiency   | Lacks short-pulse dynamics           |
+----------+---------------+--------------+------------------------------------+--------------------------------------+
| 6        | 9             | 19           | Most detailed dynamics             | Computationally expensive            |
+----------+---------------+--------------+------------------------------------+--------------------------------------+

Three-state model
-----------------

Model description
~~~~~~~~~~~~~~~~~

The model was originally described in Nikolic et al. (2009) but has
since been extended.

Konstantin Nikolic, Nir Grossman, Matthew S. Grubb, Juan Burrone, Chris
Toumazou & Patrick Degenaar (2009). Photocycles of Channelrhodopsin-2.
*Photochemistry and Photobiology*, 85(1), 400-411.
`doi:10.1111/j.1751-1097.2008.00460.x <http://dx.doi.org/10.1111/j.1751-1097.2008.00460.x>`_

The three-state model shown in figure (a) consists of a closed state,
:math:`C`, an open state, :math:`O` and a deactivated state, :math:`D`.
The :math:`C\rightarrow O` transition is activated by light, starting
the photocycle and producing a photocurrent since the open state is
conductive. Additionally, the recovery process (:math:`D \rightarrow C`)
is assissted by light.

The system of differential equations describing the Markov state
variables are given below, where all state occupancies sum to 1.

.. raw:: latex

   \begin{align*}
   \frac{dC}{dt} &= G_{r}(\phi)D - G_{a}(\phi)C \\
   \frac{dO}{dt} &= G_{a}(\phi)C - G_{d}O \\
   \frac{dD}{dt} &= G_{d}O - G_{r}(\phi)D
   \end{align*}

.. math::  C + O + D = 1 

The light dependent transition rates for this system of ODEs are given
by the following equations.

.. raw:: latex

   \begin{align*}
   G_a(\phi) &= k_a\frac{\phi^p}{\phi^p + \phi_m^p} \\
   G_r(\phi) &= k_r\frac{\phi^q}{\phi^q + \phi_m^q} + G_{r0}
   \end{align*}

Analytic Solution
~~~~~~~~~~~~~~~~~

The analytic solution for arbitrary initial conditions
(:math:`C_0, O_0, D_0`) is as follows:

.. raw:: latex

   \begin{align}
   C(t) =& ~\frac{G_d G_r}{\lambda_1 \lambda_2} + \frac{(\lambda_1-G_d-G_r)}{G_d \lambda_1 \xi} \cdot \zeta_1 \cdot e^{-\lambda_1\cdot t} &&- \frac{(\lambda_2-G_d-G_r)}{G_d \lambda_2 \xi}\cdot \zeta_2 \cdot e^{-\lambda_2\cdot t} \\
   O(t) =& ~\frac{G_a G_r}{\lambda_1 \lambda_2} - \frac{(\lambda_1-G_r)}{G_d \lambda_1 \xi}\cdot \zeta_1 \cdot e^{-\lambda_1\cdot t} &&- \frac{(\lambda_2-G_r)}{G_d \lambda_2 \xi}\cdot \zeta_2 \cdot e^{-\lambda_2\cdot t} \\
   D(t) =& ~\frac{G_a G_d}{\lambda_1 \lambda_2} + \frac{1}{\lambda_1 \xi} \cdot \zeta_1 \cdot e^{-\lambda_1\cdot t} &&- \frac{1}{\lambda_2 \xi} \cdot \zeta_2 \cdot e^{-\lambda_2\cdot t} 
   \end{align}

with the following definitions:

.. math:: \xi = \sqrt{G_a^2 + G_d^2 + G_r^2 - 2\cdot(G_aG_d + G_aG_r + G_dG_r) }

.. math::  \lambda_{1,2} = \frac{1}{2}\cdot\left((G_a+G_d+G_r) \pm \xi\right) 

.. math::  \zeta_{1,2} = C_0\cdot \left[G_d\cdot G_a\right] + O_0\cdot \left[G_d \cdot (G_a - \lambda_{1,2})\right] + D_0\cdot \left[G_r \cdot (\lambda_{1,2} - G_a - G_d)\right] 

Note also that $ \ *1 *\ 2 = G\_a G\_d + G\_a G\_r + G\_d G\_r $ i.e.
the sum of products.

Steady-State Dynamics
~~~~~~~~~~~~~~~~~~~~~

An explicit solution for the model's dynamic variables at steady-state
can also be determined:

.. math::


   \begin{align}
   C_{SS} &= \frac{G_d G_r}{G_a G_d + G_a G_r + G_d G_r} \\
   O_{SS} &= \frac{G_a G_r}{G_a G_d + G_a G_r + G_d G_r} = \frac{\tau_d}{\tau_d + \tau_r + \tau_a} \\
   D_{SS} &= \frac{G_a G_d}{G_a G_d + G_a G_r + G_d G_r}
   \end{align}

Steady-State Current
~~~~~~~~~~~~~~~~~~~~

The steady-state (plateau) current, :math:`I_{ss}`, is therefore given by:

.. math:: I_{SS} = g_0\cdot\frac{\tau_d}{\tau_d + \tau_r + \tau_a}\cdot f_v(v) \cdot (v-E)

Remarks
~~~~~~~

The system still converges to a steady state, but oscillates as it does
if:
:math:`G_r + G_d - 2 \sqrt{G_r G_d} \lt G_a \lt G_r + G_d + 2 \sqrt{G_r G_d}`

Four-state model
----------------

The model was originally described in `Nikolic et al. (2009) <http://dx.doi.org/10.1111/j.1751-1097.2008.00460.x>`_ but has
since been extended.


The system of differential equations describing the Markov state
variables are given below.

.. math::


   \begin{align*} 
   \dot{C_1} &= G_{d1}O_1 + G_{r0}C_2 - G_{a1}(\phi)C_1 \\
   \dot{O_1} &= G_{a1}(\phi)C_1 + G_{b}(\phi)O_2 - \left(G_{d1}+G_{f}(\phi)\right)O_1 \\
   \dot{O_2} &= G_{a2}(\phi)C_2 + G_{f}(\phi)O_1 - \left(G_{d2}+G_{b}(\phi)\right)O_2 \\
   \dot{C_2} &= G_{d2}O_2 - \left(G_{r0}+G_{a2}(\phi)\right)C_2
   \end{align*}

.. math::  C_1 + O_1 + O_2 + C_2 = 1 

The light dependent transition rates for this system of ODEs are given
by the following equations.

.. math::


   \begin{align*}
   G_{a1}(\phi) &= k_1\frac{\phi^p}{\phi^p + \phi_m^p} \\
   G_{f}(\phi) &= k_{f} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{f0} \\
   G_{b}(\phi) &= k_{b} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{b0} \\
   G_{a2}(\phi) &= k_2\frac{\phi^p}{\phi^p + \phi_m^p}
   \end{align*}

Steady-State Dynamics
~~~~~~~~~~~~~~~~~~~~~

.. math::  C_{1ss} = \frac{G_{d1} \cdot (G_b \cdot (G_{a2} + G_{r0}) + G_{d2} G_{r0}) + G_f G_{d2} G_{r0}} {G_{a1} \cdot \left( (G_f + G_b) \cdot (G_{a2} + G_{r0}) + G_{d2} \cdot (G_f + G_{r0}) \right) + G_{d1} \cdot (G_{a2} G_b + G_{d2} G_{r0} + G_{r0} G_b) + G_{d2} G_{r0} G_f} 

.. math::  O_{1ss} = \frac{G_{a1} \cdot (G_b \cdot (G_{a2} + G_{r0}) + G_{d2} G_{r0})} {G_{a1} \cdot \left( (G_f + G_b) \cdot (G_{a2} + G_{r0}) + G_{d2} \cdot (G_f + G_{r0}) \right) + G_{d1} \cdot (G_{a2} G_b + G_{d2} G_{r0} + G_{r0} G_b) + G_{d2} G_{r0} G_f} 

.. math::  O_{2ss} = \frac{G_{a1} G_f \cdot (G_{a2} + G_{r0})} {G_{a1} \cdot \left( (G_f + G_b) \cdot (G_{a2} + G_{r0}) + G_{d2} \cdot (G_f + G_{r0}) \right) + G_{d1} \cdot (G_{a2} G_b + G_{d2} G_{r0} + G_{r0} G_b) + G_{d2} G_{r0} G_f} 

.. math::  C_{2ss} = \frac{G_{a1} G_f G_{d2}} {G_{a1} \cdot \left( (G_f + G_b) \cdot (G_{a2} + G_{r0}) + G_{d2} \cdot (G_f + G_{r0}) \right) + G_{d1} \cdot (G_{a2} G_b + G_{d2} G_{r0} + G_{r0} G_b) + G_{d2} G_{r0} G_f} 

Steady-State Current
~~~~~~~~~~~~~~~~~~~~

.. math::  I_{ss} = g_0 \cdot \frac{G_{a1} \cdot \left(G_{d2} G_{r0} + (G_b + \gamma G_f) \cdot (G_{a2} + G_{r0} )\right)} {G_{a1} \cdot \left( (G_f + G_b) \cdot (G_{a2} + G_{r0}) + G_{d2} \cdot (G_f + G_{r0}) \right) + G_{d1} \cdot (G_{a2} G_b + G_{d2} G_{r0} + G_{r0} G_b) + G_{d2} G_{r0} G_f} \cdot f_v(v) \cdot (v-E) 




Six-state model
---------------

The model was originally described in `Grossman et al. (2013) <http://dx.doi.org/10.1007/s10827-012-0431-7>`_ but has
since been extended. 

Grossman et al. (2013). The spatial pattern of light determines the
kinetics and modulates backpropagation of optogenetic action potentials.
Journal of Computational Neuroscience, 34(3), 477-488.
http://dx.doi.org/10.1007/s10827-012-0431-7

.. math::


   \begin{align*}
   \dot{C_1} &= G_{d1}O_1 + G_{r0}C_2 - G_{a1}(\phi)C_1 \\
   \dot{I_1} &= G_{a1}(\phi)C_1 - G_{o1}I_1 \\
   \dot{O_1} &= G_{o1}I_1 + G_{b}(\phi)O_2 - \left(G_{d1} + G_{f}(\phi)\right)O_1 \\
   \dot{O_2} &= G_{o2}I_2 + G_{f}(\phi)O_1 - \left(G_{d2} + G_{b}(\phi)\right)O_2 \\
   \dot{I_2} &= G_{a2}(\phi)C_2 - G_{o2}I_2 \\
   \dot{C_2} &= G_{d2}O_2 - \left(G_{r0}+G_{a2}(\phi)\right)C_2
   \end{align*}

.. math::  C_1 + I_1 + O_1 + O_2 + I_2 + C_2 = 1 

The light dependent transition rates for this system of ODEs are given
by the following equations.

.. math::


   \begin{align*}
   G_{a1}(\phi) &= k_1\frac{\phi^p}{\phi^p + \phi_m^p} \\
   G_{f}(\phi) &= k_{f} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{f0}\\
   G_{b}(\phi) &= k_{b} \frac{\phi^q}{\phi^q + \phi_m^q} + G_{b0}\\
   G_{a2}(\phi) &= k_2\frac{\phi^p}{\phi^p + \phi_m^p}
   \end{align*}

Stead-State dynamics
~~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align*}
   C_{1ss} &= \frac{G_{o1} G_{o2} (G_{a2} G_b G_{d1} + G_b G_{d1} G_{r0} + G_{d1} G_{d2} G_{r0} + G_f G_{d2} G_{r0})} {N} \\
   I_{1ss} &= \frac{G_{a1} G_{o2} (G_{a2} G_b G_{d1} + G_b G_{d1} G_{r0} + G_{d1} G_{d2} G_{r0} + G_f G_{d2} G_{r0})} {N} \\
   O_{1ss} &= \frac{G_{a1} G_{o1} G_{o2} (G_b G_{a2} + G_{r0} G_b + G_{r0} G_{d2})}{N} \\
   O_{2ss} &= \frac{G_{a1} G_{o1} G_f G_{o2}(G_{a2} + G_{r0})}{N} \\
   I_{2ss} &= \frac{G_{a1} G_{o1} G_f G_{d2} G_{a2}}{N} \\
   C_{2ss} &= \frac{G_{a1} G_{o1} G_f G_{d2} G_{o2}}{N}
   \end{align*}

Where the normalisation constant is as follows:

.. math::


   \begin{align*}
   N &= G_{a1} G_{o1} (G_f(G_{o2} (G_{a2}+G_{d2})+G_{d2} G_{a2}) + G_b G_{o2} G_{a2}) \\ &+ G_{d1}(G_{o1} G_b G_{o2} G_{a2}+G_{a1} G_b G_{o2} G_{a2}+G_{r0}(G_{o1}(G_b G_{o2}+G_{d2} G_{o2})+G_{a1}(G_b G_{o2}+G_{d2} G_{o2}))) \\ &+ G_{r0}(G_{a1}(G_{o1}(G_b G_{o2} + G_{d2} G_{o2} + G_f G_{o2}) + G_f G_{d2} G_{o2}) + G_{o1} G_f G_{d2} G_{o2})
   \end{align*}

Steady-State current
~~~~~~~~~~~~~~~~~~~~

.. math:: I_{ss} = g_0 \cdot \frac{G_{a1} G_{o1} G_{o2} \cdot \left(G_{d2} G_{r0} + (G_b + \gamma G_f) \cdot (G_{a2} + G_{r0})\right)}{N} \cdot f_v(v) \cdot(v-E)



Package: ``models`` API
=======================

.. automodule:: pyrho.models
	:members: