Fitting
=======

PyRhO incorporates novel fitting algorithms with which each of the opsin models may be parameterised when given an appropriate set of data. The fitting algorithms use the ``lmfit`` module (Newville, 2014) which in addition to providing access to a variety of optimisation algorithms, enables numerical bounds to be placed on parameter values as well as algebraic constraints. Parameters may also be manually fixed if, for example, they have been directly measured experimentally. Once the algorithm has finished, a ``Parameters`` object is returned, populated with the extracted parameter values which may then be saved or used as the basis for simulations. Plots of photocurrents simulated with these derived parameters are drawn over the experimental data (with residual error) to allow for visual comparison of the resultant fits. 


Characterisation data
---------------------

In order to characterise each model, a set of voltage-clamped photocurrents are required, ideally collected from HEK cells to eliminate the confounding effects of other ion channels which may be present in neurons. To capture all currently modelled variable dependencies, data from three stimulation protocols are necessary, listed below by the model properties which they reveal. 

In the event of scarce data or uncharacterised variables which the user does not intend to vary, we describe the minimum set of data for the fitting procedure below and discuss the implications for the resultant model. 

* **Voltage dependence**: :math:`{E, v_0, v_1}`. Long light pulse (fixed flux) to steady-state, vary voltage clamp potential (e.g. in steps of 30 mV: :math:`V_{\rm clamp} = \{-100, -70, -40, -10, 20, 50, 80\}\ \mathrm{mV}`, :math:`n \ge 5`). Voltage clamp values should not be too close to :math:`E` as this may cause distortions in the fitting algorithms. The software will automatically find the plateau values :math:`I_{\rm ss}`, plot :math:`I_{\rm ss}` vs. :math:`V_{\rm clamp}`, and find the fitting parameters for the function :math:`f_{\rm v}(v)` given by Equation~(\ref{eq:fv}). An example is shown in Figure~\ref{fig:rectifier}.

* **Recovery rate**: :math:`\{G_{r0} \}`. Two long light pulses with varying inter-pulse-interval (IPI), Voltage clamp: :math:`-70\ \mathrm{mV}`. Light on (first pulse) -- light off for e.g. :math:`t_{\rm IPI}=\{0.5, 1, 2.5, 5, 10\}\,` seconds -- light on (second pulse). The software will automatically find the peak values for each recording, align the data to the end of the first pulse and fit an appropriate exponential curve of the form :math:`I_{\rm peak}(t) = I_{peak0} - a\cdot e^{-G_{r0}\cdot t_{\rm IPI}}`. We note here that this expression is strictly speaking correct only when both :math:`O_1` and :math:`O_2` states are empty. Consequently :math:`a` is left as a free fitting parameter and very short values for :math:`t_{\rm IPI}` should be avoided to prevent the distortions caused by the faster transitions. An example for wild-type ChR2 is given in Figure~\ref{fig:recovery}, where :math:`t_{\rm IPI} \gtrsim 100 \mathrm{ms}`. 

* **Flux dependence**: `Off-curve`: :math:`\{ G_{d(1,2)}, [G_{f0}, G_{b0}]\}`; `On-curve`: :math:`\{$All other parameters$\}`. Voltage clamp (preferably): :math:`-70\ \mathrm{mV}`, long pulse to steady-state, (e.g. T :math:`\approx` 500ms) plus decay of off-curve. Vary light intensity from near threshold to saturation (e.g. :math:`\phi = \{0.1, 0.5, 1, 5, 10, 50, 100\}\ \mathrm{mW/mm}^2`, :math:`n \ge 5`). The recorded off- and on-curves are automatically fitted. An example set is shown in Figure~\ref{fig:fit6statesSet} with more details of the algorithm given in Appendix~\ref{sec:model-dep}.


Additionally the six-state model requires one or more very short pulses in order to characterise the opsin activation rates which model the lag in transitioning to conductive states upon light stimulation:

* **Opsin activation rate**: :math:`\{ G_{o1}, G_{o2}\}` One or more `short` pulses, voltage clamp: :math:`-70\ \mathrm{mV}`. Vary pulse length, e.g. 0.5ms, 1ms, 2ms, 3ms, varied up to 10ms. PyRhO will automatically find the time of the peak current and use an iterative formula to estimate :math:`G_{o1}`. We initially assume :math:`G_{o2}=G_{o1}`. Further details of the algorithm are given in Appendix~\ref{sec:Go}. 

All light pulses should be `rectangular' (step functions) in that they have a sharp onset and offset. Examples of each protocol are included in PyRhO with illustrations provided in Figures~\ref{fig:rectifier}--\ref{fig:fit6statesSet}. The duration of the on- and off-phases should also be kept approximately equal since the optimiser will effectively weight the contributions of each according to the relative numbers of data points. Additional parameter dependencies will be added in the future which may require additional data sets for a full characterisation of the models. 

Minimal data requirements
-------------------------

In general, the most important data are those described for characterising the flux dependence, which may be considered to be the `minimal set'. 
If this set consists of only a single photocurrent, the fitting algorithms will fix the parameters which model the flux dependence (:math:`\phi_m`, :math:`p` and :math:`q`) to the initial values supplied (along with fixing those describing other uncharacterised variables) and tune the remaining parameters to return a model fit for that specific flux. This is not recommended however, as the model is under-constrained by the data (typically resulting in a poorer fit than when using a whole set of photocurrents) and is unlikely to generalise well to new experimental conditions. For best results, the flux dependence photocurrents should be measured at light intensities spanning several orders of magnitude as described above. 

If variations in other parameters or short pulses are of interest then the additional data should (ideally) be collected as described. However, if obtaining the data for a full characterisation of the model is not possible, the pre-set default values should be adequate for most practical purposes. 


Fitting procedure and algorithms
--------------------------------

Once the data have been loaded into the appropriate structures, (see :doc:`expdata`) the fitting algorithms may be called with the ``fitModels()`` function. 

.. code-block::
    
    fp = fitModels(ChR2dataSet, nStates=6, params=initialParams)

This procedure returns a ``Parameters`` object (from the ``lmfit`` module) with the calculated values and plots the resultant model fits over the experimental photocurrents. 

In general terms, the fitting algorithm first finds the model-independent variables such as the dark recovery rate and voltage dependence factors, proceeding through `off-curve' parameters by fitting a double exponential decay function, optionally fitting opsin activation rates for the six-state model and finally optimising across a set of `on-curves' to find any remaining parameters. Due to the inherent variability and imprecision in experimental measurements there is an optional second optimisation phase over the entire set of photocurrents simultaneously. The values found for the dark parameters {:math:` G_{d(1,2)}` , :math:`[G_{f0}, G_{b0}]`} (and opsin activation rates :math:`G_{o(1,2)}` if relevant) are used as the initial values, lower and upper bounds are calculated as 50% and 200% of these values respectively (set by a hyperparameter) and the model is then re-optimised to achieve an overall better fit. The main sub-routines of the algorithm are given in Algorithm~\ref{alg:fitting} with more detail for each process given in the Appendix. 

\begin{algorithm}
\caption{Fitting algorithm
    \label{alg:fitting}}
  \begin{algorithmic}[1]
    \Function{fitModel}{$\{DataSet\}$, $\{initParams\}$}
    	\If{$``rectifier"$ in $\{DataSet\}$}
			\State $(E, v_0, v_1) \gets$ fitVoltageRectifier($\mathbf{V_{c}}, \mathbf{I_{ss}}$) %\\
        \EndIf
        \State $g_0 \gets$ fitConductance($v, E, \max_\phi(\mathbf{I_p})$) %\\ %\underset{1\le \phi \le N}{\max}
    	\If{$``recovery"$ in $\{DataSet\}$}
	        \State $G_{r0} \gets$ fitPeakRecovery($\mathbf{t_p}, \mathbf{I_p}$)
        \EndIf
        \State $(G_{d(1,2)}, [G_{f0}, G_{b0}]) \gets$ fitOffCurves($\{I_\phi[t_{off}:]\}$)
        \State $(\phi_m,  k_{(1,2)}, p, [G_{r1},G_{f0},k_{f},q,G_{b0},k_{b},\gamma,G_{o(1,2)}]) \gets$ fitOnCurves($\{I_\phi[t_{on}:t_{off}]\}$)
        \If{postFitOptimisation is \texttt{True}}
	        \State $(\{All\ parameters\}) \gets$ fitCurves($\{I_\phi\}$)
        \EndIf
    \EndFunction    
  \end{algorithmic}
\end{algorithm}

When fitting the three-state model, a double exponential is fit (with two corresponding decay rates :math:`G_{d1}` and :math:`G_{d2}`) which are then weighted by their coefficients (:math:`I_{\rm slow}` and :math:`I_{\rm fast}`) and combined to form a single exponential. The mean of these values is then calculated across a set of :math:`N` photocurrents (Equation~\ref{eq:Gd}) and this value is then used in subsequent parts of the fitting algorithm. 

.. math::

    G_d = \frac{1}{N}\sum_{n=1}^N{\frac{I_{\mathrm{slow}_n}\times G_{d1_n} + I_{\mathrm{fast}_n}\times G_{d2_n}}{I_{\mathrm{slow}_n}+I_{\mathrm{fast}_n}}}~. \label{eq:Gd}


Package: ``fitting`` API
========================

.. automodule:: pyrho.fitting
    :members: