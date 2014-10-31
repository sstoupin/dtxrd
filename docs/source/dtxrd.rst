
.. _dtxrd:

************
dtxrd
************

x-ray diffraction calculator 
(dynamical theory of x-ray diffraction for perfect crystals)

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>

SYNOPSIS
============

::

       dtxrd [options] crystal h k l eta phi T d ["a" | "e"] [theta | Ex]


DESCRIPTION
============

A  program to calculate parameters of a Bragg or Laue reflection for 
a monochromatic incident wave using dynamical theory of x-ray diffraction for perfect crystals in the 
2-beam approximation

For a brief summary run::

    dtxrd -h

INPUT PARAMETERS
=================

:crystal:
       crystal type: C (diamond), Si (silicon), Ge (germanium) or Al2O3 (sapphire)

:h k l:  Miller indicies of a chosen reflection

:eta:    asymmetry angle (:math:`\eta` [degrees])

:phi:    azimuthal angle of incidence (:math:`\phi` [degrees])

:T:      crystal temperature [K]

:d:      crystal thickness [mm]

:flag: =====   =================================================================
       flag    description
       =====   =================================================================
       a       perform calculation at a given glancing angle of incidence theta
       e       perform calculation at a given photon energy Ex
       =====   =================================================================

:theta: glancing angle of incidence, theta (:math:`\theta`)

:Ex: photon energy, Ex (:math:`E_{\mathrm X}`)


OPTIONS
============

:-v, --version:
       show version of program.

:-h, --help:
       show summary of options.

:-o F, --output=F:
       write results to file F (default to stdout)

:-w D, --write=D:
       write data to file D (default - no action)

:-p, --pi:
       :math:`\pi` polarization for incident wave (default - :math:`\sigma` polarization)

:-c, --conv:
       convolve data with a virtual instrumental resolution function having FWHM of 1/10 of  the  Darwin  width
       and report the resulting FWHM of the reflectivity curve


OUTPUT PARAMETERS
======================

Basic parameters of the chosen h k l reflection: 

:d[A]:     :math:`d` [Angstrom] interplanar distance (d-spacing) of the chosen h k l reflection
       
:Eb[keV]:  :math:`E_B = \frac{hc}{2d}` [keV] Bragg energy

:thr[deg]: :math:`\theta_R` [degrees] incident glancing angle for the exact backscattering
	   (a wave with photon energy :math:`E_R` incident at this angle is reflected exactly backwards)

:Er[keV]:  :math:`E_R` [keV] photon energy for the exact backscattering

:bh:       :math:`b_{H}` asymmetry factor in the chosen scattering geometry 
           for symmetric reflection :math:`\eta = 0` and :math:`b_{H} = - 1`

Susceptibilities and refraction corrections:

:chi_{0}:  :math:`\chi_0` susceptibility 

:chi_{h}:  :math:`\chi_{H}` susceptibility 

:chi_{-h}: :math:`\chi_{\bar{H}}` susceptibility 

:wh(s):    :math:`\omega_{H}^s` refraction correction for symmetric reflection  

:wh:       :math:`\omega_{H} = \omega_{H}^s \frac{b_{H}-1}{2b_{H}}` refraction correction for the chosen reflectoin  

Central energy and angle:

:Ec[keV]:  :math:`E_c` [keV] central energy of the chosen reflection

:thc[deg]: :math:`\theta_c` [deg] central glancing angle of incidence of the chosen reflection 

Energy intrinsic (Darwin) widths (thick non-absorbing crystal) at fixed glancing angle of incidence :math:`\theta_c`:

:eps_s:   :math:`\varepsilon^s` relative energy width of symmetric h k l reflection (same for entrance and exit)
 
:eps:     :math:`\varepsilon` relative entrance energy width of the chosen h k l reflection  

:eps_pr:  :math:`\varepsilon'` relative exit energy width of the chosen h k l reflection 

:Delta_E_s[meV]:   :math:`\Delta E^s` [meV] absolute energy width of symmetric h k l reflection (same for entrance and exit)

:Delta_E[meV]:     :math:`\Delta E` [meV] absolute entrance energy width of the chosen h k l reflection 

:DeltaE_pr[meV]:   :math:`\Delta E'` [meV] absolute exit energy width of symmetric reflection 

Angular intrinsic (Darwin) widths (thick non-absorbing crystal) at fixed photon energy :math:`E_c`:

:dth_s[urad]:      :math:`\Delta \theta^s` [microradian] angular width of the symmetric h k l reflection  (same for entrance and exit)

:dth[urad]:        :math:`\Delta \theta` [microradian] angular entrance width of the chosen h k l reflection  

:dth_s[urad]:      :math:`\Delta \theta'` [microradian] angular exit width of the chosen h k l reflection 

Additional characteristics of the chosen h k l reflection:

:dE/dth[meV/urad]: :math:`\frac{dE}{d\theta}` [meV/microradian] tangent of the Bragg's Law

:Dr[urad/meV]:     :math:`D_r` [microradian/meV] intrinsic angular dispersion rate of the chosen h k l reflection 

:de[um]:           :math:`d_e` [micrometer] extinction length of the chosen h k l reflection

Reflectivity and Transmissivity:

:Rc[%]:            :math:`R_c` [%] reflectivity at center

:Tc[%]:            :math:`T_c` [%] transmissivity at center


EXAMPLES
===========

to calculate a rocking curve of a 1-mm-thick Si (111) crystal at 8 keV (111 reflection, Bragg case) run::

       dtxrd Si 1 1 1 0 0 300 1 e 8

.. image:: ../../examples/snapshots/Si111_8keV.png
            :width: 90 %
	    :alt: Si111 at 8keV

to calculate a rocking curve of a 0.1-mm-thick C (001) crystal at 12 keV (220 reflection, Laue case) run::

       dtxrd C 2 2 0 45 0 300 0.1 e 12 

.. image:: ../../examples/snapshots/C220_Laue.png
            :width: 90 %
	    :alt: C220 Laue at 12keV


SEE ALSO
============

* :ref:`throughput`
* :ref:`rcpeak`

:author: Stanislav Stoupin
:email:  <sstoupin@gmail.com>
:date: |today|
