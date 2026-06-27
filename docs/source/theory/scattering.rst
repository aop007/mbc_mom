Scattering & Incident Fields
============================

While many MoM applications are driven by localized voltage gaps (antennas), ``mbc_mom`` is fully equipped for Radar Cross Section (RCS) and electromagnetic scattering analysis.

The Incident Voltage Vector
---------------------------

In a scattering problem, the excitation vector :math:`[V]` is populated by integrating an external, incident field :math:`\vec{E}^{inc}` along the path of each basis function:

.. math::

   V_i = \int_{L_i} \vec{E}^{inc}(\vec{r}) \cdot \vec{f}_i(\vec{r}) dl

Using ``get_incident_eval_points()``, the solver returns the exact 3D Cartesian coordinates where the incident field must be sampled.

Coupling External Fields
------------------------

Once you have evaluated your arbitrary incident field (e.g., a plane wave, a Gaussian beam, or near-field coupling from an adjacent system) at the returned spatial coordinates, you pass the complex :math:`E_x`, :math:`E_y`, and :math:`E_z` arrays to ``compute_incident_v_matrix()``. 

The Rust backend performs an analytical integration against the underlying MBC basis functions to yield the exact :math:`[V]` matrix.

.. note::
   The solver defaults to 7-point Gaussian quadrature (``points_per_seg=7``) for coupling incident fields to the mesh segments. This ensures highly accurate tracking of rapidly varying spatial fields across large segments.