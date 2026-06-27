Ground Models & Boundary Conditions
===================================

``mbc_mom`` supports modeling structures in free space, over Perfect Electric Conductors (PEC), or over lossy, real ground using exact formulations. The presence of a ground plane modifies the free-space Green's function :math:`G_0(\vec{r}, \vec{r}')`.

PEC Ground
----------

When ``Mesh.set_pec_ground()`` is invoked, the solver uses classical Image Theory. The Green's function becomes:

.. math::

   G(\vec{r}, \vec{r}') = \frac{e^{-jkR_1}}{R_1} \pm \frac{e^{-jkR_2}}{R_2}

Where :math:`R_1`` is the distance to the source, and :math:`R_2`` is the distance to its virtual image.

Reflection Coefficient Approximation (RCA)
------------------------------------------

For real ground, RCA provides a fast approximation by scaling the image contribution using the Fresnel reflection coefficients (:math:`\Gamma_{TE}` and :math:`\Gamma_{TM}`).

.. math::

   G(\vec{r}, \vec{r}') \approx \frac{e^{-jkR_1}}{R_1} + \Gamma(\theta) \frac{e^{-jkR_2}}{R_2}

This is computationally cheap and highly accurate for antennas elevated at least $0.2\lambda$ above the ground, but fails for buried wires or surface-hugging structures.

Exact Sommerfeld-Norton Surface Waves
-------------------------------------

For geometries close to or penetrating a lossy half-space, ``mbc_mom`` relies on exact Sommerfeld integral evaluation. Activated via ``Mesh.set_real_ground(sigma, eps_r, use_sommerfeld=True)``, the engine evaluates the rigorous Sommerfeld identity:

.. math::

   G_{Sommerfeld} = 2 \int_{0}^{\infty} J_0(\lambda \rho) \frac{e^{-u z}}{u + u_E} \lambda d\lambda

Because this highly oscillatory integration is computationally paralyzing, the Rust backend pre-computes these exact integrals into a 2D interpolative look-up table during the initialization phase. The parallelization via ``rayon`` ensures this table is generated in milliseconds, providing rigorous accuracy including Norton surface waves without the traditional MoM bottleneck.