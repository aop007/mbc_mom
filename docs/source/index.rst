====================================
mbc_mom: High-Performance MoM Solver
====================================

``mbc_mom`` is a high-performance Method of Moments (MoM) computational electromagnetics solver. It pairs a highly parallelized Rust backend (powered by ``rayon``) with an intuitive Python interface. 

The engine features Multiradius Bridge-Currents (MBC) for seamless junction handling, analytical near-field integration, exact Sommerfeld-Norton surface wave formulations, and incident field scattering capabilities.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Theoretical Background

   theory/method_of_moments
   theory/ground_models
   theory/scattering

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/geometry
   api/solvers

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`