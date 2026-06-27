Installation
============

``mbc_mom`` requires a compiled Rust backend to achieve its high performance. This guide covers how to set up your development environment, build the Rust extension, and install the Python package locally.

Prerequisites
-------------

Since the core engine is written in Rust and wrapped in Python, you will need both toolchains. This guide assumes you are developing on a Linux environment (such as WSL Ubuntu 24).

1. **Rust Toolchain**: You need ``rustc`` and ``cargo``. The easiest way to install this is via `rustup <https://rustup.rs/>`_:

   .. code-block:: bash

       curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

2. **Pixi**: We use ``pixi`` for fast, reproducible Python environment management. 

   .. code-block:: bash

       curl -fsSL https://pixi.sh/install.sh | bash


Setting up the Environment
--------------------------

1. **Clone the Repository**

   Clone the source code to your local machine and navigate into the project directory:

   .. code-block:: bash

       git clone git@github.com:aop007/mbc_mom.git
       cd mbc_mom

2. **Initialize the Pixi Environment**

   ``mbc_mom`` is targeted for Python 3.14. Your repository should contain a ``pixi.toml`` file that manages dependencies like ``maturin`` (the build system for PyO3) and ``numpy``. 
   
   If you haven't already initialized your environment, you can install the dependencies by running:

   .. code-block:: bash

       pixi install

Building the Rust Engine
------------------------

We use `Maturin <https://www.maturin.rs/>`_ to compile the Rust backend (utilizing ``rayon`` for parallelization) and build the Python package.

To build the extension in release mode and install it directly into your active ``pixi`` environment, run:

.. code-block:: bash

    pixi run maturin develop --release

.. note::
   Omitting the ``--release`` flag will build the package in debug mode. While compilation is much faster, the resulting solver will execute significantly slower due to unoptimized physics loops. Always use ``--release`` when running realistic electromagnetic simulations.

Verifying the Installation
--------------------------

Once compiled, you can verify that the Python frontend is successfully communicating with the Rust backend by launching a Python shell within the ``pixi`` environment:

.. code-block:: bash

    pixi run python -c "import mbc_mom.geometry; print('mbc_mom installed successfully!')"

If no errors are raised, your high-performance MoM solver is ready to use. Proceed to the :doc:`quickstart` guide to run your first simulation.