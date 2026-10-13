============
Getting help
============

The :cxx:`info` function
========================

The :cxx:`info` function allows displaying information about an object
in an output stream.

Example of usage
----------------

.. code:: c++

   const auto& fed = problem.getFiniteElementDiscretization();
   const auto success = mfem_mgis::info(ctx, fed);

The :cxx:`info` function is declared in the header
:code:`MFEMMGIS/Info.hxx`. A first overload prints the information in a
given output stream, :cxx:`info(ctx, os, t)`, the second one, used above,
in the log stream of the execution context.

The information is given by specializations of the :cxx:`getInformation`
function, available for the following types:

- :cxx:`MeshDiscretization`,
- :cxx:`FiniteElementDiscretization`,
- :cxx:`PartialQuadratureSpace` and
  :cxx:`PartialQuadratureSpaceInformation`.

For the other types, :cxx:`getInformation` does nothing: :cxx:`info`
prints nothing.
