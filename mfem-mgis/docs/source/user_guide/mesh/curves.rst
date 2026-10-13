===============
Defining curves
===============

This section describes the utility functions available in ``MFEM/MGIS``
for creating points and curves from parameters (header
``MFEMMGIS/Geometry.hxx``, which requires ``MGIS`` built with ``TFEL``), and
how to export values along them.

.. contents::
     :depth: 3
     :local:


Use Cases
^^^^^^^^^

The geometry utilities are particularly useful for:

- Defining boundary conditions or loads along specific paths in a mesh.
- Creating custom geometric features for post-processing or visualization.
- Generating points for interpolation or sampling along curves.

Overview
^^^^^^^^

The ``makePoint<N>`` function (``N`` being 2 or 3) creates a point from a
parameter holding its coordinates, for instance
``mfem_mgis::makePoint<2>(ctx, std::vector<mfem_mgis::Parameter>{0, 0.2})``.

The ``makePointsOnCurve`` function generates a set of points along a curve
based on a geometric description and a discretization strategy.

+--------------------------+-------------------------------------------------------------------+
| Function                 | Description                                                       |
+==========================+===================================================================+
| ``makePointsOnCurve<2>`` | Generates points along a curve in 2D space.                       |
+--------------------------+-------------------------------------------------------------------+
| ``makePointsOnCurve<3>`` | Generates points along a curve in 3D space.                       |
+--------------------------+-------------------------------------------------------------------+

**Parameters:**

The function accepts a ``Parameters`` object structured as follows:

- **Curve Type**: Currently, only the ``Line`` type is supported.

  - ``InitialPoint``: A vector of real numbers representing the starting point of the line.
  - ``FinalPoint``: A vector of real numbers representing the ending point of the line.
  - ``Discretization``: A nested ``Parameters`` object describing how the line is discretized.

**Discretization Strategies:**

1. **Uniform Discretization**

   - ``NumberOfPoints``: An integer specifying the number of points to generate along the line.
     Must be greater than or equal to 2.

   **Example:**

   .. code-block:: c++

      const auto params = mfem_mgis::Parameters{
          {"Line",
           mfem_mgis::Parameters{
               {"InitialPoint", std::vector<mfem_mgis::Parameter>{-2, 0.5}},
               {"FinalPoint", std::vector<mfem_mgis::Parameter>{6, 8}},
               {"Discretization",
                mfem_mgis::Parameters{
                    {"Uniform", mfem_mgis::Parameters{{"NumberOfPoints", 10}}}}}
           }}};
      const auto points = mfem_mgis::makePointsOnCurve<2>(ctx, params);

2. **Geometric Progression Discretization**

   - ``NumberOfPoints``: An integer specifying the number of points to generate along the line.
     Must be greater than or equal to 2.
   - ``NormalizedInitialDensity``: A positive real number representing the density of points at the start of the line. This density is normalised for a curve of length 1.
   - ``NormalizedFinalDensity``: A positive real number representing the density of points at the end of the line. This density is normalised for a curve of length 1.

   **Example:**

   .. code-block:: c++

      const auto params = mfem_mgis::Parameters{
          {"Line",
           mfem_mgis::Parameters{
               {"InitialPoint", std::vector<mfem_mgis::Parameter>{2, 0}},
               {"FinalPoint", std::vector<mfem_mgis::Parameter>{17, 0}},
               {"Discretization",
                mfem_mgis::Parameters{
                    {"GeometricProgression",
                     mfem_mgis::Parameters{
                         {"NumberOfPoints", 10},
                         {"NormalizedInitialDensity", mfem_mgis::real{1} / 150},
                         {"NormalizedFinalDensity", mfem_mgis::real{1} / 3}}}}}
           }}};
      const auto points = mfem_mgis::makePointsOnCurve<2>(ctx, params);

.. note::

   The geometric progression discretization generates points with varying
   densities along the line, allowing for finer discretization at specific
   regions of the curve.

Error Handling
^^^^^^^^^^^^^^

Both ``makePoint`` and ``makePointsOnCurve`` functions return an
``std::optional`` type. If the input parameters are invalid or inconsistent
(e.g., wrong number of coordinates, missing parameters, invalid types),
the functions return an empty optional, and an error message is registered
in the context.

Common errors include:

- Missing required parameters (e.g., ``InitialPoint``, ``FinalPoint``, ``Discretization``).
- Incorrect parameter types (e.g., passing a string where a number is expected).
- Inconsistent space dimensions (e.g., providing 3 coordinates for a 2D point).
- Invalid discretization parameters (e.g., ``NumberOfPoints`` less than 2).

**Example of Error Handling:**

.. code-block:: c++

   const auto invalid_params = mfem_mgis::Parameters{
       {"Line",
        mfem_mgis::Parameters{
            {"InitialPoint", std::vector<mfem_mgis::Parameter>{-2, 0.5, 2}},  // Invalid: 3D point for 2D curve
            {"FinalPoint", std::vector<mfem_mgis::Parameter>{6, 8}},
            {"Discretization",
             mfem_mgis::Parameters{{"Uniform", mfem_mgis::Parameters{{"NumberOfPoints", 10}}}}}
        }}};
   const auto points = mfem_mgis::makePointsOnCurve<2>(ctx, invalid_params);
   if (!points) {
       // Handle error: points is empty
   }

Exporting values along points sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Named points and points sets can be declared with the mesh, by the
``Points`` and ``PointsSets`` parameters of the finite element
discretization. The points of a line may then be given by their names.
The ``PointsSetCurvesWriter`` class (header
``MFEMMGIS/PostProcessing/PointsSetCurvesWriter.hxx``) writes the values of
grid functions at the points of a points set to a file. Its parameters are
``File`` (required), ``PointsSet`` (name or definition of the points set,
required), ``Precision``, ``ExportCurvilinearAbscissa`` and
``ExportCoordinates`` (``true`` by default).

**Example** (adapted from the test ``tests/PointsSetCurveTest.cxx``):

.. code-block:: c++

   using namespace mfem_mgis;
   const auto line = dict{
       {"Line", dict{{"InitialPoint", "C"},
                     {"FinalPoint", "B"},
                     {"Discretization",
                      dict{{"Uniform", dict{{"NumberOfPoints", 20}}}}}}}};
   // fed and writer are std::optional objects, empty on failure
   auto fed = construct<FiniteElementDiscretization>(
       ctx, dict{{"MeshFileName", "mesh.msh"},
                 {"FiniteElementFamily", "H1"},
                 {"FiniteElementOrder", 1},
                 {"UnknownsSize", 2},
                 {"Parallel", true},
                 {"Points", dict{{"B", list{1, 0.2}}, {"C", list{0, 0.2}}}},
                 {"PointsSets", dict{{"lCB", line}}}});
   auto writer = construct<PointsSetCurvesWriter>(
       ctx, *fed, dict{{"File", "values.txt"}, {"PointsSet", "lCB"}});
   // f is a grid function defined on the finite element space of fed
   writer->add(ctx, "f", f) | or_die;
   writer->writeFileHeader(ctx) | or_die;
   writer->writeValues(ctx, {.begin = 0, .end = 0, .dt = 0}, ets) | or_die;

The ``PointsSetCurves`` and ``Curves`` post-processings of the physical
systems (classes ``PointsSetCurvesPostProcessing`` and
``CurvesPostProcessing``, header directory ``MFEMMGIS/PostProcessing``)
write such values at each time step of a simulation (see
:ref:`simulation_post_processings`).
