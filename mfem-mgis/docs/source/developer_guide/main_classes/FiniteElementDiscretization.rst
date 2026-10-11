.. _mfem_mgis_FiniteElementDiscretization:

===============================
FiniteElementDiscretization
===============================

The ``FiniteElementDiscretization`` class extends ``MeshDiscretization`` to handle the lifetime of finite element collections and spaces. It provides a high-level interface for creating and managing finite element spaces on meshes or submeshes.

Main Features
-------------

- **Finite Element Space Management**: Creates and manages finite element spaces for both sequential and parallel computations.
- **Finite Element Collection Handling**: Stores and provides access to finite element collections.
- **Submesh Support**: Allows creation of finite element spaces on submeshes defined by material or boundary identifiers.
- **Reusability**: Reuses existing finite element spaces when possible to optimize memory and performance.
- **Integration with MeshDiscretization**: Inherits all mesh management capabilities from ``MeshDiscretization``.

Key Methods
-----------

- ``getFiniteElementSpace<parallel>()``: Returns the finite element space (parallel or sequential).
- ``getFiniteElementCollection()``: Returns the underlying finite element collection.
- ``getFiniteElementSpacesManager()``: Returns the manager for finite element spaces.
- ``getVSize()`` / ``getTrueVSize()``: Returns the total number of unknowns.

Usage Example
-------------

The following example demonstrates basic usage of ``FiniteElementDiscretization``:

.. code-block:: c++

  using namespace mfem_mgis;
  auto ctx = mgis::Context{};
  auto ofed = construct<FiniteElementDiscretization>(
      ctx, dict{{"MeshFileName", "mesh.mesh"},
                {"FiniteElementFamily", "H1"},
                {"FiniteElementOrder", 2},
                {"UnknownsSize", 2},
                {"NumberOfUniformRefinements", 0},
                {"Materials", dict{{"Attr1", 1}, {"Attr2", 2}}},
                {"Boundaries", dict{{"left", 1}, {"right", 2}}},
                {"Parallel", false}});
  
  // Get the finite element space
  auto& fes = ofed->getFiniteElementSpace<false>();
  
  // Get the finite element spaces manager
  auto fespaces = ofed->getFiniteElementSpacesManager();
  
  // Create a finite element space on a submesh
  auto fes1 = fespaces.getFiniteElementSpace<false>(
      ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
            .identifiers = list{"Attr1"},
            .number_of_components = 3});

Parameters
----------

In addition to the parameters supported by ``MeshDiscretization``, ``FiniteElementDiscretization`` supports:

- ``FiniteElementFamily`` (string): Name of the finite element family (e.g., ``H1``).
- ``FiniteElementOrder`` (int): Order of the polynomial approximation.
- ``FiniteElementSpaceOrdering`` (string): Ordering of the degrees of freedom, ``byNODES`` (all the values of the first component, then all the values of the second component, etc., default) or ``byVDIM`` (all the components of the first node, then all the components of the second node, etc.). The ``Elasticity`` strategy of the ``HypreBoomerAMG`` preconditioner requires ``byVDIM``.
- ``UnknownsSize`` (int): Number of components of the unknowns.
- ``IntegrationRuleOrder`` (int, optional): Order of the integration rules of the behaviour integrators. By default, each behaviour integrator uses its own rule, of order ``2 * t.OrderGrad(&e)``, i.e. 27 points on the hexahedra of order 1; the order 3 gives 8 points, which integrate exactly the elastic stiffness of affine hexahedra of order 1 and divide by about 3.4 the memory and the work at the integration points. The quadrature spaces of the behaviour integrators and their integrations use the same rule. The order is the same for all the geometries and all the orders of the elements: the order 3 under-integrates the hexahedra of order 2, and the orders 0 and 1 give one point on the hexahedra (hourglass modes). A rule having a weight which is not positive (rules of MFEM for some simplices) is refused.

.. note::

  The ``FiniteElementDiscretization`` class is designed to be lightweight and integrates seamlessly with the ``FiniteElementSpacesManager`` for advanced finite element space management.
