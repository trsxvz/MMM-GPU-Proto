/*!
 * \file   include/MFront/GenericMaterialProperty/Types.h
 * \brief
 * \author Thomas Helfer
 * \date 03/05/2022
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with linking exception
 * or the CECILL-A licence. A copy of thoses licences are delivered
 * with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICMATERIALPROPERTY_TYPES_H
#define LIB_MFRONT_GENERICMATERIALPROPERTY_TYPES_H

#ifdef __cplusplus
#include <cstddef>
#else  /* __cplusplus */
#include <stddef.h>
#endif /* __cplusplus */

// NOLINTBEGIN(modernize-use-using)
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*!
 * \brief numerical type used by the interface: the type given by the
 * `MFRONT_GMP_REAL_TYPE` macro if it is defined (`float`, for instance), or
 * `double`. `mfront` defines this macro at the beginning of the generated
 * sources of the `generic` and `generic-parallel` interfaces when the
 * `real_type` DSL option is given.
 */
#ifdef MFRONT_GMP_REAL_TYPE
typedef MFRONT_GMP_REAL_TYPE mfront_gmp_real;
#else  /* MFRONT_GMP_REAL_TYPE */
typedef double mfront_gmp_real;
#endif /* MFRONT_GMP_REAL_TYPE */
/*!
 * \brief offset of the versions of the application programming interfaces
 * exported by a library whose numerical type is not `double` (see
 * `MFRONT_GB_REAL_TYPE_API_VERSION_OFFSET`)
 */
#define MFRONT_GMP_REAL_TYPE_API_VERSION_OFFSET \
  (sizeof(mfront_gmp_real) == sizeof(double) ? 0 : 1000)
/*!
 * \brief numerical type used by the interface
 */
typedef size_t mfront_gmp_size_type;

#ifdef __cplusplus
} // end of extern "C"
#endif /* __cplusplus */
// NOLINTEND(modernize-use-using)

#endif /* LIB_MFRONT_GENERICMATERIALPROPERTY_TYPES_H */
