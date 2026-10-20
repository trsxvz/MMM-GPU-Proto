/*!
 * \file   include/MFront/GenericBehaviour/Types.h
 * \brief
 * \author Thomas Helfer
 * \date 02/07/2018
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with linking exception
 * or the CECILL-A licence. A copy of thoses licences are delivered
 * with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICBEHAVIOUR_TYPES_H
#define LIB_MFRONT_GENERICBEHAVIOUR_TYPES_H

#ifdef __cplusplus
#include <cstddef>
#else  /* __cplusplus */
#include <stddef.h>
#endif /* __cplusplus */

/*!
 * \brief numerical type used by the interface: the type given by the
 * `MFRONT_GB_REAL_TYPE` macro if it is defined (`float`, for instance), or
 * `double`. `mfront` defines this macro at the beginning of the generated
 * sources of the `generic` and `generic-parallel` interfaces when the
 * `real_type` DSL option is given.
 */
#ifdef MFRONT_GB_REAL_TYPE
typedef MFRONT_GB_REAL_TYPE mfront_gb_real;
#else  /* MFRONT_GB_REAL_TYPE */
typedef double mfront_gb_real;
#endif /* MFRONT_GB_REAL_TYPE */
/*!
 * \brief offset of the versions of the application programming interfaces
 * exported by a library whose numerical type is not `double`: such a
 * library is refused by a caller expecting `double`, even if this caller
 * does not read the size of the numerical type (symbol `<f>_real_size`)
 */
#define MFRONT_GB_REAL_TYPE_API_VERSION_OFFSET \
  (sizeof(mfront_gb_real) == sizeof(double) ? 0 : 1000)
/*!
 * \brief numerical type used by the interface
 */
typedef size_t mfront_gb_size_type;

#endif /* LIB_MFRONT_GENERICBEHAVIOUR_TYPES_H */
