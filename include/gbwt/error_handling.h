#ifndef GBWT_ERROR_HANDLING_H
#define GBWT_ERROR_HANDLING_H

/*
  error_handling.h: GBWT_THROW, the macro GBWT uses to report unrecoverable
  errors, so that the error-reporting mechanism can be selected at build time.

  GBWT_THROW(exception) throws the given exception object, unless
  GBWT_NO_EXCEPTIONS is defined, in which case it logs exception.what() and
  calls std::abort() instead (via Abseil if GBWT_USE_ABSEIL_LOGGING is set).
*/

// GBWT_THROW always constructs its argument, even when not thrown.
#include <stdexcept>

#if defined(GBWT_NO_EXCEPTIONS)

#if defined(GBWT_USE_ABSEIL_LOGGING)
#include "absl/log/absl_log.h"
#define GBWT_THROW(exception) ABSL_LOG(FATAL) << (exception).what()
#else
#include <cstdlib>
#include <iostream>
#define GBWT_THROW(exception) \
    do { std::cerr << (exception).what() << std::endl; std::abort(); } while (0)
#endif

#else

#define GBWT_THROW(exception) throw (exception)

#endif

#endif // GBWT_ERROR_HANDLING_H
