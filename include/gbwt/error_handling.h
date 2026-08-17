#ifndef GBWT_ERROR_HANDLING_H
#define GBWT_ERROR_HANDLING_H

/*
  error_handling.h: GBWT_THROW, the macro GBWT uses to report unrecoverable
  errors.
*/

#include <stdexcept>

#define GBWT_THROW(exception) throw (exception)

#endif // GBWT_ERROR_HANDLING_H
