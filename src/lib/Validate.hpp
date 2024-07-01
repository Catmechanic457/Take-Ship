#ifndef VALIDATE_H
#define VALIDATE_H

#define NOTE "NOTE"
#define WARN "WARNING"
#define ERROR "ERROR"

#include <iostream>
#include <format>
#include <string.h>

std::string messages[] = {
    "value is out of range.",
    "value has an incorrect format.",
    "value has an invalid type.",
    "value is not normalized.",
    "value must be positive.",
    "value must be negative.",
    "value must be non-zero.",
    "array has too many elements.",
    "array has too few elements.",
    "data is missing a key.",
    "value is invalid."
};

enum error_type {
    RANGE_ERROR,
    FORMAT_ERROR,
    TYPE_ERROR,
    NOT_NORMALIZED,
    EXPECTED_POSITIVE,
    EXPECTED_NEGATIVE,
    IS_ZERO,
    TOO_MANY_ELEMENTS,
    TOO_FEW_ELEMENTS,
    MISSING_KEY,
    GENERIC_INVALID
};

/**
 * \brief Send an error message to console.
 * \param var_name Name of issue variable
 * \param type Type of error message to send
 * \param warn Type of warning to issue
*/
void send_error_message(std::string var_name, error_type type, std::string warn) {
    std::cout << '[' << warn << "]" << " Issue for input `" << var_name << "` : " << messages[type];
}


#endif