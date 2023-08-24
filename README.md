# BigInteger

- Header file: `biginteger.h`
- Implementation file: `biginteger.cpp`

## Implemented
* Default constructor (initializes a number with zero).
* Constructors from numeric types.
* Explicit constructor from `std::string`.
* Copy constructor.
* Move constructor.
* Operators:
  * assignments by copying,
  * assignment by moving,
  * comparison
* Operations:
  * arithmetic operations: addition, subtraction, multiplication, division, unary minus, unary plus, increments and decrements.
  * Bitwise operations: and, or, excluding or, not, bit shifts.
* External function `std::string to_string(BigInteger const&)`.

