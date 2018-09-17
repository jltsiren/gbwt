/*
  Copyright (c) 2018 Jouni Siren

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <gbwt/variants.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "GBWT unit tests";

//------------------------------------------------------------------------------

int
main(int, char**)
{
  Version::print(std::cout, tool_name);

  TestResults result;
  Verbosity::set(Verbosity::SILENT);
  
  std::cout << "variants.cpp" << std::endl;
  result += testVariants();

  std::cout << std::endl;

  std::cout << "Executed " << result.tests << " tests";
  if(result.failures > 0) { std::cout << " (" << result.failures << " failures)"; }
  std::cout << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
