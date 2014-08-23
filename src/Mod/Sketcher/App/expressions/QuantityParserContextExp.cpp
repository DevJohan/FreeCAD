/*
 * QuantityParserContext.cpp
 *
 *  Created on: 22 aug 2014
 *      Author: johan
 */


//The implementation of the driver is straightforward. The parse member function deserves some attention. The error functions are simple stubs, they should actually register the located error messages and set error state.

#include "QuantityParserContextExp.h"
#include "QuantityParserExp.hpp"

namespace QuantityParserExp{


QuantityParserContext::QuantityParserContext ()
: trace_scanning (false), trace_parsing (false)
{
}

QuantityParserContext::~QuantityParserContext ()
{
}

int QuantityParserContext::parse (const std::string &f)
{
	file = f;
	scan_begin ();
	QuantityParser parser (*this);
	parser.set_debug_level( trace_parsing );
	int res = parser.parse ();
	scan_end ();
	return res;
}

void
QuantityParserContext::error (const std::string& m)
{
	std::cerr << m << std::endl;
}

} // namespace QuantityParserExp
