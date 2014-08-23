/*
 * QuantityParserContext.h
 *
 *  Created on: 22 aug 2014
 *      Author: johan
 */

#ifndef QUANTITYPARSERCONTEXT_H_
#define QUANTITYPARSERCONTEXT_H_

#include <string>
#include <map>
#include <Mod/Sketcher/App/expressions/QuantityParserExp.hpp>

#define YY_DECL QuantityParser::symbol_type yylex (QuantityParserContext& context)

namespace QuantityParserExp{

class QuantityParserContext {
public:
	QuantityParserContext ();
	virtual ~QuantityParserContext ();

	std::map<std::string, int> variables;

	Quantity result;

	//To encapsulate the coordination with the Flex scanner, it is useful to have member functions to open and close the scanning phase.

	// Handling the scanner.
	void scan_begin ();
	void scan_end ();
	bool trace_scanning;

	//Similarly for the parser itself.

	// Run the parser on file F.
	// Return 0 on success.
	int parse (const std::string& f);
	// The name of the file being parsed.
	// Used later to pass the file name to the location tracker.
	std::string file;
	// Whether parser traces should be generated.
	bool trace_parsing;

	//To demonstrate pure handling of parse errors, instead of simply dumping them on the standard error output, we will pass them to the compiler driver using the following two member functions. Finally, we close the class declaration and CPP guard.

	// Error handling.
	void error (const std::string& m);
};

} // namespace QuantityParserExp

#endif /* QUANTITYPARSERCONTEXT_H_ */
